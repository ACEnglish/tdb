"""
Merge putatively redundant alleles per-locus
"""
import sys
import tdb
import truvari
import numpy as np


def build_matrix(alleles):
    """
    Compare all-against-all to build a boolean adjacency matrix
    """
    n_entries = len(alleles)
    match_matrix = np.zeros((n_entries, n_entries), dtype=bool)
    for i in range(n_entries - 1):
        for j in range(i + 1, n_entries):
            szsim, _ = truvari.sizesim(len(alleles[i]), len(alleles[j]))
            state = True
            if szsim < THRESHOLD:
                state = False
            else:
                sqsim = truvari.seqsim(alleles[i], alleles[j])
                if sqsim < THRESHOLD:
                    state = False
            match_matrix[i, j] = state
            match_matrix[j, i] = state
    return match_matrix


def find_matching_sets(matrix):
    """
    Creates a lookup of which alleles match
    returns the new allele numbers lookup as dict and a list of original
    allele numbers to keep
    """
    n = len(matrix)
    visited = [False] * n
    matched_sets = []

    def dfs(item, current_set):
        """
        Depth first search to find chain of matches
        """
        visited[item] = True
        current_set.append(item + 1)  # alt alleles start at number 1
        for other in range(n):
            if matrix[item][other] and not visited[other]:
                dfs(other, current_set)

    for i in range(n):
        if not visited[i]:
            current_set = []
            dfs(i, current_set)
            matched_sets.append(current_set)

    # This just keeps the first allele
    # A better strategy would be to keep the most frequently observed allele
    to_keep = [idx[0] for idx in matched_sets]
    # Create a lookup of old allele number to the new allele number
    to_rename = {old_num: new_num + 1
                 for new_num, entry_set in enumerate(matched_sets)
                 for old_num in entry_set}

    return to_rename, to_keep


if __name__ == '__main__':
    THRESHOLD = 0.98

    d = tdb.load_tdb(sys.argv[1])
    # This currently just collects stats on how many loci/alleles would benefit from merging
    # Needs to be updated to create a new tdb
    talleles = 0
    ralleles = 0
    tloci = 0
    rloci = 0
    for grp, alleles in d['allele'].groupby(['LocusID']):
        matrix = build_matrix(list(alleles['sequence']))
        to_rename, to_keep = find_matching_sets(matrix)
        a1 = len(alleles)
        a2 = len(to_keep)
        talleles += a1
        ralleles += a2
        tloci += 1
        if a1 != a2:
            rloci += 1

    print('loci', tloci, rloci, rloci/tloci)
    print('alleles', talleles, ralleles, ralleles/talleles)
