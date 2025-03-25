import sys
from concurrent.futures import ThreadPoolExecutor

import polars as pl
from tqdm import tqdm
import pyarrow.parquet as pq

import tdb

def process_sample(path):
    """
    Open a sample and generate its allele count
    """
    pf = pq.ParquetFile(path)
    sample_data = pl.from_arrow(pf.read(columns=['LocusID', 'allele_number']))
    return sample_data.group_by(['LocusID', 'allele_number']).len().rename({"len": "AC"})

def merge_in(data, results):
    """
    Given a batch of samples' counts, consolidate into main data
    """
    combined = pl.concat(results, how="vertical").group_by(['LocusID', 'allele_number']).sum()
    data = data.join(combined, on=['LocusID', 'allele_number'], how="left").fill_null(0)
    data = data.with_columns([
        (pl.col("AC") + pl.col("AC_right")).alias("AC"),
    ])
    return data.drop(["AC_right"])


if __name__ == '__main__':
    # fn = "../../AoU_TRs.v0.1.tdb/"
    fn = sys.argv[1]
    batch_size = 250
    out_prefix = "result"
    min_af = 0.01
    lps_norm = 100
    names = tdb.get_tdb_filenames(fn)
    columns = ["LocusID", "allele_number"]

    # Read base allele table efficiently
    counts = pl.from_arrow(
        pq.read_table(names['allele'],
                      columns=columns)
    )

    counts = counts.with_columns(
        pl.lit(0).alias('AC'),
    )

    # Use multiple threads for reading + batch merge
    sample_paths = list(names['sample'].values())

    with ThreadPoolExecutor(max_workers=4) as executor:  # Adjust thread count
        results = []
        for sample_counts in tqdm(executor.map(process_sample, sample_paths), 
                                  total=len(sample_paths), desc="Processing"):
            results.append(sample_counts)
            
            # Merge in batches
            if len(results) >= batch_size:
                counts = merge_in(counts, results)
                results = []  # Reset batch

    # Merge remaining results
    if results:
        counts = merge_in(counts, results)

    AN = counts.group_by("LocusID").agg(pl.col('AC').sum().alias('AN'))
    counts = counts.join(AN, on="LocusID", how='left')
    counts = counts.with_columns((pl.col('AC') / pl.col('AN')).alias('AF'))

    counts.write_csv(f"{out_prefix}.allele_seq.txt", separator='\t')

    columns = ["LocusID", "allele_number", 'allele_length']
    bylen = pl.from_arrow(pq.read_table(names['allele'], columns=columns))

    bylen = bylen.join(counts, on=["LocusID", "allele_number"], how="left")  # Use "inner" if necessary
    bylen = bylen.group_by(['LocusID', 'allele_length']).agg(pl.col('AC').sum(), pl.col('AN').first())
    bylen = bylen.with_columns((pl.col('AC') / pl.col('AN')).alias('AF'))

    bylen.write_csv(f"{out_prefix}.allele_len.txt", separator='\t')

    lps = (bylen.filter(bylen['AF'] >= min_af)
        .group_by('LocusID')
        .agg((pl.col('allele_length').n_unique() /
             (pl.col('AN').first() / lps_norm)
            ).alias('LPS')
        )
    ).select(['LocusID', 'LPS'])

    lps.write_csv(f"{out_prefix}.length_polymorphism.txt", separator='\t')
