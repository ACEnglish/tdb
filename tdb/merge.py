"""
Merge multiple tdb files together
"""
import os
import sys
import shutil
import logging
import argparse
import duckdb
import truvari
import numpy as np
import pandas as pd

import tdb
from tdb.create import get_samples

GLOBAL_DUCK_SET=[]
def setup_duck(con):
    global GLOBAL_DUCK_SET
    for i in GLOBAL_DUCK_SET:
        con.execute(i)

def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False

    if os.path.exists(args.output):
        logging.error(f"output {args.output} already exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"output {args.output} must end with `.tdb`")
        check_fail = True
    seen_samples = {}
    for i in args.inputs:
        if not os.path.exists(i):
            logging.error(f"input {i} does not exist")
            check_fail = True
        if not i.rstrip('/').endswith(".tdb"):
            logging.error(f"unrecognized file extension on {i}")
            logging.error("expected .tdb")
            check_fail = True
        else: # can only check sample of valid file names
            for s in get_samples(i):
                if s in seen_samples:
                    logging.error(f"input {i} has redundant sample with {seen_samples[s]}")
                    check_fail = True
                seen_samples[s] = i
    return check_fail


def join_loci_tables(original_loci, second_loci):
    """
    Creates a union of two loci tables

    Returns a temporary filename holding the joined LocusIDs as
    update_LocusID to_LocusID
    """
    con = duckdb.connect()
    setup_duck(con)

    loci_lookup_parquet_path = truvari.make_temp_filename(suffix=".pq")

    query = f"""
    COPY (
        SELECT
            second.locusid AS update_LocusID,
            original.locusid AS to_LocusID
        FROM
            read_parquet('{second_loci}') AS second
        FULL JOIN
            read_parquet('{original_loci}') AS original
        ON
            second.chrom = original.chrom
            AND second.start = original.start
            AND second.end = original.end
        ORDER BY update_LocusID, to_LocusID
    ) TO '{loci_lookup_parquet_path}' (FORMAT PARQUET)
    """

    con.execute(query)
    con.close()

    # WARNING!! I'm not doing the disjoint solving, yet
    return loci_lookup_parquet_path

def update_allele_locusid(second_allele, loci_lookup):
    """
    Makes a temporary allele table with updated LocusID
    
    returns the path of the temporary allele table
    """
    con = duckdb.connect()
    setup_duck(con)

    second_updated_locusid = truvari.make_temp_filename(suffix=".pq")

    create_updated_allele = f"""
    COPY (
        SELECT
            lookup.to_LocusID as LocusID,
            second.allele_number as allele_number,
            second.allele_length as allele_length,
            second.sequence as sequence,
        FROM
            read_parquet('{second_allele}') AS second
        LEFT JOIN
            read_parquet('{loci_lookup}') AS lookup
        ON
            second.LocusID == lookup.update_LocusID
        ORDER BY LocusID, allele_number
    ) TO '{second_updated_locusid}' (FORMAT PARQUET)
    """

    con.execute(create_updated_allele)
    con.close()
    
    return second_updated_locusid

def create_allele_lookup(original_allele, second_allele):
    """
    create an allele lookup between two db's allele tables
    
    Note that the second allele table must already have its loci updated
    
    Returns the path to the allele lookup
    """
    con = duckdb.connect()
    setup_duck(con)

    partial_lookup = truvari.make_temp_filename(suffix=".pq")

    query = f"""
    COPY (
        SELECT
            original.LocusID as LocusID_orig,
            second.LocusID as LocusID_second,
            second.allele_number as update_allele_number,
            original.allele_number as to_allele_number,
        FROM
            read_parquet('{original_allele}') AS original
        FULL JOIN
            read_parquet('{second_allele}') AS second
        ON
            second.LocusID == original.LocusID
            AND second.allele_length == original.allele_length
            AND second.sequence == original.sequence
        ORDER BY LocusID_orig, LocusID_second, update_allele_number, to_allele_number
    ) TO '{partial_lookup}' (FORMAT PARQUET)
    """
    
    con.execute(query)
    con.close()
    
    return partial_lookup

def update_allele_numbers(partial_lookup):
    """
    Given the create_allele_lookup table, renumber the second allele tables allele_numbers

    Writes a temporary file which holds all the new alleles in second relative to first
    And a temporary file with allele lookup information which will be used for updating the
    second database's sample tables
    """
    new_alleles_path = truvari.make_temp_filename(suffix=".pq")
    allele_lookup_path = truvari.make_temp_filename(suffix=".pq")

    data = pd.read_parquet(partial_lookup)
    data['LocusID'] = data['LocusID_orig'].combine_first(data['LocusID_second'])
    data.drop(columns=["LocusID_orig", "LocusID_second"], inplace=True)
    data.sort_values(["LocusID", "to_allele_number", "update_allele_number"], inplace=True)
    data['to_allele_number_new'] = data.groupby(['LocusID']).cumcount()

    lktypes = {"LocusID": np.uint32, "update_allele_number": np.uint16, "to_allele_number_new": np.uint16}
    new_alleles = data[data['to_allele_number'].isna()]
    new_alleles = new_alleles.drop(columns=["to_allele_number"]).astype(lktypes)
    new_alleles.to_parquet(new_alleles_path, index=False)
    
    logging.info("new alleles: %d", len(new_alleles))
    # I need to translate all the alleles... in the sample table. This is only for sample table
    allele_lookup = data[["LocusID", "update_allele_number", "to_allele_number_new"]].dropna().astype(lktypes)
    allele_lookup.to_parquet(allele_lookup_path, index=False)

    return new_alleles_path, allele_lookup_path


def consolidate_alleles(original_allele, second_allele, new_alleles, compress):
    """
    This consolidates the original alleles with the second alleles.
    Only the subset of new alleles in second alleles are pulled.
    The new alleles are renumbered according to `update_allele_numbers`

    This overwrites original_allele
    """
    con = duckdb.connect()
    setup_duck(con)

    tmp = truvari.make_temp_filename(suffix=".pq")

    subset_query = f"""
    COPY (
        SELECT
            new_alleles.LocusID AS LocusID,
            new_alleles.to_allele_number_new AS allele_number,
            second.allele_length AS allele_length,
            second.sequence AS sequence
        FROM
            read_parquet('{new_alleles}') AS new_alleles
        LEFT JOIN
            read_parquet('{second_allele}') AS second
        ON
            second.LocusID == new_alleles.LocusID
            AND second.allele_number == new_alleles.update_allele_number
        ORDER BY LocusID, allele_number, allele_length, sequence
    ) TO '{tmp}' (FORMAT PARQUET)
    """
    con.execute(subset_query)

    tmp2 = truvari.make_temp_filename(suffix=".pq")
    # And concatenate (UNION)
    comp = ", COMPRESSION GZIP" if compress else ""
    concatenate_query = f"""
    COPY (
        SELECT 
            CAST(LocusID AS UINTEGER) AS LocusID,
            CAST(allele_number AS USMALLINT) AS allele_number,
            CAST(allele_length AS USMALLINT) AS allele_length,
            CAST(sequence AS BLOB) AS sequence
        FROM read_parquet('{original_allele}')
        UNION ALL
        SELECT
            CAST(LocusID AS UINTEGER) AS LocusID,
            CAST(allele_number AS USMALLINT) AS allele_number,
            CAST(allele_length AS USMALLINT) AS allele_length,
            CAST(sequence AS BLOB) AS sequence
        FROM read_parquet('{tmp}')
        ORDER BY LocusID, allele_number, allele_length, sequence
    ) TO '{tmp2}' (FORMAT PARQUET{comp})
    """

    con.execute(concatenate_query)
    con.close()

    shutil.move(tmp2, original_allele)
    shutil.os.remove(tmp)

def create_sample_lookup(loci_lookup, allele_lookup):
    """
    Join the loci lookup and allele lookup so that we can translate the samples' LocusID and allele_number columns
    
    Returns the temporary sample_lookup file
    """
    con = duckdb.connect()
    setup_duck(con)

    sample_lookup = truvari.make_temp_filename(suffix=".pq")

    query = f"""
    COPY (
        SELECT
            *
        FROM
            read_parquet('{allele_lookup}') AS allele_lookup
        LEFT JOIN
            read_parquet('{loci_lookup}') as loci_lookup
        ON
            loci_lookup.to_LocusID == allele_lookup.LocusID
    ) TO '{sample_lookup}' (FORMAT PARQUET)
    """
    con.execute(query)
    con.close()
    return sample_lookup

def update_sample_table(second_sample, sample_lookup, output_path, compress):
    """
    Update the LocusID and allele_number of a second_sample to the consolidated ids.

    Writes directly to the destination output_path
    """
    con = duckdb.connect()
    setup_duck(con)

    comp = ", COMPRESSION GZIP" if compress else ""
    query = f"""
    COPY(
        SELECT
            CAST(lookup.to_LocusID AS UINTEGER) AS LocusID,
            CAST(lookup.to_allele_number_new AS USMALLINT) AS allele_number,
            CAST(sample.spanning_reads AS USMALLINT) AS spanning_reads,
            CAST(sample.length_range_lower AS USMALLINT)  AS length_range_lower,
            CAST(sample.length_range_upper AS USMALLINT) AS length_range_upper,
            CAST(sample.average_methylation AS FLOAT4) AS average_methylation
        FROM
            read_parquet('{second_sample}') AS sample
        LEFT JOIN
            read_parquet('{sample_lookup}') AS lookup
        ON
            lookup.update_LocusID = sample.LocusID
            AND lookup.update_allele_number == sample.allele_number
        ORDER BY LocusID, allele_number
    ) TO '{output_path}' (FORMAT PARQUET{comp})
    """
    con.execute(query)
    con.close()

def tdb_consolidate(tdb_1, tdb_2, allele_gz=True, samp_gz=True):
    """
    Consolidates two tdbs into a destination
    Assumes tdb_1 is the base and its ids won't change.
    tdb_1's samples aren't touched. Those should have already been handled.

    Automatically removes the temporary flies it creates
    The tdbs should be get_tdb_filenames or whatever it's called
    """
    loci_lookup = join_loci_tables(tdb_1['locus'], tdb_2['locus'])
    # I'm still not writing joined locus.
    second_allele_locus_updated = update_allele_locusid(tdb_2['allele'], loci_lookup)
    partial_lookup = create_allele_lookup(tdb_1['allele'], second_allele_locus_updated)
    new_alleles, allele_lookup = update_allele_numbers(partial_lookup)
    consolidate_alleles(tdb_1['allele'], second_allele_locus_updated, new_alleles, allele_gz)
    sample_lookup = create_sample_lookup(loci_lookup, allele_lookup)

    output_dir = os.path.dirname(tdb_1['locus'])
    for name, second_sample in tdb_2['sample'].items():
        output_name = os.path.join(output_dir, f"sample.{name}.pq")
        update_sample_table(second_sample, sample_lookup, output_name, samp_gz)

    # Clean up after yourself
    # You know, if you give them static names you won't need to clean but once...
    shutil.os.remove(loci_lookup)
    shutil.os.remove(second_allele_locus_updated)
    shutil.os.remove(partial_lookup)
    shutil.os.remove(new_alleles)
    shutil.os.remove(allele_lookup)
    shutil.os.remove(sample_lookup)

def merge_main(args):
    global GLOBAL_DUCK_SET
    parser = argparse.ArgumentParser(prog="tdb merge", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--into", merge into the first tdb listed instead of copying
    # it into the output. This will replace append
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output tdb directory")
    parser.add_argument("--mem", default=None,
                        help="Maximum memory in GB (off)")
    parser.add_argument("--threads", default=1, type=int,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--no-compress", action="store_false",
                        help="Skip compression (faster merge, bigger output)")
    parser.add_argument("inputs", metavar="IN", nargs="+",
                        help="tdb files")
    args = parser.parse_args(args)

    truvari.setup_logging()

    if check_args(args):
        logging.error("cannot create database. exiting")
        sys.exit(1)
    
    GLOBAL_DUCK_SET.append(f"SET threads = {args.threads};")
    if args.mem:
        GLOBAL_DUCK_SET.append(f"SET memory_limit = '{args.mem}GB';")

    num_tdbs = len(args.inputs)
    first_tdb_name = args.inputs.pop(0)
    logging.info("Consolidating %s (1/%d)", first_tdb_name, num_tdbs)
    shutil.copytree(first_tdb_name, args.output)
    dest_tdb = tdb.get_tdb_filenames(args.output)
    
    for pos, i in enumerate(args.inputs):
        pos += 1 
        logging.info("Consolidating %s (%d/%d)", i, pos, num_tdbs)
        update_tdb = tdb.get_tdb_filenames(i)
        acomp = args.no_compress and pos == num_tdbs
        tdb_consolidate(dest_tdb, update_tdb,
                         allele_gz=acomp,
                         samp_gz=args.no_compress)
    logging.info("Finished")
