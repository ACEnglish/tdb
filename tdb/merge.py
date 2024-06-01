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

GLOBAL_DUCK_SET = ["SET default_null_order ='NULLS LAST';"]


def setup_duck(con):
    """
    Set duckdb threads/memory
    """
    for i in GLOBAL_DUCK_SET:
        con.execute(i)


def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False

    if os.path.exists(args.output):
        logging.error(f"Output {args.output} already exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"Output {args.output} must end with `.tdb`")
        check_fail = True
    seen_samples = {}
    for i in args.inputs:
        if not os.path.exists(i):
            logging.error(f"Input {i} does not exist")
            check_fail = True
        if not i.rstrip('/').endswith(".tdb"):
            logging.error(f"Unrecognized file extension on {i} expected .tdb")
            check_fail = True
        else:  # can only check sample of valid file names
            for s in tdb.get_tdb_samplenames(i):
                if s in seen_samples:
                    logging.error(
                        f"Input {i} has redundant sample {s} with {seen_samples[s]}")
                    check_fail = True
                seen_samples[s] = i
    return check_fail


def join_loci_tables(original_loci, second_loci, compress):
    """
    Creates a union of two loci tables

    Returns a temporary filename holding the joined LocusIDs as
    update_LocusID to_LocusID
    as well as a path the the updated locus table which should be moved
    """
    con = duckdb.connect()
    setup_duck(con)

    loci_lookup_parquet_path = truvari.make_temp_filename(suffix=".pq")
    up_locus = truvari.make_temp_filename(suffix=".pq")

    comp = ", COMPRESSION GZIP" if compress else ""
    do_order = "ORDER BY update_LocusID, to_LocusID, chrom, start" if compress else ""

    query = f"""
    COPY (
        SELECT
            second.locusid AS update_LocusID,
            original.locusid AS to_LocusID,
            COALESCE(original.chrom, second.chrom) AS chrom,
            COALESCE(original.start, second.start) AS start,
            COALESCE(original.end, second.end) AS end
        FROM
            read_parquet('{second_loci}') AS second
        FULL JOIN
            read_parquet('{original_loci}') AS original
        ON
            second.chrom = original.chrom
            AND second.start = original.start
            AND second.end = original.end
        {do_order}
    ) TO '{loci_lookup_parquet_path}' (FORMAT PARQUET{comp})
    """
    con.execute(query)
    con.close()

    union = pd.read_parquet(loci_lookup_parquet_path)
    mask = union['to_LocusID'].isna()
    start_id = union['to_LocusID'].max() + 1
    union.loc[mask, 'to_LocusID'] = np.arange(start_id, start_id + mask.sum())
    union['to_LocusID'] = union['to_LocusID'].astype(np.uint32)
    new_loci = mask.sum()
    if new_loci:
        logging.info("new loci: %d", new_loci)

    # Write the new locus table
    to_write = union.rename(columns={"to_LocusID": "LocusID"})[
        ["LocusID", "chrom", "start", "end"]]
    to_write.to_parquet(up_locus, index=False,
                        compression='gzip' if compress else None)

    # And write the updated lookup
    union[["update_LocusID", "to_LocusID"]].to_parquet(
        loci_lookup_parquet_path, index=False)

    return loci_lookup_parquet_path, up_locus


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
    ) TO '{second_updated_locusid}' (FORMAT PARQUET)
    """
    con.execute(create_updated_allele)
    con.close()

    return second_updated_locusid


def create_allele_lookup(original_allele, second_allele):
    """
    create an allele lookup between two db's allele tables


    Note that the second allele table must already have its loci updated

    Returns the path of a new alleles table and the full allele lookup
    """
    con = duckdb.connect()
    setup_duck(con)

    new_alleles_path = truvari.make_temp_filename(suffix=".pq")
    allele_lookup_path = truvari.make_temp_filename(suffix=".pq")

    query = f"""
    CREATE TABLE lookup AS SELECT
        COALESCE(original.LocusID, second.LocusID) AS LocusID,
        second.allele_number AS update_allele_number,
        original.allele_number AS to_allele_number,
        CAST(
            ROW_NUMBER() OVER (PARTITION BY COALESCE(original.LocusID, second.LocusID) ORDER BY original.allele_number) - 1 
        AS USMALLINT) AS to_allele_number_new
    FROM
        read_parquet('{original_allele}') AS original
    FULL JOIN
        read_parquet('{second_allele}') AS second
    ON
        second.LocusID == original.LocusID
        AND second.allele_length == original.allele_length
        AND second.sequence == original.sequence
    ORDER BY LocusID, update_allele_number, to_allele_number;

    COPY (
        SELECT
            LocusID,
            update_allele_number,
            to_allele_number_new
        FROM
            lookup
        WHERE
            update_allele_number IS NOT NULL
    ) TO '{allele_lookup_path}' (FORMAT PARQUET);

     COPY (
        SELECT
            LocusID,
            update_allele_number, 
            to_allele_number_new,
        FROM 
            lookup
        WHERE
            to_allele_number IS NULL
    ) TO '{new_alleles_path}' (FORMAT PARQUET)
    """
    con.execute(query)

    new_allele_count = con.execute("""
        SELECT COUNT(*) AS row_count
        FROM lookup
        WHERE to_allele_number IS NULL;
        """).fetchone()[0]
    logging.info("new alleles: %d", new_allele_count)

    con.close()

    return new_alleles_path, allele_lookup_path


def consolidate_alleles(original_allele, second_allele, new_alleles, compress):
    """
    This consolidates the original alleles with the second alleles.
    Only the subset of new alleles in second alleles are pulled.
    The new alleles are renumbered according to `update_allele_numbers`

    Returns a path to a new allele table which should be moved when ready
    """
    spillover = truvari.make_temp_filename(suffix=".db")
    con = duckdb.connect(spillover)
    setup_duck(con)

    up_allele = truvari.make_temp_filename(suffix=".pq")

    comp = ", COMPRESSION GZIP" if compress else ""
    do_order = "ORDER BY LocusID, allele_number, allele_length, sequence" if compress else ""

    query = f"""
    COPY (
        WITH subset_alleles AS (
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
                second.LocusID = new_alleles.LocusID
                AND second.allele_number = new_alleles.update_allele_number
        )
        SELECT *
        FROM read_parquet('{original_allele}')
        UNION ALL
        SELECT *
        FROM subset_alleles
        {do_order}
    ) TO '{up_allele}' (FORMAT PARQUET{comp});
    """
    con.execute(query)
    con.close()

    shutil.os.remove(spillover)
    return up_allele


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


def update_sample_table(second_sample, sample_lookup, compress):
    """
    Update the LocusID and allele_number of a second_sample to the consolidated ids.

    Writes directly to the destination output_path
    """
    con = duckdb.connect()
    setup_duck(con)
    
    output_path = truvari.make_temp_filename(suffix=".pq")
    comp = ", COMPRESSION GZIP" if compress else ""
    do_order = "ORDER BY LocusID, allele_number" if compress else ""

    query = f"""
    COPY(
        SELECT
            lookup.to_LocusID AS LocusID,
            lookup.to_allele_number_new AS allele_number,
            sample.spanning_reads AS spanning_reads,
            sample.length_range_lower AS length_range_lower,
            sample.length_range_upper AS length_range_upper,
            sample.average_methylation AS average_methylation
        FROM
            read_parquet('{second_sample}') AS sample
        LEFT JOIN
            read_parquet('{sample_lookup}') AS lookup
        ON
            lookup.update_LocusID = sample.LocusID
            AND lookup.update_allele_number == sample.allele_number
        {do_order}
    ) TO '{output_path}' (FORMAT PARQUET{comp})
    """
    con.execute(query)
    con.close()
    return output_path


def tdb_consolidate(tdb_1, tdb_2, allele_gz=True, samp_gz=True):
    """
    Consolidates two tdbs into a destination
    Assumes tdb_1 is the base and its ids won't change.
    tdb_1's samples aren't touched. Those should have already been handled.

    Automatically removes the temporary flies it creates
    The tdbs should be get_tdb_filenames or whatever it's called
    """
    loci_lookup, up_locus = join_loci_tables(tdb_1['locus'], tdb_2['locus'], allele_gz)
    second_allele_locus_updated = update_allele_locusid(tdb_2['allele'],
                                                        loci_lookup)
    new_alleles, allele_lookup = create_allele_lookup(tdb_1['allele'],
                                                      second_allele_locus_updated)
    logging.debug("Consolidating alleles")
    up_allele = consolidate_alleles(tdb_1['allele'], second_allele_locus_updated,
                        new_alleles, allele_gz)
    logging.debug("Creating sample lookup")
    sample_lookup = create_sample_lookup(loci_lookup, allele_lookup)

    output_dir = os.path.dirname(tdb_1['locus'])
    logging.debug("Updating samples")
    up_samples = {}
    for name, second_sample in tdb_2['sample'].items():
        up_samples[name] = update_sample_table(second_sample, sample_lookup, samp_gz)

    logging.debug("putting into output")
    shutil.move(up_locus, tdb_1['locus'])
    shutil.move(up_allele, tdb_1['allele'])
    for name, tmp_path in up_samples.items():
        output_name = os.path.join(output_dir, f"sample.{name}.pq")
        shutil.move(tmp_path, output_name)

    # Clean up after yourself
    # You know, if you give them static names you won't need to clean but once...
    shutil.os.remove(loci_lookup)
    shutil.os.remove(second_allele_locus_updated)
    shutil.os.remove(new_alleles)
    shutil.os.remove(allele_lookup)
    shutil.os.remove(sample_lookup)


def merge_main(args):
    """
    Merge multiple tdbs
    """
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
                        help="Don't sort/compress output (faster merge, bigger output)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument("inputs", metavar="IN", nargs="+",
                        help="tdb files")
    args = parser.parse_args(args)

    truvari.setup_logging(args.debug)

    if check_args(args):
        logging.error("Cannot create database. Exiting")
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
        pos += 2
        logging.info("Consolidating %s (%d/%d)", i, pos, num_tdbs)
        update_tdb = tdb.get_tdb_filenames(i)
        acomp = args.no_compress and pos == num_tdbs
        tdb_consolidate(dest_tdb, update_tdb,
                        allele_gz=acomp,
                        samp_gz=args.no_compress)
    logging.info("Finished")
