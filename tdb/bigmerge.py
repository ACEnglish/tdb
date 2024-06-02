"""
Merge multiple tdb files together.

Faster than tdb when there's more than 10 files.
"""
import os
import sys
import shutil
import logging
import argparse
import concurrent.futures

import duckdb
import truvari

import tdb

def consolidate_locus(con, db_paths, output_dir, compress=False):
    """
    con - duckdb connection
    base - the tdb we're merging to
    """
    logging.info("Consolidating locus tables")
    base = tdb.get_tdb_filenames(db_paths[0])

    query = """
    CREATE TABLE loci_lookup (
        dbname TEXT,
        on_Lhash UBIGINT,
        update_LocusID UINTEGER,
        to_LocusID UINTEGER,
        to_LocusID_new UINTEGER,
    );

    CREATE TABLE new_loci (
        LocusID UINTEGER,
        chrom TEXT,
        start UINTEGER,
        "end" UINTEGER,
    );
    """
    con.execute(query)

    for dbname in db_paths:
        names = tdb.get_tdb_filenames(dbname)
        locus_pq = names['locus']
        query = f"""
            INSERT INTO loci_lookup (dbname, on_Lhash, update_LocusID, to_LocusID)
            SELECT
                '{dbname}' AS dbname,
                hash(locus.chrom || '-' || locus.start || '-' || locus.end) AS on_Lhash,
                locus.LocusID AS update_LocusID,
                NULL AS to_LocusID,
            FROM
                read_parquet('{locus_pq}') AS locus;
        """
        con.execute(query)

    # Set what the LocusIDs should be set to_ by looking at the destination database
    bl = base['locus']
    query = f"""
        UPDATE loci_lookup
        SET to_LocusID = dest.LocusID
        FROM read_parquet('{bl}') AS dest
        WHERE loci_lookup.on_Lhash = hash(dest.chrom || '-' || dest.start || '-' || dest.end)
    """
    con.execute(query)

    seen_loci = con.execute("""
        SELECT COALESCE(MAX(to_LocusID), -1) + 1FROM loci_lookup;
    """).fetchone()[0]

    # Should check how many new loci there are
    # logging.info("no new loci to consolidate")

    # For Loci which aren't in the destination database, give them a to_LocusID_new
    query = f"""
        WITH unique_hashes AS (
            SELECT DISTINCT on_Lhash
            FROM loci_lookup
            WHERE to_LocusID IS NULL
        ),
        distinct_rows AS (
            SELECT
                on_Lhash,
                ROW_NUMBER() OVER (ORDER BY on_Lhash) + {seen_loci} AS to_LocusID_new
            FROM unique_hashes
        )
        UPDATE loci_lookup
        SET to_LocusID_new = distinct_rows.to_LocusID_new
        FROM distinct_rows
        WHERE loci_lookup.on_Lhash = distinct_rows.on_Lhash;
    """
    con.execute(query)

    con.execute("UPDATE loci_lookup SET to_LocusID_new = COALESCE(to_LocusID, to_LocusID_new);")

    # Now I have my lookup, lets see which loci I need to pull
    query = """
    CREATE TABLE loci_pull AS
    SELECT DISTINCT ON (on_Lhash)
        dbname,
        on_Lhash,
    FROM loci_lookup
    WHERE to_LocusID IS NULL;

    SELECT DISTINCT ON (dbname) dbname FROM loci_pull
    """
    to_pull = con.execute(query).fetchall()

    logging.info("Merging locus")
    for dbname, in to_pull:
        names = tdb.get_tdb_filenames(dbname)
        m_locus = names['locus']
        query = f"""
            INSERT INTO new_loci (LocusID, chrom, start, "end")
            SELECT
                loci_lookup.to_LocusID_new AS LocusID,
                original.chrom,
                original.start,
                original."end"
            FROM
                loci_pull
            JOIN
                loci_lookup
                ON loci_pull.dbname = loci_lookup.dbname AND loci_pull.on_Lhash = loci_lookup.on_Lhash
            JOIN
                read_parquet('{m_locus}') AS original
                ON loci_lookup.on_Lhash = hash(original.chrom || '-' || original.start || '-' || original."end")
            WHERE
                loci_pull.dbname = '{dbname}';
        """
        con.execute(query)

    olocus = os.path.join(output_dir, "locus.pq")
    comp = ""
    do_order = ""
    if compress:
        comp = ", COMPRESSION GZIP"
        do_order = 'ORDER BY chrom, start, "end"'

    query = f"""
    COPY (
        SELECT *
        FROM read_parquet('{base["locus"]}')
        UNION ALL
        SELECT *
        FROM new_loci
        {do_order}
    ) TO '{olocus}' (FORMAT PARQUET{comp})
    """
    con.execute(query)


def allele_pusher(con, dbname):
    """
    Read input tdb alleles for consolidation
    """
    local_con = con.cursor()

    names = tdb.get_tdb_filenames(dbname)
    allele_pq = names['allele']
    query = f"""
    INSERT INTO allele_lookup (dbname, update_LocusID, to_LocusID, update_allele_number, to_allele_number, to_allele_number_new, on_Ahash)
    SELECT
        '{dbname}' AS dbname,
        allele.LocusID AS update_LocusID,
        loci_lookup.to_LocusID_new as to_LocusID,
        allele.allele_number as update_allele_number,
        NULL AS to_allele_number,
        NULL AS to_allele_number_new,
        hash(CAST(allele.sequence AS TEXT)) AS on_Ahash
    FROM
        read_parquet('{allele_pq}') AS allele
    JOIN
        loci_lookup
    ON
        loci_lookup.dbname = '{dbname}'
        AND loci_lookup.update_LocusID = allele.LocusID
    """
    local_con.execute(query).fetchall()

def allele_puller(con, dbname, num_loci):
    logging.debug("pulling %d alleles from %s", num_loci, dbname)
    local_con = con.cursor()
    names = tdb.get_tdb_filenames(dbname)
    m_allele = names['allele']
    query = f"""
        INSERT INTO new_alleles (LocusID, allele_number, allele_length, sequence)
        SELECT
            allele_pull.to_LocusID AS LocusID,
            allele_pull.to_allele_number_new AS allele_number,
            original.allele_length,
            original.sequence
        FROM
            allele_pull
        JOIN
            read_parquet('{m_allele}') AS original
        ON allele_pull.update_LocusID = original.LocusID
            AND allele_pull.update_allele_number = original.allele_number
        WHERE
            allele_pull.dbname = '{dbname}';
    """
    local_con.execute(query).fetchall()
    logging.debug("pulled alleles from %s", dbname)


def consolidate_allele(con, db_paths, output_dir, compress=False, threads=1):
    """
    Consolidates alleles using db_paths[0] as the baseline
    """
    logging.info("Consolidating allele tables")
    base = tdb.get_tdb_filenames(db_paths[0])
    query = """
    CREATE TABLE allele_lookup (
        dbname TEXT,
        update_LocusID UINTEGER,
        to_LocusID UINTEGER,
        update_allele_number USMALLINT,
        to_allele_number USMALLINT,
        to_allele_number_new USMALLINT,
        on_Ahash UBIGINT,
    );

    CREATE TABLE new_alleles (
        LocusID UINTEGER,
        allele_number USMALLINT,
        allele_length USMALLINT,
        sequence BLOB,
    );
    """
    con.execute(query)

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(allele_pusher, con, dbname) for dbname in db_paths]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()  # This will raise an exception if the task failed
            except Exception as e:
                logging.error(f"An error occurred: {e}")

    # Set what the allele_number should be set to_ by looking at the destination database
    ba = base['allele']
    query = f"""
    UPDATE allele_lookup
    SET to_allele_number = dest.allele_number
    FROM read_parquet('{ba}') AS dest
    WHERE
        allele_lookup.to_LocusID = dest.LocusID
        AND allele_lookup.on_Ahash = hash(CAST(dest.sequence AS TEXT))
    """
    con.execute(query)

    con.execute("""
    CREATE TEMP TABLE temp_distinct_rows AS
    WITH unique_hashes AS (
        SELECT DISTINCT
            on_Ahash,
            to_LocusID,
            to_allele_number
        FROM allele_lookup
    )
    SELECT
        on_Ahash,
        to_LocusID,
        ROW_NUMBER() OVER (PARTITION BY to_LocusID ORDER BY to_allele_number) - 1 AS to_allele_number_new
    FROM unique_hashes;
    """)

    logging.debug("updating allele_number")
    query = """
    UPDATE allele_lookup
    SET to_allele_number_new = COALESCE(allele_lookup.to_allele_number, temp_distinct_rows.to_allele_number_new)
    FROM temp_distinct_rows
    WHERE
        allele_lookup.to_LocusID = temp_distinct_rows.to_LocusID
        AND allele_lookup.on_Ahash = temp_distinct_rows.on_Ahash;
    """
    con.execute(query)

    # Now I need to figure out which alleles are new
    logging.debug("figuring out alleles to pull")
    query = """
    CREATE TABLE allele_pull AS
    SELECT DISTINCT ON (to_LocusID, to_allele_number_new)
        dbname,
        to_LocusID,
        to_allele_number_new,
        update_LocusID,
        update_allele_number,
    FROM allele_lookup
    WHERE to_allele_number IS NULL;

    SELECT dbname, COUNT(*) as occurrences
    FROM allele_pull
    GROUP BY dbname
    ORDER BY occurrences DESC;
    """
    to_pull = con.execute(query).fetchall()

    logging.info("Merging allele")
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(allele_puller, con, dbname, num_loci) for dbname, num_loci in to_pull]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()  # This will raise an exception if the task failed
            except Exception as e:
                logging.error(f"An error occurred: {e}")

    comp = ""
    do_order = ""
    if compress:
        comp = ", COMPRESSION GZIP"
        do_order = "ORDER BY LocusID, allele_number"

    alocus = os.path.join(output_dir, "allele.pq")
    query = f"""
    COPY (
        SELECT *
        FROM read_parquet('{base["allele"]}')
        UNION ALL
        SELECT *
        FROM new_alleles
        {do_order}
    ) TO '{alocus}' (FORMAT PARQUET{comp})
    """
    con.execute(query)

def sample_puller(con, dbname, output_dir, compress):
    """
    Translates samples to new ids and moves
    """
    logging.debug("pulling samples from %s", dbname)
    comp = ""
    do_order = ""
    if compress:
        comp = ", COMPRESSION GZIP"
        do_order = "ORDER BY LocusID, allele_number"

    files = tdb.get_tdb_filenames(dbname)
    for _, sample_pq in files['sample'].items():
        out_name = os.path.join(output_dir, os.path.basename(sample_pq))
        query = f"""
            COPY (
                SELECT
                    allele_lookup.to_LocusID as LocusID,
                    allele_lookup.to_allele_number_new as allele_number,
                    sample.spanning_reads,
                    sample.length_range_lower,
                    sample.length_range_upper,
                    sample.average_methylation,
                FROM read_parquet('{sample_pq}') as sample
                JOIN allele_lookup
                ON
                    allele_lookup.dbname == '{dbname}'
                    AND allele_lookup.update_LocusID == sample.LocusID
                    AND allele_lookup.update_allele_number == sample.allele_number
                {do_order}
            ) TO '{out_name}' (FORMAT PARQUET{comp})
        """
        con.execute(query)
    logging.debug("pulled %d samples from %s", len(files['sample']), dbname)


def consolidate_sample(con, db_names, output_dir, compress=False, threads=1):
    """
    Translate each sample table to the new database's keys
    """
    logging.info("Updating samples")
    # would need to copy from the base
    base = tdb.get_tdb_filenames(db_names[0])
    for sample_pq in base['sample'].values():
        out_name = os.path.join(output_dir, os.path.basename(sample_pq))
        shutil.copy(sample_pq, out_name)

    # And then update the rest
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(sample_puller, con, dbname, output_dir, compress) for dbname in db_names[1:]]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()  # This will raise an exception if the task failed
            except Exception as e:
                logging.error(f"An error occurred: {e}")


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


def bigmerge_main(args):
    """
    bigmerge main entrypoint
    """
    parser = argparse.ArgumentParser(prog="tdb bigmerge", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT",
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
        logging.error("Cannot merge database. Exiting")
        sys.exit(1)

    os.mkdir(args.output)

    temp_db = truvari.make_temp_filename(suffix=".duckdb")
    con = duckdb.connect(temp_db)
    con.execute("SET default_null_order ='NULLS LAST';")
    con.execute(f"SET threads = {args.threads};")
    if args.mem:
        con.execute(f"SET memory_limit = '{args.mem}GB';")

    consolidate_locus(con, args.inputs, args.output, args.no_compress)
    consolidate_allele(con, args.inputs, args.output, args.no_compress, args.threads)
    consolidate_sample(con, args.inputs, args.output, args.no_compress, args.threads)

    con.close()
