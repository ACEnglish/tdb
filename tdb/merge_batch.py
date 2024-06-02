"""
Merge multiple tdb files together.

Faster than tdb when there's more than 10 files.
"""
import os
import tdb
import duckdb
import glob
import shutil
import truvari
import logging

def consolidate_locus(con, db_paths, output_dir):
    """
    con - duckdb connection
    base - the tdb we're merging to
    """
    loging.info("Consolidating locus tables")
    base = tdb.get_tdb_filenames(db_paths[0])

    #TODO: Back to worrying about types
    query = """
    CREATE TABLE loci_lookup (
        dbname TEXT,
        on_Lhash TEXT,
        update_LocusID INTEGER,
        to_LocusID INTEGER,
        to_LocusID_new INTEGER,
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
                md5(locus.chrom || '-' || locus.start || '-' || locus.end) AS on_Lhash,
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
        WHERE loci_lookup.on_Lhash = md5(dest.chrom || '-' || dest.start || '-' || dest.end)
    """
    con.execute(query)

    seen_loci, new_loci = con.execute(f"""
        SELECT COALESCE(MAX(to_LocusID), -1) + 1, COUNT(*) FROM loci_lookup;
    """).fetchone()[0]
    
    if new_loci == 0:
        logging.info("no loci to consolidate")

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
    query = f"""
    CREATE TABLE loci_pull AS
    SELECT DISTINCT ON (on_Lhash)
        dbname,
        on_Lhash,
    FROM loci_lookup
    WHERE to_LocusID IS NULL;

    SELECT DISTINCT ON (dbname) dbname FROM loci_pull
    """
    to_pull = con.execute(query).fetchall()

    # Make a temporary file holding the new loci entries
    query = """
    CREATE TABLE new_loci (
        LocusID INTEGER,
        chrom TEXT,
        start INTEGER,
        "end" INTEGER,
    );
    """
    con.execute(query)

    logging.info("Updating locus")
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
                ON loci_lookup.on_Lhash = md5(original.chrom || '-' || original.start || '-' || original."end")
            WHERE
                loci_pull.dbname = '{dbname}';
        """
        con.execute(query)

    # TODO: needs a do_order/compress
    olocus = os.path.join(output_dir, "locus.pq")
    query = f"""
    COPY (
        SELECT * 
        FROM read_parquet('{base["locus"]}')
        UNION ALL
        SELECT *
        FROM new_loci
    ) TO '{olocus}' (FORMAT PARQUET)
    """
    con.execute(query)

def consolidate_allele(con, db_paths, output_dir):
    """
    Consolidates alleles using db_paths[0] as the baseline
    """
    logging.info("Consolidating allele")
    base = tdb.get_tdb_filenames(db_paths[0])
    query = """
    CREATE TABLE allele_lookup (
        dbname TEXT,
        update_LocusID INTEGER,
        to_LocusID INTEGER,
        update_allele_number INTEGER,
        to_allele_number INTEGER,
        to_allele_number_new INTEGER,
        on_Ahash TEXT,
    );
    """
    con.execute(query);

    for dbname in db_paths:
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
            md5(CAST(allele.sequence AS TEXT)) AS on_Ahash
        FROM
            read_parquet('{allele_pq}') AS allele
        JOIN
            loci_lookup
        ON
            loci_lookup.dbname = '{dbname}'
            AND loci_lookup.update_LocusID = allele.LocusID
        """
        con.execute(query)

    # Set what the allele_number should be set to_ by looking at the destination database
    ba = base['allele']
    query = f"""
        UPDATE allele_lookup
        SET to_allele_number = dest.allele_number
        FROM read_parquet('{ba}') AS dest
        WHERE 
            allele_lookup.to_LocusID = dest.LocusID
            AND allele_lookup.on_Ahash = md5(CAST(dest.sequence AS TEXT))
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
    query = f"""
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
    query = f"""
    CREATE TABLE allele_pull AS
    SELECT DISTINCT ON (to_LocusID, to_allele_number_new)
        dbname,
        to_LocusID,
        to_allele_number_new,
        update_LocusID,
        update_allele_number,
    FROM allele_lookup
    WHERE to_allele_number IS NULL;

    SELECT DISTINCT ON (dbname) dbname FROM allele_pull
    """
    to_pull = con.execute(query).fetchall()

    # This might be helpful for MASSIVE merges. But probably not
    #con.execute("""
    #CREATE INDEX 
        #idx_allele_pull 
    #ON allele_pull(dbname, update_LocusID, update_allele_number);
    #
    #CREATE INDEX 
        #idx_allele_lookup 
    #ON allele_lookup(dbname, update_LocusID, update_allele_number);
    #""")

    # Make a temporary file holding the new allele entries
    query = """
    CREATE TABLE new_alleles (
        LocusID INTEGER,
        allele_number INTEGER,
        allele_length INTEGER,
        sequence BLOB,
    );
    """
    con.execute(query)

    logging.info("Updating allele")
    for dbname, in to_pull:
        logging.debug(dbname)
        names = tdb.get_tdb_filenames(dbname)
        m_allele = names['allele']
        query = f"""
            INSERT INTO new_alleles (LocusID, allele_number, allele_length, sequence)
            SELECT
                allele_pull.to_LocusID AS LocusID,
                allele_pull.to_allele_number_new AS allele_number,
                original.allele_length,
                original.sequence,
            FROM
                allele_pull
            JOIN
                read_parquet('{m_allele}') AS original 
                ON allele_pull.update_LocusID = original.LocusID
                AND allele_pull.update_allele_number = original.allele_number
            WHERE
                allele_pull.dbname = '{dbname}';
        """
        con.execute(query)

    # TODO: needs a do_order/compress
    alocus = os.path.join(output_dir, "allele.pq")
    query = f"""
    COPY (
        SELECT * 
        FROM read_parquet('{base["allele"]}')
        UNION ALL
        SELECT *
        FROM new_alleles
    ) TO '{alocus}' (FORMAT PARQUET)
    """
    con.execute(query)

def consolidate_samples(con, db_names):
    logging.info("Updating samples")
    # would need to copy from the base
    base = tdb.get_tdb_filenames(db_names[0])
    for sample_pq in base['sample'].values():
        out_name = os.path.join(output_dir, os.path.basename(sample_pq))
        shutil.copy(sample_pq, out_name)
        
    # And then update the rest
    for dbname in db_names[1:]:
        files = tdb.get_tdb_filenames(dbname)
        for sample, sample_pq in files['sample'].items():
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
                ) TO '{out_name}' (FORMAT PARQUET)
            """
            con.execute(query)

def merge_batch_main(args):
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

    os.mkdir(args.output)

    temp_db = truvari.make_temp_filename(suffix=".duckdb")

    con = duckdb.connect(temp_db)
    
    con.execute("SET default_null_order ='NULLS LAST';")
    con.execute("SET threads = 8;")
    con.execute("SET memory_limit = '40GB';")

    consolidate_locus(con, base, args.inputs, args.output)
    consolidate_alleles(con, base, args.inputs, args.output)
    consolidate_samples(con, base, args.inputs, args.output)
    con.close()
