"""
De-identify a tdb
"""
import os
import sys
import shutil
import logging
import argparse
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import tdb
from tdb.create import S_SCHEMA


def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False
    if not os.path.exists(args.input):
        logging.error(f"input {args.input} does not exists")
        check_fail = True
    if not args.input.rstrip('/').endswith(".tdb"):
        logging.error(f"input {args.input} must end with `.tdb`")
        check_fail = True
    if os.path.exists(args.output):
        logging.error(f"output {args.output} exists")
        check_fail = True
    if not args.output.rstrip('/').endswith(".tdb"):
        logging.error(f"output {args.output} must end with `.tdb`")
        check_fail = True
    return check_fail


def deid_main(args):
    """
    Remove genotypes from a tdb
    """
    parser = argparse.ArgumentParser(prog="tdb append", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # these could be positional arguments, but a little bit of user friction
    # will help prevent unintentional overwriting
    parser.add_argument("-i", "--input", metavar="IN", type=str,
                        help="Input DB")
    parser.add_argument("-o", "--output", metavar="OUT", type=str,
                        help="Output DB")
    parser.add_argument("-s", "--remove-sequences", action="store_true",
                        help="Remove sequences from allele table")
    parser.add_argument("-S", "--shuffle-samples", action="store_true",
                        help="Shuffle sample tables together and split (experimental)")
    args = parser.parse_args(args)
    tdb.setup_logging()

    if check_args(args):
        logging.error("argument error. exiting")
        sys.exit(1)

    in_file_names = tdb.get_tdb_filenames(args.input)
    out_file_names = tdb.get_tdb_filenames(args.output)
    os.mkdir(args.output)
    logging.info("Copying locus")
    shutil.copy(in_file_names['locus'], out_file_names['locus'])

    alleles = pd.read_parquet(in_file_names['allele'])
    if args.remove_sequences:
        alleles["sequence"] = ""

    logging.info("Changing allele_numbers")

    ref = alleles[alleles["allele_number"] == 0].copy()
    alt = alleles[alleles["allele_number"] != 0].copy().sort_values(["LocusID", "allele_length"])

    alt["allele_number"] = alt.groupby(["LocusID"]).cumcount() + 1
    out = pd.concat([ref, alt]).sort_values(["LocusID", "allele_number"])
    out.to_parquet(out_file_names['allele'], index=False, compression='gzip')

    if args.shuffle_samples:
        # For testing, we want deterministic shuffling
        if "TDB_SEED" in os.environ and os.environ["TDB_SEED"] == "123":
            logging.critical('seed')
            seed = 123
        else:
            seed = None
        logging.info("Shuffling samples (experimental)")
        parts = []
        n_samps = 0
        for i in in_file_names['sample'].values():
            parts.append(pd.read_parquet(i))
            n_samps += 1

        samps = pd.concat(parts).sample(
            frac=1, random_state=seed).sort_values(["LocusID"])

        for idx in range(n_samps):
            o_fn = os.path.join(args.output, f"sample.{idx}.pq")
            value = samps.iloc[idx:len(samps):n_samps]
            m_table = pa.Table.from_pandas(value, schema=S_SCHEMA, preserve_index=False)
            writer = pq.ParquetWriter(o_fn, S_SCHEMA, compression='gzip')
            writer.write(m_table)
            writer.close()

    logging.info("Finished")
