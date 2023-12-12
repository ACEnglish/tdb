"""
Add new VCF to existing database
"""
import os
import sys
import logging
import argparse
import truvari

import tdb

def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False
    if not os.path.exists(args.to):
        logging.error(f"output {args.to} does not exists")
        check_fail = True
    if not args.to.endswith(".tdb"):
        logging.error(f"output {args.to} must end with `.tdb`")
        check_fail = True
    if not os.path.exists(args.fr):
        logging.error(f"input {args.fr} does not exist")
        check_fail = True
    if not args.fr.endswith((".vcf", ".vcf.gz", ".tdb")):
        logging.error("unrecognized file extension on {i}")
        logging.error("expected .vcf .vcf.gz or .tdb")
        check_fail = True
    return check_fail

def append_main(args):
    """
    Add new vcf or tdb to an existing one
    """
    parser = argparse.ArgumentParser(prog="tdb append", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    # these could be positional arguments, but a little bit of user friction
    # will help prevent unintentional overwriting
    parser.add_argument("--to", metavar="DB", type=str,
                        help="Existing DB to pull into")
    parser.add_argument("--fr", metavar="FR", type=str,
                        help="Input vcf or db to pull from")
    args = parser.parse_args(args)
    truvari.setup_logging()

    if check_args(args):
        logging.error("cannot append to database. exiting")
        sys.exit(1)

    exist_db = tdb.load_tdb(args.to)
    new_db = tdb.load_tdb(args.fr) if args.fr.endswith(".tdb") else tdb.vcf_to_tdb(args.fr)
    result = tdb.tdb_consolidate(exist_db, new_db)
    tdb.dump_tdb(result, args.to)
    logging.info("Finished")
