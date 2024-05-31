"""
Create a tdb
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

    if os.path.exists(args.output):
        logging.error(f"Output  {args.output} already exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"Output {args.output} must end with `.tdb`")
        check_fail = True
    if not os.path.exists(args.input):
        logging.error(f"Input {args.input} does not exist")
        check_fail = True
    if not args.input.rstrip('/').endswith((".vcf", ".vcf.gz")):
        logging.error(f"Unrecognized file extension on {args.input}. Expected .vcf .vcf.gz")
        check_fail = True
    return check_fail

def create_main(args):
    """
    Create a new tdb from multiple input calls
    """
    parser = argparse.ArgumentParser(prog="tdb create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output tdb directory")
    parser.add_argument("input", metavar="IN",
                        help="VCF file")
    args = parser.parse_args(args)

    truvari.setup_logging()
    if check_args(args):
        logging.error("Cannot create database. Exiting")
        sys.exit(1)

    logging.info("Loading %s", args.input)
    m_data = tdb.vcf_to_tdb(args.input)
    logging.info("Writing parquet files")
    tdb.dump_tdb(m_data, args.output)
    logging.info("Finished")
