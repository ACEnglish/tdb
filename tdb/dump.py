"""
Join all tables in a tdb and dump to a flat tsv file
"""
import logging
import truvari
import argparse
import pandas as pd

import tdb

def dump_main(args):
    """
    Main entrypoint
    """
    parser = argparse.ArgumentParser(prog="tdb dump", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", default="/dev/stdout",
                        help="Output file (stdout)")
    parser.add_argument("input", metavar="IN",
                        help="Input tdb")
    args = parser.parse_args(args)

    truvari.setup_logging()

    data = tdb.load_tdb(args.input)
    view = pd.merge(data['locus'], data['allele'], how='right', on="LocusID")
    use_header = True
    for sample, table in data['sample'].items():
        view2 = pd.merge(view, table, on=['LocusID', 'allele_number'])
        view2['sample'] = sample
        view2.to_csv(args.output, mode='w' if use_header else 'a', sep='\t', index=False, header=use_header)
        use_header=False
    logging.info("Finished")
