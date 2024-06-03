"""
Join all tables in a tdb and dump to a flat tsv file
"""
import logging
import argparse
import pandas as pd

import tdb

def dump_tdb(data, file=None):
    """
    Full join of samples. If a file is provided, it
    is iteratively written to. If a file isn't provided,
    the DataFrame is returned, otherwise an empty dataframe
    """
    data['allele']['sequence'] = data['allele']['sequence'].str.decode('utf8')
    view = pd.merge(data['locus'], data['allele'], how='right', on="LocusID")
    if len(data['sample']) == 0:
        if file:
            view.to_csv(file, mode='w', sep='\t', index=False)
        return view

    use_header = file is not None
    parts = []
    for sample, table in data['sample'].items():
        view2 = pd.merge(view, table, on=['LocusID', 'allele_number'])
        view2['sample'] = sample
        if file:
            view2.to_csv(file, mode='w' if use_header else 'a', sep='\t', index=False, header=use_header)
        parts.append(view2)
        use_header=False

    if parts:
        return pd.concat(parts)
    return view

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

    tdb.setup_logging()

    data = tdb.load_tdb(args.input)
    dump_tdb(data, args.output)
    logging.info("Finished")
