"""
Development debugging utility to check if two tdbs are identical
"""
import argparse
import sys
import pandas as pd
import tdb


def dataframes_equal(df1, df2, path):
    """Check if two DataFrames are equal."""
    df1 = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
    df2 = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)
    if not df1.equals(df2):
        print(f"DataFrames differ at {path}", file=sys.stderr)
        return False
    return True


def dictionaries_equal(dict1, dict2, path):
    """Recursively check if two dictionaries are equal."""
    if dict1.keys() != dict2.keys():
        print(
            f"Dictionary keys differ at {path}: {dict1.keys()} != {dict2.keys()}", file=sys.stderr)
        return False

    for key in dict1:
        new_path = f"{path}.{key}"
        if isinstance(dict1[key], pd.DataFrame) and isinstance(dict2[key], pd.DataFrame):
            if not dataframes_equal(dict1[key], dict2[key], new_path):
                return False
        elif isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
            if not dictionaries_equal(dict1[key], dict2[key], new_path):
                return False
        else:
            if dict1[key] != dict2[key]:
                print(
                    f"Values differ at {new_path}: {dict1[key]} != {dict2[key]}", file=sys.stderr)
                return False

    return True


def check_dicts_equal(dict1, dict2):
    """
    Check if two dictionaries with keys 'locus', 'allele', and 'sample' are the same.
    'locus' and 'allele' values are DataFrames.
    'sample' value is a dictionary of sample_names: DataFrames.
    """
    if not isinstance(dict1, dict) or not isinstance(dict2, dict):
        print("Either dict1 or dict2 is not a dictionary", file=sys.stderr)
        return False

    required_keys = {'locus', 'allele', 'sample'}

    if set(dict1.keys()) != required_keys or set(dict2.keys()) != required_keys:
        print(
            f"Dictionary keys do not match required keys at root: {dict1.keys()} != {dict2.keys()}", file=sys.stderr)
        return False

    # Check locus and allele
    if not dataframes_equal(dict1['locus'], dict2['locus'], 'locus'):
        return False
    if not dataframes_equal(dict1['allele'], dict2['allele'], 'allele'):
        return False

    # Check sample
    if not dictionaries_equal(dict1['sample'], dict2['sample'], 'sample'):
        return False

    return True


def join_cmp(db1, db2, strip):
    """
    Dump the databases and compare
    """
    dump1 = tdb.dump_tdb(db1)
    dump2 = tdb.dump_tdb(db2)
    if strip:
        dump1.drop(columns=["LocusID", "allele_number"], inplace=True)
        dump2.drop(columns=["LocusID", "allele_number"], inplace=True)

    dump1.sort_values(by=list(dump1.columns), inplace=True)
    dump1.reset_index(drop=True, inplace=True)
    dump2.sort_values(by=list(dump2.columns), inplace=True)
    dump2.reset_index(drop=True, inplace=True)

    if not dump1.equals(dump2):
        print("Dumps differ", file=sys.stderr)
        return False
    return True


def debug_check_main(args):
    """
    Tool for checking if two tdbs are equal
    """
    parser = argparse.ArgumentParser(prog="tdb equal", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dbA", metavar="A",
                        help="First tdb")
    parser.add_argument("dbB", metavar="B",
                        help="Second tdb")
    parser.add_argument("--strip", action="store_true",
                        help="Don't compare LocusID, allele_number keys")
    parser.add_argument("--join", action="store_true",
                        help="Do a full join comparison")
    args = parser.parse_args(args)
    db1 = tdb.load_tdb(args.dbA)
    db2 = tdb.load_tdb(args.dbB)

    if args.join:
        join_cmp(db1, db2, args.strip)

    if args.strip:
        db1['locus'].drop(columns=['LocusID'], inplace=True)
        db2['locus'].drop(columns=['LocusID'], inplace=True)
        db1['allele'].drop(columns=['LocusID', "allele_number"], inplace=True)
        db2['allele'].drop(columns=['LocusID', "allele_number"], inplace=True)
        for sample in db1['sample']:
            db1['sample'][sample].drop(
                columns=['LocusID', "allele_number"], inplace=True)
        for sample in db2['sample']:
            db2['sample'][sample].drop(
                columns=['LocusID', "allele_number"], inplace=True)

    if check_dicts_equal(db1, db2):
        print("Databases equal", file=sys.stderr)
        sys.exit(0)
    else:
        sys.exit(1)
