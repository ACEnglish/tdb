"""
Development debugging utility to check if two tdbs are identical
"""
import sys
import pandas as pd
import tdb

def dataframes_equal(df1, df2, path):
    """Check if two DataFrames are equal."""
    keys = ["LocusID"]
    if path != 'locus':
        keys.append("allele_number")
    df1 = df1.sort_values(by=keys).set_index(keys)
    df2 = df2.sort_values(by=keys).set_index(keys)
    if not df1.equals(df2):
        print(f"DataFrames differ at {path}", file=sys.stderr)
        return False
    return True

def dictionaries_equal(dict1, dict2, path):
    """Recursively check if two dictionaries are equal."""
    if dict1.keys() != dict2.keys():
        print(f"Dictionary keys differ at {path}: {dict1.keys()} != {dict2.keys()}", file=sys.stderr)
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
                print(f"Values differ at {new_path}: {dict1[key]} != {dict2[key]}", file=sys.stderr)
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
        print(f"Dictionary keys do not match required keys at root: {dict1.keys()} != {dict2.keys()}", file=sys.stderr)
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

def debug_check_main(args):
    """
    Tool for checking if tdbs are equal
    """
    db1 = tdb.load_tdb(args[0])
    db2 = tdb.load_tdb(args[1])
    if check_dicts_equal(db1, db2):
        sys.exit(0)
    else:
        sys.exit(1)
