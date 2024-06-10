"""
Checks if two parquet files are equal
"""
import sys
import pandas as pd

def compare_parquet_files(file1, file2):
    df1 = pd.read_parquet(file1)
    df2 = pd.read_parquet(file2)

    if df1.shape != df2.shape:
        return False

    if set(df1.columns) != set(df2.columns):
        return False

    return df1.equals(df2)

# Example usage
file1 = sys.argv[1]
file2 = sys.argv[2]

if compare_parquet_files(file1, file2):
    print("The Parquet files are equal.")
    sys.exit(0)
else:
    print("The Parquet files are not equal.")
    sys.exit(1)

