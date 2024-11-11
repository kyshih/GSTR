# shrinkage_analysis/data_utils.py

import pandas as pd

def get_excluded_samples_from_file(file_path: str):
    """Reads a text file and returns a list of excluded sample IDs."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def df_to_dict(df, key_col, value_col):
    """Converts two columns of a DataFrame to a dictionary."""
    return df.set_index(key_col)[value_col].to_dict()

def filter_df(df: pd.DataFrame, column: str, values: list, exclude: bool = True):
    """
    Filters a DataFrame based on a list of values in a specific column.
    If `exclude` is True, it removes rows where `column` is in `values`.
    Default returns df after excluding unwanted rows
    """
    if exclude:
        return df[~df[column].isin(values)]
    else:
        return df[df[column].isin(values)]

def find_unique_values(df: pd.DataFrame, column: str):
    """Returns a set of unique values from a specified column in the DataFrame."""
    return set(df[column].unique())
