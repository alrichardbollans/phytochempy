import re
import urllib
from typing import List
from urllib.error import HTTPError

import cirpy
import numpy as np
import pandas as pd

COMPOUND_NAME_COLUMN = 'example_compound_name'

def resolve_cas_to_smiles(cas_id: str):
    """
    Resolves a given CAS ID to its corresponding SMILES representation.

    :param cas_id: A string representing the CAS ID.
    :return: A string representing the SMILES representation of the given CAS ID.

    """
    try:
        print(cas_id)
        out = cirpy.resolve(cas_id, 'smiles')
    except (urllib.error.HTTPError, urllib.error.URLError):
        out = np.nan
        if cas_id == cas_id and cas_id != '':
            print(f'WARNING: cas id not resolved: {cas_id}')
    return out


def resolve_cas_to_inchikey(cas_id: str):
    """
    Resolve a CAS ID to its corresponding InChIKey.

    :param cas_id: The CAS ID to resolve.
    :return: The resolved InChIKey or np.nan if not found.
    """
    try:
        print(cas_id)
        inch = cirpy.resolve(cas_id, 'stdinchikey')
        if inch is not None:
            out = inch.replace('InChIKey=', '')
        else:
            out = np.nan
    except (urllib.error.HTTPError, urllib.error.URLError):
        out = np.nan
        if cas_id == cas_id and cas_id != '':
            print(f'WARNING: cas id not resolved: {cas_id}')
    return out


def simplify_inchi_key(inch: str):
    """
    :param inch: The original InChIKey string.
    :return: The simplified InChIKey string, using only the first 14 characters of the original InChIKey.
    """
    # Using the connectivity layer of the InChIKey, i.e. the first 14 characters, to simplify.
    # As in e.g. https://www.sciencedirect.com/science/article/abs/pii/S2352007822002372 https://pubs.acs.org/doi/abs/10.1007/s13361-016-1589-4
    if inch == inch and inch is not None:
        return inch[:14]
    else:
        return inch

def sanitize_filename(pathway: str, replace_spaces: bool = True)->str:
    """
    Sanitize given pathway/compound names as they are currently not very nice for handling filenames and/or importing exporting in R etc..

    Args:
        pathway (str): The filename to sanitize.
        replace_spaces (bool): Whether to replace spaces with underscores.

    Returns:
        str: The sanitized filename.
    """
    if pathway == pathway:
        # Remove leading and trailing whitespace
        pathway = pathway.strip()

        # Replace spaces with underscores (optional)
        if replace_spaces:
            pathway = pathway.replace(" ", "_")

        # Remove illegal characters
        sanitized_filename = re.sub(r'[\\/*?:"<>|]', "", pathway)

        # Limit filename length (optional, for example, to 255 characters)
        sanitized_filename = sanitized_filename[:255]

        return sanitized_filename
    else:
        return pathway


def filter_rows_containing_compound_keyword(df: pd.DataFrame, cols: List[str], keyword: str):
    """
    :param df: A pandas DataFrame containing the data to filter.
    :param cols: A list of column names to search for the keyword.
    :param keyword: The keyword to search for in the specified columns.
    :return: A tuple of three pandas DataFrames:
            - filtered_df: A DataFrame containing the rows where the keyword was found in at least one of the specified columns.
            - negative_filtered_df: A DataFrame containing the rows where none of the specified columns matched the keyword.
            - nan_df: A DataFrame containing the rows where the values in all specified columns are empty strings or NaN.

    """
    # Initialize an empty DataFrame to store the filtered rows
    filtered_df = pd.DataFrame()
    negative_filtered_df = pd.DataFrame()
    nan_df = pd.DataFrame()
    # Iterate through the rows of the original DataFrame
    for index, row in df.iterrows():
        # Iterate through the specified columns
        if any(keyword in str(row[col]).lower() for col in cols):
            # If the keyword is found, add the entire row to the filtered DataFrame
            filtered_df = pd.concat([filtered_df, df.iloc[[index]]])

        elif all((str(row[col]) == '' or str(row[col]) == 'nan' or row[col] != row[col]) for col in cols):
            nan_df = pd.concat([nan_df, df.iloc[[index]]])
        else:
            negative_filtered_df = pd.concat([negative_filtered_df, df.iloc[[index]]])
    # Reset the index of the filtered DataFrame
    filtered_df.reset_index(drop=True, inplace=True)
    negative_filtered_df.reset_index(drop=True, inplace=True)
    nan_df.reset_index(drop=True, inplace=True)

    assert len(df.index) == (len(filtered_df.index) + len(negative_filtered_df.index)) + len(nan_df.index)

    return filtered_df, negative_filtered_df, nan_df
