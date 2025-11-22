import os
import pathlib
import re
import time
import uuid
from typing import List

import cirpy
import numpy as np
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from standardiser import standardise
from tqdm import tqdm

COMPOUND_NAME_COLUMN = 'example_compound_name'


def standardise_smiles_to_MAIP_smiles(smiles: str) -> str:
    '''
    Molecular standardisation following MAIP procedure: https://chembl.gitbook.io/malaria-project/molecule-standardisation
    using the https://github.com/flatkinson/standardiser package.

    Should only be used when only interested in parent (bioactive) components.

    Process:
    - Break bonds to Group I or II metals
    - Neutralize charges by adding/removing protons
    - Apply standardization rules
    - Neutralise any charges exposed by rules
    - Discard any salt/solvate compounds
    - Return standardised parent
    :param smiles:
    :return:
    '''

    try:

        parent = standardise.run(smiles)

    except (standardise.StandardiseException, TypeError) as e:

        # print(f'SMILES standardisation warning: {e.message}')
        return None

    return parent


def standardise_SMILES(smiles: str) -> str:
    '''
    Use rdkit sanitzation procedure to standardise molecule and resolve to parent fragment.
    :param smiles:
    :return:
    '''
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        Chem.SanitizeMol(mol)
        parent_clean_mol = rdMolStandardize.FragmentParent(mol)

        return Chem.MolToSmiles(parent_clean_mol, isomericSmiles=True)
    except:
        return None


def resolve_cas_to_smiles(cas_id: str):
    """
    Resolves a given CAS ID to its corresponding SMILES representation.

    :param cas_id: A string representing the CAS ID.
    :return: A string representing the SMILES representation of the given CAS ID.

    """

    out = cirpy.resolve(cas_id, 'smiles')
    if out is None:
        out = None
    return out


def resolve_cas_to_inchikey(cas_id: str):
    """
    Resolve a CAS ID to its corresponding InChIKey.

    :param cas_id: The CAS ID to resolve.
    :return: The resolved InChIKey or np.nan if not found.
    """

    inch = cirpy.resolve(cas_id, 'stdinchikey')
    if inch is not None:
        out = inch.replace('InChIKey=', '')
    else:
        out = None
    return out


def get_smiles_and_inchi_from_cas_ids(cas_ids: List[str], tempout_dir: str = None):
    """
    :param cas_ids: A list of CAS IDs for which SMILES and InChIKey are to be retrieved.
    :param tempout_dir: The directory where temporary output files are stored. Default is None.
    :return: A pandas DataFrame containing CAS ID, SMILES, and InChIKey for each CAS ID in the input list.

    """
    print(f'Warning this method predominantly relies on CAS ID translations based on the Knapsack database.'
          f'This is likely appropriate for data downloaded from Knapsack but may not be ideal for other sources. ')
    unique_cas_ids = list(set(cas_ids))

    # First save time by reading previous temp outputs
    temp_file_tag = 'cirpy_cas_cache_salt_'
    existing_df = None
    if tempout_dir is not None:
        existing_info = []
        for temp_cl_file in os.listdir(tempout_dir):
            if temp_cl_file.startswith(temp_file_tag):
                existing_info.append(pd.read_csv(os.path.join(tempout_dir, temp_cl_file), index_col=0))
        if len(existing_info) > 0:
            existing_df = pd.concat(existing_info)

            already_known_cas_ids = set(existing_df['CAS ID'].tolist())
            unique_cas_ids = [s for s in unique_cas_ids if s not in already_known_cas_ids]
    cache_csv = os.path.join(tempout_dir, temp_file_tag + str(uuid.uuid4()) + '.csv')
    out_df = pd.DataFrame()

    for c_id in tqdm(unique_cas_ids, desc='Resolving CAS IDs..'):
        inch_result, smiles_result = get_compound_ids_from_CAS_ID_from_knapsack(c_id)
        if smiles_result is None:
            smiles_result = resolve_cas_to_smiles(c_id)
        if inch_result is None:
            inch_result = resolve_cas_to_inchikey(c_id)
        result = {'CAS ID': c_id, 'SMILES': smiles_result, 'InChIKey': inch_result}
        ent_df = pd.DataFrame(result, index=[0])
        out_df = pd.concat([out_df, ent_df])

        if tempout_dir is not None:
            if len(out_df.index) > 0:
                out_df.to_csv(cache_csv)
    if existing_df is not None:
        out_df = pd.concat([out_df, existing_df])
    out_df = out_df.sort_values(by='CAS ID')
    out_df = out_df.reset_index(drop=True)
    return out_df


def add_CAS_ID_translations_to_df(df: pd.DataFrame, cas_id_col: str, tempout_dir: str = None) -> pd.DataFrame:
    """
    :param df: DataFrame containing the data
    :param cas_id_col: Name of the column in the DataFrame that contains the CAS IDs
    :param tempout_dir: Optional directory path to store temporary files (default: None)
    :return: DataFrame with CAS ID translations added

    This method takes a DataFrame and adds CAS ID translations to it. The translations are obtained by calling the 'get_smiles_and_inchi_from_cas_ids' method.

    """
    pathlib.Path(tempout_dir).mkdir(parents=True, exist_ok=True)
    _info = get_smiles_and_inchi_from_cas_ids(df[cas_id_col].dropna(), tempout_dir)

    _info[[cas_id_col, 'SMILES', 'InChIKey']].drop_duplicates(keep='first').dropna(subset=cas_id_col)

    all_metabolites_with_info = pd.merge(df, _info, how='left', on=cas_id_col)

    return all_metabolites_with_info


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


def sanitize_filename(pathway: str, replace_spaces: bool = True) -> str:
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


def fill_match_ids(df: pd.DataFrame, given_col: str) -> pd.DataFrame:
    """
    :param df: A pandas DataFrame containing chemical data.
    :param given_col: The name of the column in the DataFrame that needs to be filled with matching compound IDs.
    :return: A pandas DataFrame with the given_col filled with matching compound IDs.

    This method takes a DataFrame and a column name as parameters.
    It removes cases with no identifiers and then fills in the matching details for the given_col column from other columns
    * in the DataFrame.

    """
    chem_id_cols = ['SMILES', 'InChIKey', 'CAS ID']
    # Remove cases with no identifiers
    df = df.copy(deep=True).dropna(subset=chem_id_cols, how='all')
    ### Fill in matching details for compound ID columns from other dataframes
    cols_to_do = [cl for cl in chem_id_cols if cl != given_col]
    # Step 1: Create a mappings to match known pairs of IDs
    mapping1 = {}
    mapping2 = {}
    for index, row in df.iterrows():
        if row[given_col] == row[given_col]:  # If not nan
            if row[cols_to_do[0]] == row[cols_to_do[0]]:
                mapping1[row[cols_to_do[0]]] = row[given_col]
            if row[cols_to_do[1]] == row[cols_to_do[1]]:
                mapping2[row[cols_to_do[1]]] = row[given_col]

    if any(x in mapping1.keys() for x in ['', np.nan, 'nan']):
        raise ValueError
    if any(x in mapping2.keys() for x in ['', np.nan, 'nan']):
        raise ValueError

    # Step 2: Iterate through the DataFrame to fill in the empty given_col values
    for index, row in df.iterrows():
        if row[given_col] != row[given_col]:
            m1 = row[cols_to_do[0]]
            m2 = row[cols_to_do[1]]
            # prioritise inchikey and smiles
            if m1 in mapping1:
                v = mapping1[m1]
                df.at[index, given_col] = v
            elif m2 in mapping2:
                v = mapping2[m2]
                df.at[index, given_col] = v

    return df


def get_compound_ids_from_CAS_ID_from_knapsack(cas_id: str):
    """
    Fetches compound identifiers (InChIKey and SMILES) associated with a given CAS ID from the Knapsack database.

    This function constructs a URL to query the Knapsack database with a specified CAS ID. It retrieves
    data in the form of HTML tables, parses them, and extracts the relevant compound information such as
    InChIKey and SMILES. The function is robust to cases where the underlying web structure may have
    encoding issues and handles such scenarios appropriately.

    Raises an IndexError if the expected data structure is not present in the HTML content or the
    specific compound information cannot be found.

    Parameters:
    cas_id: str
        The CAS ID of the compound for which the InChIKey and SMILES need to be fetched.

    Returns:
    tuple[str | None, str | None]
        A tuple containing:
        - The InChIKey of the compound, or None if not available.
        - The SMILES of the compound, or None if not available.
    """
    if cas_id == '' or cas_id is None:
        return None, None
    url = f'http://www.knapsackfamily.com/knapsack_core/information.php?sname=CAS_ID&word={cas_id}'
    time.sleep(.01)
    try:

        tables = pd.read_html(url, flavor='html5lib')

    except UnicodeEncodeError:
        response = requests.get(url)
        decoded = response.content.decode()
        tables = pd.read_html(decoded, flavor='html5lib')

    except ValueError:
        return None, None
    meta_table = tables[0]

    inchikey_df = meta_table[meta_table['Metabolite Information'] == 'InChIKey']
    try:
        inchikey = inchikey_df['Structural formula'].values[0]
    except IndexError:
        inchikey = None
    smiles_df = meta_table[meta_table['Metabolite Information'] == 'SMILES']
    try:
        smiles = smiles_df['Structural formula'].values[0]
    except IndexError:
        smiles = None
    return inchikey, smiles
