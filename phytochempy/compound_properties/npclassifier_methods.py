import json
import os
import time
import urllib
import uuid

import numpy as np
import pandas as pd
from typing import List

import requests
from tqdm import tqdm

_NP_CLASSIFIER_COLUMNS = [
    'NPclassif_class_results', 'NPclassif_superclass_results',
    'NPclassif_pathway_results', 'NPclassif_isglycoside']


def get_npclassifier_result_columns_in_df(df: pd.DataFrame) -> List[str]:
    out = []
    for c in df.columns.tolist():
        if any(x in c for x in _NP_CLASSIFIER_COLUMNS):
            out.append(c)
    return out



def npclassify_smiles(smiles: str) -> dict:
    # From https://ccms-ucsd.github.io/GNPSDocumentation/api/
    # Function to classify a single SMILES string
    if smiles == smiles and smiles is not None:
        time.sleep(0.1)  # check rate limiting
        safe_string = urllib.parse.quote(smiles)
        url = f"https://npclassifier.gnps2.org/classify?smiles={safe_string}"
        response = requests.get(url)
        if response.status_code == 200:
            return_dict = json.loads(response.text)

            for k in ['class_results', 'superclass_results', 'pathway_results']:
                if len(return_dict[k]) == 0:
                    return_dict[k] = np.nan
                else:
                    # Some compounds are given multiple values for classes etc. Output these in separate columns
                    for i, val in enumerate(return_dict[k]):
                        return_dict[k + '_' + str(i)] = val
                    return_dict[k] = ':'.join(return_dict[k])

            for k in list(return_dict.keys())[:]:
                return_dict['NPclassif_' + k] = return_dict[k]
                del return_dict[k]

            return_dict['SMILES'] = smiles
            return return_dict
        else:

            return None
    else:
        return None


def get_npclassif_classes_from_smiles(smiles: List[str], npclassifier_cache_dir: str = None):
    unique_smiles = list(set(smiles))

    # First save time by reading previous temp outputs
    temp_file_tag = 'npclassifierinfo_cache_salt_'
    manual_file_tag = 'npclassifierinfo_manual_'
    existing_df = None
    if npclassifier_cache_dir is not None:
        existing_info = []
        for temp_cl_file in os.listdir(npclassifier_cache_dir):
            if temp_cl_file.startswith(temp_file_tag):
                existing_info.append(pd.read_csv(os.path.join(npclassifier_cache_dir, temp_cl_file), index_col=0))
            if temp_cl_file.startswith(manual_file_tag):
                existing_info.append(read_manual_npclassifier_input(os.path.join(npclassifier_cache_dir, temp_cl_file)))

        if len(existing_info) > 0:
            existing_df = pd.concat(existing_info)

            already_known_smiles = existing_df['SMILES'].tolist()
            # remove the item for all its occurrences
            for alread_known in already_known_smiles:
                c = unique_smiles.count(alread_known)
                for i in range(c):
                    unique_smiles.remove(alread_known)

    out_df = pd.DataFrame()
    failed_smiles = []
    for sm in tqdm(unique_smiles):
        if sm is not None:
            result = npclassify_smiles(sm)

            if result is not None:
                ent_df = pd.DataFrame(result, index=[0])
                out_df = pd.concat([out_df, ent_df])
            else:
                failed_smiles.append(sm)
    if npclassifier_cache_dir is not None:
        if len(out_df.index) > 0:
            out_df.to_csv(os.path.join(npclassifier_cache_dir, temp_file_tag + str(uuid.uuid4()) + '.csv'))
        if existing_df is not None:
            out_df = pd.concat([out_df, existing_df])
    out_df = out_df.sort_values(by='SMILES')
    out_df = out_df.reset_index(drop=True)

    if len(failed_smiles) > 0:
        failed_file = os.path.join(npclassifier_cache_dir, 'failed_smiles_for_manual_upload.csv')
        print(
            f'WARNING: some SMILES were not resolved to NPClassifier through the API. These will be saved to {failed_file}. '
            f'You can attempt to resolve these through the gnps portal.')
        print(
            f'To avoid repeating searches for these SMILES, you can save the manual output in: {npclassifier_cache_dir} with a file name starting with: {manual_file_tag}. Then rerun this function.')
        ## Some smiles aren't resolved through the API that seem to be resolved ok through the portal:
        # https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22NPCLASSIFIER%22%7D
        # See: https://ccms-ucsd.github.io/GNPSDocumentation/api/#structure-natural-product-classifier-np-classifier
        df = pd.DataFrame(failed_smiles, columns=['SMILES'])
        df = df[['SMILES']].dropna().drop_duplicates(
            keep='first')
        df.to_csv(failed_file)
    return out_df


def get_npclassifier_classes_from_df(df: pd.DataFrame, smiles_col: str, tempout_dir: str = None) -> pd.DataFrame:
    npclassifier_info = get_npclassif_classes_from_smiles(df[smiles_col].dropna(), tempout_dir)

    npclassifier_info[
        ['SMILES'] + _NP_CLASSIFIER_COLUMNS].drop_duplicates(keep='first').dropna(subset='SMILES')

    all_metabolites_with_class_info = pd.merge(df, npclassifier_info, how='left', on='SMILES')

    return all_metabolites_with_class_info


def read_manual_npclassifier_input(npclassifier_output_file: str):
    ### Get NPclassifier info
    np_classif_results = pd.read_csv(npclassifier_output_file,
                                     sep='\t').drop_duplicates(keep='first')
    rename_dict = {'smiles': 'SMILES'}
    for c in np_classif_results.columns:
        if c != 'smiles':
            rename_dict[c] = 'NPclassif_' + c
    np_classif_results = np_classif_results.rename(columns=rename_dict).dropna(subset='SMILES')

    for col in ['NPclassif_class_results', 'NPclassif_superclass_results',
                'NPclassif_pathway_results']:
        np_classif_results['col_split'] = np_classif_results[col].str.split(':')
        max_splits = np_classif_results['col_split'].str.len().max()
        if max_splits > 1:
            # Create new columns for each substring
            for i in range(max_splits):
                np_classif_results[f'{col}_{i}'] = np_classif_results['col_split'].str.get(i)
        else:
            np_classif_results[f'{col}_{0}'] = np_classif_results['col_split'].str.get(0)
    np_classif_results = np_classif_results.drop(columns=['col_split'])
    return np_classif_results
