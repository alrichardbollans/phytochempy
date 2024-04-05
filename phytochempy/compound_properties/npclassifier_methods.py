import json
import os
import time
import uuid

import numpy as np
import pandas as pd
from typing import List

import requests
from tqdm import tqdm

NP_CLASSIFIER_COLUMNS = [
    'NPclassif_class_results', 'NPclassif_superclass_results',
    'NPclassif_pathway_results', 'NPclassif_isglycoside']


def npclassify_smiles(smiles: str) -> dict:
    # Function to classify a single SMILES string
    if smiles == smiles and smiles is not None:
        time.sleep(0.1)  # check rate limiting
        url = f"https://npclassifier.gnps2.org/classify?smiles={smiles}"
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
            print(f'Warning: No npClassifier response for {smiles}')
            return None
    else:
        return None


def get_npclassif_classes_from_smiles(smiles: List[str], tempout_dir: str = None):
    unique_smiles = list(set(smiles))

    # First save time by reading previous temp outputs
    existing_df = None
    if tempout_dir is not None:
        existing_info = []
        for temp_cl_file in os.listdir(tempout_dir):
            if temp_cl_file.startswith('npclassifierinfo_'):
                existing_info.append(pd.read_csv(os.path.join(tempout_dir, temp_cl_file), index_col=0))
        if len(existing_info) > 0:
            existing_df = pd.concat(existing_info)

            already_known_smiles = existing_df['SMILES'].tolist()
            # remove the item for all its occurrences
            for alread_known in already_known_smiles:
                c = unique_smiles.count(alread_known)
                for i in range(c):
                    unique_smiles.remove(alread_known)

    out_df = pd.DataFrame()

    for sm in tqdm(unique_smiles):
        result = npclassify_smiles(sm)
        if result is not None:
            ent_df = pd.DataFrame(result, index=[0])
            out_df = pd.concat([out_df, ent_df])

    if tempout_dir is not None:
        if len(out_df.index) > 0:
            out_df.to_csv(os.path.join(tempout_dir, 'npclassifierinfo_' + str(uuid.uuid4()) + '.csv'))
        if existing_df is not None:
            out_df = pd.concat([out_df, existing_df])
    out_df = out_df.sort_values(by='SMILES')
    out_df = out_df.reset_index(drop=True)
    return out_df


def get_npclassifier_classes_from_df(df: pd.DataFrame, smiles_col: str, tempout_dir: str = None) -> pd.DataFrame:
    npclassifier_info = get_npclassif_classes_from_smiles(df[smiles_col].dropna(), tempout_dir)

    npclassifier_info[
        ['SMILES'] + NP_CLASSIFIER_COLUMNS].drop_duplicates(keep='first').dropna(subset='SMILES')

    all_metabolites_with_class_info = pd.merge(df, npclassifier_info, how='left', on='SMILES')

    return all_metabolites_with_class_info
