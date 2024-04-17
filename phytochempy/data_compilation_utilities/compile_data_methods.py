import os.path
from typing import List

import pandas as pd
from wcvpy.wcvp_download import wcvp_accepted_columns
from wcvpy.wcvp_name_matching import get_genus_from_full_name, output_record_col_names

from phytochempy.compound_properties import simplify_inchi_key, COMPOUND_NAME_COLUMN, fill_match_ids


def merge_and_tidy_compound_datasets(datasets: List[pd.DataFrame], output_csv: str):
    """
    Given datasets from WikiData and Knapsack, merge them and tidy them.

    :param datasets: A list of datasets to merge.
    :param output_csv: The path to the output CSV file.
    :return: The merged and tidied dataset.
    """
    all_metabolites_in_taxa = pd.concat(datasets)
    all_metabolites_in_taxa = all_metabolites_in_taxa.dropna(subset=wcvp_accepted_columns['name_w_author'])
    ### Format
    all_metabolites_in_taxa = all_metabolites_in_taxa.sort_values(
        by=wcvp_accepted_columns['name']).reset_index(drop=True)

    output_columns = output_record_col_names + [COMPOUND_NAME_COLUMN, 'SMILES', 'InChIKey', 'CAS ID', 'Source']

    all_metabolites_in_taxa = all_metabolites_in_taxa[output_columns]

    start_cols = ['accepted_name_w_author', COMPOUND_NAME_COLUMN, 'InChIKey', 'SMILES', 'CAS ID', 'Source']

    all_metabolites_in_taxa = all_metabolites_in_taxa[
        start_cols + [col for col in all_metabolites_in_taxa.columns if
                      col not in start_cols]]

    # Tidy final list

    for c_id in ['SMILES', 'InChIKey', 'CAS ID']:
        all_metabolites_in_taxa = fill_match_ids(all_metabolites_in_taxa, c_id)

    all_metabolites_in_taxa['InChIKey_simp'] = all_metabolites_in_taxa['InChIKey'].apply(simplify_inchi_key)
    all_metabolites_in_taxa.to_csv(output_csv)
    return all_metabolites_in_taxa


def get_manual_MAIP_to_upload(df: pd.DataFrame, _temp_output_path):
    ''' This generates files to upload to obtain information, which can be merged in later steps.'''

    ### Get MAIP info
    # https://www.ebi.ac.uk/chembl/maip/
    all_metabolites_to_send_to_maip = df[['SMILES']].dropna().drop_duplicates(
        keep='first')
    all_metabolites_to_send_to_maip['id'] = all_metabolites_to_send_to_maip['SMILES']
    all_metabolites_to_send_to_maip.to_csv(os.path.join(_temp_output_path, 'smiles_for_MAIP.csv'))


def add_manual_info_files(df: pd.DataFrame, maip_output_file: str = None):
    all_metabolites_with_info = df.copy(deep=True)

    if maip_output_file is not None:
        ### Get MAIP values
        maip_results = pd.read_csv(maip_output_file)
        maip_results = maip_results.rename(columns={'smiles': 'SMILES', 'model_score': 'MAIP_model_score'})
        maip_results = maip_results[['SMILES', 'MAIP_model_score']].dropna(subset='SMILES')
        all_metabolites_with_info = pd.merge(all_metabolites_with_info, maip_results[~maip_results['SMILES'].isna()], how='left',
                                             on='SMILES')

    return all_metabolites_with_info


def _checks(df):
    ### Checks
    def has_different_maipscore(group):
        return group['MAIP_model_score'].nunique() > 1

    # Check cases with the same simplified inchikey_simp value have same SMILEs/MAIP score. This is overwhelmingly the case. Similarly for compound classes.
    inchi_problems1 = pre_final_df.groupby(compound_id_col).filter(has_different_maipscore).drop_duplicates(
        subset=[compound_id_col, 'SMILES'],
        keep='first').sort_values(
        by=compound_id_col)

    if len(inchi_problems1.index) > 0:
        inchi_problems1.to_csv(os.path.join(_temp_output_path, 'same_inchisimp_diff_smiles.csv'))
        print(
            f'WARNING: Same some compounds with same InChIsimp and different smiles. See {os.path.join(_temp_output_path, "same_inchisimp_diff_smiles.csv")}')


def tidy_final_dataset(pre_final_df: pd.DataFrame, _temp_output_path: str, final_taxa_compound_csv: str, compound_id_col: str) -> None:
    """
    :param pre_final_df: A pandas DataFrame containing the pre-final dataset.
    :param _temp_output_path: The path to temporarily store the output.
    :param final_taxa_compound_csv: The path to save the final taxa compound CSV file.
    :param compound_id_col: The column name of the compound ID in the DataFrame.
    :return: None

    Tidies the pre-final dataset by dropping rows with missing values in the compound ID column or plant accepted name column.
    Adds a 'Genus' column by applying the 'get_genus_from_full_name' function to the accepted name column.
    Drops duplicate rows based on the compound ID and compound name with author columns.
    Sorts the DataFrame by plant name with author and compound name column. Saves the sorted DataFrame to a CSV file.
    """
    ## Tidy data a bit
    pre_final_df = pre_final_df.dropna(subset=[compound_id_col, wcvp_accepted_columns['name_w_author']], how='any')
    # Add genus column
    pre_final_df['Genus'] = pre_final_df[wcvp_accepted_columns['name']].apply(
        get_genus_from_full_name)

    def drop_repeat_id(df, id_col):
        # Remove known duplicates according to ID
        df = df[df[id_col].isna() | (~df.duplicated(subset=[id_col, wcvp_accepted_columns['name_w_author']], keep='first'))]
        return df

    pre_final_df = drop_repeat_id(pre_final_df, compound_id_col)
    taxa_with_all_info = pre_final_df.sort_values(by=[wcvp_accepted_columns['name_w_author'], COMPOUND_NAME_COLUMN]).reset_index(drop=True)

    taxa_with_all_info.to_csv(final_taxa_compound_csv)
