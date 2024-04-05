import os.path

import numpy as np
import pandas as pd
from wcvp_download import wcvp_accepted_columns
from wcvp_name_matching import get_genus_from_full_name, output_record_col_names

from phytochempy.compound_properties import get_classyfire_classes_from_df, get_compound_info_from_chembl_apm_assays, \
    simplify_inchi_key, add_chembl_apm_data_to_compound_df, add_bioavailability_rules_to_df, COMPOUND_NAME_COLUMN, get_npclassifier_classes_from_df
from phytochempy.knapsack_searches import get_knapsack_compounds_in_family, tidy_knapsack_results
from phytochempy.wikidata_searches import generate_wikidata_search_query, submit_query, tidy_wikidata_output


def get_wikidata(wiki_data_id: str, temp_output_csv: str, tidied_output_csv: str, limit: int = 100000):
    # Example usage
    my_query = generate_wikidata_search_query(wiki_data_id, limit)
    submit_query(my_query, temp_output_csv, limit)
    tidy_wikidata_output(temp_output_csv, tidied_output_csv)


def get_knapsack_data(families_of_interest: list, temp_output_path: str, tidied_output_csv: str, add_smiles_and_inchi: bool = True):
    def _temp_out_for_fam(faml: str) -> str:
        return os.path.join(temp_output_path, faml + '_kn_search.csv')

    def _temp_out_for_fam_Acc(faml: str) -> str:
        return os.path.join(temp_output_path, faml + '_kn_search_accepted_info.csv')

    for fam in families_of_interest:
        get_knapsack_compounds_in_family(fam, _temp_out_for_fam(fam))
        tidy_knapsack_results(_temp_out_for_fam(fam), _temp_out_for_fam_Acc(fam), fam, cirpy_cache_dir=temp_output_path,
                              add_smiles_and_inchi=add_smiles_and_inchi)

    all_kn_dfs = pd.DataFrame()

    for fam in families_of_interest:
        new_df = pd.read_csv(_temp_out_for_fam_Acc(fam), index_col=0)
        all_kn_dfs = pd.concat([all_kn_dfs, new_df])

    all_kn_dfs.to_csv(tidied_output_csv)


def merge_and_tidy_compound_datasets(datasets: list, output_csv: str):
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
    chem_id_cols = ['SMILES', 'InChIKey', 'CAS ID']
    # Remove cases with no identifiers
    all_metabolites_in_taxa = all_metabolites_in_taxa.dropna(subset=chem_id_cols, how='all')

    ### Fill in matching details for compound ID columns from other dataframes
    def fill_match_ids(given_col: str):
        cols_to_do = [cl for cl in chem_id_cols if cl != given_col]
        # Step 1: Create a mapping of cols_to_do[0] and cols_to_do[1] to 'SMILES' values for non-empty rows
        mapping1 = {}
        mapping2 = {}
        for index, row in all_metabolites_in_taxa.iterrows():
            if row[given_col] == row[given_col]:
                if row[cols_to_do[0]] == row[cols_to_do[0]]:
                    mapping1[row[cols_to_do[0]]] = row[given_col]
                if row[cols_to_do[1]] == row[cols_to_do[1]]:
                    mapping2[row[cols_to_do[1]]] = row[given_col]

        if any(x in mapping1.keys() for x in ['', np.nan, 'nan']):
            raise ValueError
        if any(x in mapping2.keys() for x in ['', np.nan, 'nan']):
            raise ValueError

        # Step 2: Iterate through the DataFrame to fill in the empty 'SMILES' values
        for index, row in all_metabolites_in_taxa.iterrows():
            if row[given_col] != row[given_col]:
                m1 = row[cols_to_do[0]]
                m2 = row[cols_to_do[1]]
                # prioritise inchikey and smiles
                if m1 in mapping1:
                    v = mapping1[m1]
                    all_metabolites_in_taxa.at[index, given_col] = v
                elif m2 in mapping2:
                    v = mapping2[m2]
                    all_metabolites_in_taxa.at[index, given_col] = v

    for c_id in chem_id_cols:
        fill_match_ids(c_id)

    all_metabolites_in_taxa['InChIKey_simp'] = all_metabolites_in_taxa['InChIKey'].apply(simplify_inchi_key)
    all_metabolites_in_taxa.to_csv(output_csv)
    return all_metabolites_in_taxa


def add_chembl_data(df: pd.DataFrame, out_csv: str = None, update: bool = False, compound_id_column: str = 'InChiKey'):
    get_compound_info_from_chembl_apm_assays(update=update)
    df_with_assay_data = add_chembl_apm_data_to_compound_df(df, output_csv=out_csv, compound_id_col=compound_id_column)
    return df_with_assay_data


def add_classyfire_info(df: pd.DataFrame, _temp_output_path: str, output_csv: str = None):
    ### Get classyfire info
    # Server was down as of 13/12/23
    all_metabolites_with_classyfire_info = get_classyfire_classes_from_df(df, 'SMILES',
                                                                          tempout_dir=_temp_output_path)

    if output_csv is not None:
        all_metabolites_with_classyfire_info.to_csv(output_csv)
    return all_metabolites_with_classyfire_info


def add_npclassifier_info(df: pd.DataFrame, _temp_output_path: str, output_csv: str = None):
    all_metabolites_with_info = get_npclassifier_classes_from_df(df, 'SMILES',
                                                                 tempout_dir=_temp_output_path)

    if output_csv is not None:
        all_metabolites_with_info.to_csv(output_csv)
    return all_metabolites_with_info


def add_bioavailability_info(df: pd.DataFrame, out_csv: str = None):
    bio_av = add_bioavailability_rules_to_df(df, 'SMILES')

    if out_csv is not None:
        bio_av.to_csv(out_csv)
    return bio_av


def get_manual_files_to_upload(df: pd.DataFrame, _temp_output_path):
    ''' Example for data without, AFAIK, no API. This generates files to upload to obtain information'''

    ### Get NPclassifier info
    # https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22NPCLASSIFIER%22%7D
    # See: https://ccms-ucsd.github.io/GNPSDocumentation/api/#structure-natural-product-classifier-np-classifier
    all_metabolites_to_send_to_npclassifier = df[['SMILES']].dropna().drop_duplicates(
        keep='first')
    all_metabolites_to_send_to_npclassifier.to_csv(
        os.path.join(_temp_output_path, 'smiles_for_np_classifier.csv'))

    ### Get MAIP info
    # https://www.ebi.ac.uk/chembl/maip/
    all_metabolites_to_send_to_maip = df[['SMILES']].dropna().drop_duplicates(
        keep='first')
    all_metabolites_to_send_to_maip['id'] = all_metabolites_to_send_to_maip['SMILES']
    all_metabolites_to_send_to_maip.to_csv(os.path.join(_temp_output_path, 'smiles_for_MAIP.csv'))


def add_manual_info_files(df: pd.DataFrame, npclassifier_output_file: str = None, maip_output_file: str = None):
    all_metabolites_with_info = df.copy(deep=True)
    if npclassifier_output_file is not None:
        ### Get NPclassifier info
        np_classif_results = pd.read_csv(npclassifier_output_file,
                                         sep='\t').drop_duplicates(keep='first')
        rename_dict = {'smiles': 'SMILES'}
        for c in np_classif_results.columns:
            if c != 'smiles':
                rename_dict[c] = 'NPclassif_' + c
        np_classif_results = np_classif_results.rename(columns=rename_dict).dropna(subset='SMILES')
        all_metabolites_with_info = pd.merge(all_metabolites_with_info, np_classif_results[~np_classif_results['SMILES'].isna()], how='left',
                                             on='SMILES')
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


def tidy_and_check_final_dataset(pre_final_df: pd.DataFrame, _temp_output_path: str, final_taxa_compound_csv: str, compound_id_col: str) -> None:
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


if __name__ == '__main__':
    ### Example workflow

    # Define context
    comp_id_column = 'InChIKey'  # Where appropriate, which column should be used to determine compound uniqueness. This is not applicable to some properties, e.g. where SMILES must be used to generate data
    families = ['Pandanaceae']
    wiki_data_id_for_order = 'Q736182'
    temp_outputs_folder = 'temp'
    tidied_outputs_folder = 'tidied'

    ## Get compound-taxa pair data
    # get_wikidata(wiki_data_id_for_order, os.path.join(temp_outputs_folder, 'wikidata.csv'), os.path.join(tidied_outputs_folder, 'wikidata.csv'))
    get_knapsack_data(families, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'knapsack_data.csv'))

    ## Merge and tidy the data
    tidy_wiki_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'wikidata.csv'), index_col=0)
    tidy_knapsack_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'knapsack_data.csv'), index_col=0)
    all_compounds_in_taxa = merge_and_tidy_compound_datasets([tidy_wiki_data, tidy_knapsack_data],
                                                             os.path.join(tidied_outputs_folder, 'merged_data.csv'))

    get_manual_files_to_upload(all_compounds_in_taxa, temp_outputs_folder)

    ## Add extra information related to the compound properties

    # These steps can be included/removed as needed
    # For the longer processes, to avoid repeats you can simply read the associated temp_output if the step has already been run
    with_npclass_classes = add_npclassifier_info(all_compounds_in_taxa, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'npclassifier.csv'))

    with_chembl_data = add_chembl_data(with_npclass_classes, os.path.join(temp_outputs_folder, 'chembl.csv'), compound_id_column=comp_id_column)
    # with_chembl_data = pd.read_csv(os.path.join(temp_outputs_folder, 'chembl.csv'), index_col=0)
    with_bioavailibility = add_bioavailability_info(with_chembl_data, os.path.join(tidied_outputs_folder, 'bioavailibility.csv'))
    # with_bioavailibility = pd.read_csv(os.path.join(tidied_outputs_folder, 'bioavailibility.csv'), index_col=0)
    # Issues with classyfire servers
    # with_classyfire_classes = add_classyfire_info(with_chembl_data, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'classyfire.csv'))

    ## This step requires some manual input

    all_info = add_manual_info_files(with_bioavailibility, maip_output_file=os.path.join(tidied_outputs_folder, 'example_maip_file.csv'))

    ### Then tidy and output final dataset
    tidy_and_check_final_dataset(all_info, tidied_outputs_folder, os.path.join('outputs', 'all_taxa_compound_data.csv'), comp_id_column)
