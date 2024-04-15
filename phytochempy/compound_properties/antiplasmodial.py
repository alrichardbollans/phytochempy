import os

import pandas as pd
import numpy as np
from pkg_resources import resource_filename
from tqdm import tqdm

from phytochempy.compound_properties import simplify_inchi_key, sanitize_filename

_input_path = resource_filename(__name__, 'inputs')
chembl_apm_assay_info_csv = os.path.join(_input_path, 'chembl_apm_assay_info.csv')


def convert_chembl_assay_value_to_ic50(given_val: float):
    """
    Convert a ChEMBL assay value to IC50.

    :param given_val: The ChEMBL assay value to be converted.
    :return: The converted IC50 value.
    """
    if given_val is not None:
        given_val = float(given_val)
        return (10 ** -given_val) * (10 ** 6)
    else:
        return given_val


def update_compound_info_from_chembl_apm_assays(pchembl_active_threshold: float = 6, compound_id_col: str = 'InChIKey_simp',
                                                assay_csv: str = chembl_apm_assay_info_csv):

    # THis needs manually reviewing e.g.
    # https://www.ebi.ac.uk/chembl/g/#browse/activities/full_state/eyJsaXN0Ijp7InNldHRpbmdzX3BhdGgiOiJFU19JTkRFWEVTX05PX01BSU5fU0VBUkNILkFDVElWSVRZIiwiY3VzdG9tX3F1ZXJ5IjoiYXNzYXlfY2hlbWJsX2lkOkNIRU1CTDc2Mjk5MCIsInVzZV9jdXN0b21fcXVlcnkiOnRydWUsInNlYXJjaF90ZXJtIjoiIiwidGV4dF9maWx0ZXIiOiJDSEVNQkwxMTEwNzYifX0%3D
    # Is counted as active, but the ic50 value is Concentration required to reduce chloroquine IC50 by 50%
    compound_data = []
    from chembl_webresource_client.new_client import new_client
    target = new_client.target

    pf = target.filter(pref_name__startswith='Plasmodium ').only('organism', 'target_chembl_id')
    organisms_to_check = [p['organism'] for p in pf]
    print(f'Organisms to check: {organisms_to_check}')
    for p_sp in pf:
        organism = p_sp['organism']
        activity = new_client.activity
        pf_activities = activity.filter(target_chembl_id=p_sp['target_chembl_id'],
                                        pchembl_value__isnull=False,  # has some value
                                        # pchembl_value__gte=6  #See https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
                                        # pChEMBL is defined as: -Log(molar IC50, XC50, EC50, AC50, Ki, Kd or Potency).
                                        # pIC50 value for IC50 < 1μM is 6
                                        # This condition ensures that only compounds with a pIC50 value (the negative logarithm of IC50)
                                        # greater than or equal to 6 (corresponding to IC50 <= 1μM) are retrieved.
                                        ).only('molecule_chembl_id', 'pchembl_value', 'standard_type', 'standard_units', 'standard_value',
                                               'assay_chembl_id')

        # Download the compounds and collect data

        for i in tqdm(range(len(pf_activities)), desc=f'Getting compounds for {organism}', ascii=False, ncols=172):

            compound_activity = pf_activities[i]
            molecule_id = compound_activity['molecule_chembl_id']
            pchembl_value = compound_activity['pchembl_value']
            standard_type = compound_activity['standard_type']
            standard_units = compound_activity['standard_units']
            standard_value = compound_activity['standard_value']
            assay_chembl_id = compound_activity['assay_chembl_id']

            # Get compound details
            compound_details = new_client.molecule.get(molecule_id)
            inchikey = None
            smiles = None
            if compound_details['molecule_structures'] is not None:
                inchikey = compound_details['molecule_structures']['standard_inchi_key']
                smiles = compound_details['molecule_structures']['canonical_smiles']
            name = compound_details['pref_name']
            compound_data.append(
                {'Compound_Name': name, 'assay_standard_value': standard_value, 'assay_standard_units': standard_units,
                 'assay_standard_type': standard_type, 'assay_pchembl_value': pchembl_value, 'target_organism': organism,
                 'assay_chembl_id': assay_chembl_id,
                 'InChIKey': inchikey, 'Smiles': smiles,
                 'molecule_chembl_id': molecule_id, 'natural_product': compound_details['natural_product']})

        # Save the df in case of breaking
        df = pd.DataFrame(compound_data).drop_duplicates(keep='first').reset_index(drop=True)
        df.to_csv(os.path.join(_input_path, f'{sanitize_filename(organism)}.csv'))

    # Create a DataFrame from the compound data
    df = pd.DataFrame(compound_data).drop_duplicates(keep='first').reset_index(drop=True)

    # Add some more info
    df['InChIKey_simp'] = df['InChIKey'].apply(simplify_inchi_key)
    df['assay_ic50_from_pchembl'] = df['assay_pchembl_value'].apply(convert_chembl_assay_value_to_ic50)
    df['mean_ic50'] = df.groupby([compound_id_col])[
        'assay_ic50_from_pchembl'].transform('mean')

    active_chembl_compounds_assays = df[df['assay_pchembl_value'] > pchembl_active_threshold]
    inactive_chembl_compounds_assays = df[df['assay_pchembl_value'] <= pchembl_active_threshold]

    # Compound is active if found in chembl actives, else it is nan unless found in inactives in which case it is inactive.
    df['active_chembl_compound'] = df[compound_id_col].apply(
        lambda x: 1 if x in active_chembl_compounds_assays[compound_id_col].dropna().values else 0 if x in
                                                                                                      inactive_chembl_compounds_assays[
                                                                                                          compound_id_col].dropna().values else np.nan)

    df = df.sort_values(by=compound_id_col)
    df = df.reset_index(drop=True)

    df.to_csv(assay_csv)
    return df


def add_chembl_apm_data_to_compound_df(compound_df: pd.DataFrame, assay_csv: str = chembl_apm_assay_info_csv, output_csv: str = None,
                                       pchembl_active_threshold: float = 6,
                                       compound_id_col: str = 'InChIKey'):
    chembl_compound_assays = pd.read_csv(os.path.join(assay_csv),
                                         index_col=0).dropna(subset=[compound_id_col])
    chembl_compound_assays['mean_ic50_μM_by_id_col'] = chembl_compound_assays.groupby([compound_id_col])[
        'assay_ic50_from_pchembl'].transform('mean')

    active_chembl_compounds = chembl_compound_assays[chembl_compound_assays['assay_pchembl_value'] > pchembl_active_threshold]

    inactive_chembl_compounds = chembl_compound_assays[chembl_compound_assays['assay_pchembl_value'] <= pchembl_active_threshold]

    # Compound is active if found in chembl actives, else it is nan unless found in inactives in which case it is inactive.
    compound_df['active_chembl_compound'] = compound_df[compound_id_col].apply(
        lambda x: 1 if x in active_chembl_compounds[compound_id_col].dropna().values else 0 if x in
                                                                                               inactive_chembl_compounds[
                                                                                                   compound_id_col].dropna().values else np.nan)
    if output_csv is not None:
        compound_df.to_csv(output_csv)
    return compound_df


if __name__ == '__main__':
    update_compound_info_from_chembl_apm_assays()