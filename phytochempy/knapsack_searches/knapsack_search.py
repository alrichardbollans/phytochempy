import os
import time

import pandas as pd
import requests
from tqdm import tqdm
from wcvp_download import wcvp_accepted_columns
from wcvp_name_matching import get_accepted_info_from_names_in_column

from phytochempy.compound_properties import resolve_cas_to_smiles, resolve_cas_to_inchikey, COMPOUND_NAME_COLUMN

_KNAPSACK_organism_column = 'Organism'


def get_knapsack_compounds_for_taxon(name: str):
    """
    Retrieves metabolites for a given taxon.

    :param name: Name of the taxon.
    :return: The metabolite table containing the metabolites for the taxon.
    """
    # Note that this is greedy as name matches in Knapsack search include partial e.g. Cissus matches Narcissus
    url_name_format = name.replace(' ', '%20')

    url = f'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word={url_name_format}'
    try:
        time.sleep(.01)
        tables = pd.read_html(url, flavor='html5lib')

    except UnicodeEncodeError:
        response = requests.get(url)
        decoded = response.content.decode()
        tables = pd.read_html(decoded, flavor='html5lib')
    metabolite_table = tables[0]
    metabolite_table['knapsack_search_term'] = name
    metabolite_table = metabolite_table.rename(columns={'Metabolite': COMPOUND_NAME_COLUMN})
    return metabolite_table


def get_knapsack_compounds_in_family(family: str, temp_output_csv: str):
    """
    :param family: The name of the family to search for metabolites.
    :param temp_output_csv: The path of the temporary output CSV file to save the results.
    :return: None

    This method retrieves metabolite data from the specified family using the wcvp_download and wcvp_columns modules. It searches for all taxa in the given family and retrieves metabolites
    * for each genus. The results are saved in a temporary CSV file.

    Example usage:
    ```
    get_knapsack_compounds_in_family("Rosaceae", "temp_output.csv")
    ```
    """
    if not os.path.isdir(os.path.dirname(temp_output_csv)) and os.path.dirname(temp_output_csv) != '':
        os.mkdir(os.path.dirname(temp_output_csv))

    from wcvp_download import get_all_taxa, wcvp_columns
    # Note that this is greedy as name matches in Knapsack search include partial e.g. Cissus matches Narcissus
    # Account for this by removing results without name resolution
    wcvp_data = get_all_taxa(families_of_interest=[family])

    # Search family by looking for all data from genera
    genera_list = wcvp_data[wcvp_columns['genus']].unique()
    failed_genera = []
    all_genera_df = pd.DataFrame()
    for i in tqdm(range(len(genera_list)), desc=f"Searching genera in Knapsack for {family}…", ascii=False, ncols=80):
        genus = genera_list[i]
        try:
            genus_table = get_knapsack_compounds_for_taxon(genus)
            if len(genus_table.index) > 0:
                all_genera_df = pd.concat([all_genera_df, genus_table])
        except:
            failed_genera.append(genus)

    all_genera_df['Source'] = 'KNApSAcK'

    all_genera_df.to_csv(temp_output_csv)
    if len(failed_genera) > 0:
        print('WARNING: Searching for the following genera raised an error and should be manually checked:')
        print(failed_genera)
        raise ValueError


def tidy_knapsack_results(knapsack_results_csv: str, output_csv: str, family: str, manual_resolution_csv: str = None,
                          wcvp_version_number: str = None, add_smiles_and_inchi: bool = True):
    """
    :param add_smiles_and_inchi:
    :param knapsack_results_csv: The file path of the CSV containing the raw knapsack results.
    :param output_csv: The file path of the CSV to save the tidied results to.
    :param family: The name of the family to filter the results by.
    :param manual_resolution_csv: Optional. The file path of a CSV containing manual resolutions for some compounds.
    :param wcvp_version_number: Optional. The version number of the WCVP database to use for resolving compound names.
    :return: The DataFrame containing the tidied knapsack results.

    """
    all_genera_df = pd.read_csv(knapsack_results_csv, index_col=0)

    if not os.path.isdir(os.path.dirname(output_csv)) and os.path.dirname(output_csv) != '':
        os.mkdir(os.path.dirname(output_csv))
    acc_df = get_accepted_info_from_names_in_column(all_genera_df, _KNAPSACK_organism_column, wcvp_version=wcvp_version_number,
                                                    manual_resolution_csv=manual_resolution_csv)
    acc_df = acc_df[acc_df[wcvp_accepted_columns['family']] == family]
    if add_smiles_and_inchi:
        print('Getting SMILES from CAS ID')
        acc_df['SMILES'] = acc_df['CAS ID'].fillna('').apply(resolve_cas_to_smiles)
        print('Getting InChIKey from CAS ID')
        acc_df['InChIKey'] = acc_df['CAS ID'].fillna('').apply(resolve_cas_to_inchikey)
    acc_df.to_csv(output_csv)

    return acc_df


def get_knapsack_formulas_for_compound(metabolite: str):
    """
    :param metabolite: The name of the metabolite for which to retrieve the knapsack formulas.
    :return: A list of molecular formulas associated with the given metabolite.
    """
    url_name_format = metabolite.replace(' ', '%20')
    url_name_format = url_name_format.replace('+', 'plus')
    url = f'http://www.knapsackfamily.com/knapsack_core/result.php?sname=metabolite&word={url_name_format}'
    try:

        tables = pd.read_html(url, flavor='html5lib')

    except UnicodeEncodeError:
        response = requests.get(url)
        decoded = response.content.decode()
        tables = pd.read_html(decoded, flavor='html5lib')
    meta_table = tables[0]
    try:
        formulas = meta_table['Molecular formula'].values.tolist()
        if formulas == []:
            print(f'Warning: No info found for {metabolite}')
        return formulas
    except (KeyError, IndexError):
        print(f'Warning: No info found for {metabolite}')
        return []


if __name__ == '__main__':
    ## Example usage
    fam = 'Gentianaceae'
    temp_csv = 'knapsack_gentian.csv'
    final_csv = 'knapsack_tidied_gentian.csv'
    get_knapsack_compounds_in_family(fam, temp_csv)
    tidy_knapsack_results(temp_csv, final_csv, fam)