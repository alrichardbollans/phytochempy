import time

import pandas as pd
import requests
from wcvpy.wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns
from wcvpy.wcvp_name_matching import get_accepted_wcvp_info_from_ipni_ids_in_column, output_record_col_names, get_accepted_info_from_names_in_column

from phytochempy.compound_properties import COMPOUND_NAME_COLUMN


def generate_wikidata_search_query(taxon_id: str, limit: int, language: str = 'en') -> str:
    '''Input used query to search all compounds in given taxon: https://query.wikidata.org.
    Example `generate_wikidata_search_query('Q1073514', 10)` will return 10 compounds in Ophiorrhiza
    '''
    query_string = (f'SELECT DISTINCT ?structure ?structureLabel ?structure_smiles ?structure_cas ?structure_inchikey ?organism ?organism_name '
                    f'?ipniID ?chembl_id WHERE {{VALUES ?taxon {{ wd:{taxon_id}}}?organism (wdt:P171*) ?taxon;'
                    f'wdt:P225 ?organism_name.?structure (p:P703/ps:P703) ?organism. OPTIONAL {{?structure wdt:P235 '
                    f'?structure_inchikey.}}OPTIONAL {{?structure wdt:P233 ?structure_smiles.}}OPTIONAL {{?structure wdt:P231 '
                    f'?structure_cas.}}OPTIONAL {{?organism wdt:P961 ?ipniID.}}OPTIONAL {{?structure wdt:P592 ?chembl_id.}}      '
                    f'SERVICE wikibase:label {{bd:serviceParam wikibase:language "{language}".}}}}    LIMIT {str(limit)}')
    print(query_string)
    return query_string


def submit_query(query_string: str, out_csv: str, query_limit: int):
    """
    :param query_string: A string containing the SPARQL query to be submitted.
    :param out_csv: The name of the output CSV file where the query results will be saved.
    :param query_limit: An integer specifying the maximum number of results to retrieve from the query.

    :return: None
    """
    time.sleep(5)  # Rate limiting
    # Define the SPARQL endpoint URL
    endpoint_url = "https://query.wikidata.org/sparql"

    # Prepare the headers and parameters for the POST request
    headers = {"Accept": "application/sparql-results+json"}
    params = {"query": query_string}

    # Send the POST request to the SPARQL endpoint
    response = requests.post(endpoint_url, headers=headers, params=params)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()

        # Extract the results
        results = data['results']['bindings']

        # Convert the results to a DataFrame
        df = pd.json_normalize(results)

        # rename columns to match outputs through csv downloads of manual queries.
        rename_dict = {}
        for c in df.columns.tolist():
            if '.value' in c:
                rename_dict[c] = c.replace('.value', '')
        df = df.rename(columns=rename_dict)

        if len(df) == query_limit:
            print(f'WARNING: Wikidata query returned {len(df)}, which is equal to query limit. To obtain all data, rerun with higher limit.')

        df.to_csv(out_csv)
    else:
        raise ValueError("Wikidata response error. Code:", response.status_code)


def tidy_wikidata_output(wikidata_results_csv: str, output_csv: str, wcvp_version_number: str = None):
    '''

    Resolve names in wikidata query output and tidy.
    
    :return:
    '''

    wiki_df = pd.read_csv(wikidata_results_csv, index_col=0)
    wiki_df = wiki_df[~(
            wiki_df['structure_inchikey'].isna() & wiki_df['structure_cas'].isna() & wiki_df[
        'chembl_id'].isna())]  # only use those structures with some useful identifier
    wiki_df['Source'] = 'WikiData'
    # rename to match other data sources
    wiki_df = wiki_df.rename(
        columns={'structureLabel': COMPOUND_NAME_COLUMN, 'structure_cas': 'CAS ID',
                 'structure_inchikey': 'InChIKey',
                 'structure_smiles': 'SMILES', 'organism_name': 'wikidata_name_snippet',
                 'ipniID': 'wikidata_ipniID'})
    all_taxa = get_all_taxa(version=wcvp_version_number)
    ipni_matches = get_accepted_wcvp_info_from_ipni_ids_in_column(wiki_df,
                                                                  'wikidata_ipniID',
                                                                  all_taxa)[
        wiki_df.columns.tolist() + output_record_col_names]
    unmatched = ipni_matches[ipni_matches[wcvp_columns['wcvp_id']].isna()][wiki_df.columns]

    ipni_matched = ipni_matches[~ipni_matches[wcvp_columns['wcvp_id']].isna()]
    ipni_matched['matched_by'] = 'ipni_id'
    name_matched = get_accepted_info_from_names_in_column(unmatched, 'wikidata_name_snippet', wcvp_version=wcvp_version_number)

    acc_df = pd.concat([ipni_matched, name_matched])
    acc_df = acc_df.dropna(subset=wcvp_accepted_columns['name_w_author'])
    acc_df = acc_df.sort_values(by=wcvp_accepted_columns['name'])
    acc_df.to_csv(output_csv)


def get_wikidata(wiki_data_id: str, temp_output_csv: str, tidied_output_csv: str, limit: int = 100000):
    """
    Fetches data from WikiData based on the provided ID and stores the result in a CSV file after tidying it.

    :param wiki_data_id: The ID of the WikiData entry to fetch data from.
    :param temp_output_csv: The file path to store the temporary output CSV file.
    :param tidied_output_csv: The file path to store the tidied output CSV file.
    :param limit: The maximum number of results to fetch (default is 100000).
    :return: None
    """
    my_query = generate_wikidata_search_query(wiki_data_id, limit)
    submit_query(my_query, temp_output_csv, limit)
    tidy_wikidata_output(temp_output_csv, tidied_output_csv)