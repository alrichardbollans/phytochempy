import numpy as np
import pandas as pd


def get_pathway_based_diversity_measures(measure_df: pd.DataFrame, pathways: list, taxon_name_col: str = 'Genus') -> pd.DataFrame:
    """

    This method calculates various diversity measures for pathways based on a given DataFrame.

    This could be refactored to be more user-friendly..

    Parameters:
    - measure_df: A pandas DataFrame containing the data for diversity measures. The DataFrame must have the following columns:
        - taxon_name_col: The names of the taxa, not used in these calculations but kept for output
        - 'mean_identified_as_{pathway}': The mean number of identified compounds for each pathway.
        - 'identified_{pathway}_count': The count of identified compounds for each pathway.
        - 'identified_compounds_count': The total count of identified compounds.

        Note: {pathway} in the column names should be replaced with the actual pathway name.

    - pathways: A list of pathway names to calculate diversity measures for.

    Returns:
    - measure_df: The updated pandas DataFrame containing the calculated diversity measures.

    Note:
    - The following diversity measures are calculated:
        - Shannon index (H)
        - Bias corrected Shannon index (Hbc)
        - Simpson index (G)
        - Pielou index (J)

    - The method also adds min-max scaled versions of the diversity measures to the DataFrame.

    Example usage:
    df = get_pathway_based_diversity_measures(measure_df, pathways)
    """
    ## Read data for all pathways into

    ### Begin with Shannon index
    measure_df['H'] = 0
    for pathway in pathways:
        measure_df[f'ln_mean_identified_as_{pathway}'] = np.log(measure_df[f'mean_identified_as_{pathway}']).replace(-np.inf, 0)
        measure_df['H'] = measure_df['H'] + measure_df[f'mean_identified_as_{pathway}'] * measure_df[
            f'ln_mean_identified_as_{pathway}']

    measure_df['H'] = -measure_df['H']
    ## Bias corrected shannon
    # From chao_nonparametric_2003, following beck_comparing_2010.
    # Note that there are updated metrics for calculating coverage e.g. chao_coveragebased_2012
    measure_df['number_singletons'] = 0
    for pathway in pathways:
        measure_df.loc[measure_df[f'identified_{pathway}_count'] == 1, 'number_singletons'] += 1
    measure_df['sample_coverage'] = 1 - (measure_df['number_singletons'] / measure_df['identified_compounds_count'])

    measure_df['Hbc'] = 0
    for pathway in pathways:
        ###New log
        measure_df[f'ln_Cmean_identified_as_{pathway}'] = np.log(
            measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'])
        addition = ((measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'] * measure_df[f'ln_Cmean_identified_as_{pathway}']) /
                    (1 - (1 - measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage']) ** measure_df[
                        'identified_compounds_count'])).fillna(0)
        measure_df['Hbc'] = measure_df['Hbc'] + addition
    measure_df['Hbc'] = -measure_df['Hbc']

    # Simpson index also refered to as gini-simpson
    # Used in e.g. corre_evaluation_2023
    measure_df['G'] = 0
    for pathway in pathways:
        measure_df['G'] = measure_df['G'] + (
                measure_df[f'mean_identified_as_{pathway}'] * measure_df[f'mean_identified_as_{pathway}'])
    measure_df['G'] = 1 - measure_df['G']

    measure_df['number_of_apparent_categories'] = 0
    for pathway in pathways:
        measure_df[f'binary_identified_as_{pathway}'] = 0
        measure_df.loc[measure_df[f'identified_{pathway}_count'] > 0, f'binary_identified_as_{pathway}'] = 1
        measure_df['number_of_apparent_categories'] = measure_df['number_of_apparent_categories'] + measure_df[
            f'binary_identified_as_{pathway}']

    ## Pielou index measures evenness of the classes
    ## It normalises H by the 'richness' for the given genus, i.e. the number of different pathways present
    ## as discussed in corre_evaluation_2023.
    ## Where number_of_apparent_categories =1, this is left undefined
    measure_df['J'] = measure_df['H'] / (np.log(measure_df['number_of_apparent_categories']))

    measure_df = measure_df[[taxon_name_col, 'H', 'Hbc', 'G', 'J']]
    from sklearn.preprocessing import MinMaxScaler
    for index in ['H', 'Hbc', 'G', 'J']:
        scaler = MinMaxScaler()
        measure_df[index + '_minmax'] = scaler.fit_transform(measure_df[[index]])
        print(index)

    return measure_df
