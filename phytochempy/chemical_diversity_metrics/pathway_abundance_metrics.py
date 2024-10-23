import numpy as np
import pandas as pd
from pandas import DataFrame
from typing import Any

from phytochempy.compound_properties import get_npclassifier_pathway_columns_in_df
from phytochempy.data_compilation_utilities import get_pathway_version_resolved_at_taxon_level

NP_PATHWAYS = ['Terpenoids', 'Fatty_acids', 'Polyketides', 'Carbohydrates', 'Amino_acids_and_Peptides', 'Shikimates_and_Phenylpropanoids',
               'Alkaloids']


def split_multiple_pathways_into_duplicate_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Resolve cases with multiple assinged compounds by separating into different rows

    :param df: A pandas DataFrame containing multiple pathways encoded as binary values.
    :return: A pandas DataFrame with duplicate rows created for each pathway that has multiple occurrences.
    """
    issues_to_resolve = df[df[NP_PATHWAYS].sum(axis=1) > 1]
    if len(issues_to_resolve) > 0:
        # List to store new rows
        new_rows = []

        # Iterate through each row in the dataframe
        for index, row in issues_to_resolve.iterrows():
            # Count the number of 1s in the row
            count_ones = row[NP_PATHWAYS].sum()

            # If more than one 1 is found
            if count_ones > 1:
                one_cols = []
                for col in NP_PATHWAYS:
                    if row[col] == 1:
                        one_cols.append(col)
                # Iterate through each column in the row
                for col in one_cols:

                    new_row = row.copy()
                    for p in NP_PATHWAYS:
                        new_row[p] = 0
                    new_row[col] = 1
                    new_rows.append(new_row.T)

        # Add new rows to dataframe
        resolution_df = pd.DataFrame()._append(new_rows)
        assert len(resolution_df[resolution_df[NP_PATHWAYS].sum(axis=1) > 1]) == 0
        # Drop duplicate rows
        df = df[df[NP_PATHWAYS].sum(axis=1) < 2]

        out_df = pd.concat([df, resolution_df])

        return out_df
    else:
        return df


def get_group_level_version_for_all_pathways(df: pd.DataFrame, taxon_grouping, use_distinct: bool = False) -> pd.DataFrame:
    '''## Generate group data for all pathways

    ## 'Distinct' Pathway data splits compounds into multiple rows if they are associated with multiple compounds
    ## This is only used for calculations of diversity indices
    ## Else 'nondistinct' pathway data should be used.

    - returns A pandas DataFrame containing the data for diversity measures. The DataFrame has the following columns:
        - taxon_grouping: The names of the groups, not used in these calculations but kept for output
        - 'mean_identified_as_{pathway}': The mean number of identified compounds for each pathway.
        - 'identified_{pathway}_count': The count of identified compounds for each pathway.
        - 'identified_compounds_count': The total count of identified compounds.

    '''

    if use_distinct:
        new_df = split_multiple_pathways_into_duplicate_rows(df)
    else:
        new_df = df.copy()

    out_df = pd.DataFrame()
    out_df[taxon_grouping] = new_df[taxon_grouping].unique()
    original_length = len(out_df)

    for pathway in NP_PATHWAYS:
        group_pathway_df = get_pathway_version_resolved_at_taxon_level(new_df, pathway, taxon_grouping_col=taxon_grouping)
        if 'identified_compounds_count' not in out_df.columns:
            out_df = pd.merge(out_df, group_pathway_df, on=[taxon_grouping], how='left')
        else:
            out_df = pd.merge(out_df, group_pathway_df, on=[taxon_grouping, 'identified_compounds_count'], how='left')
    assert len(out_df) == original_length

    return out_df


def separate_into_pathway(df: pd.DataFrame, pathway: str) -> tuple[DataFrame, Any]:
    """
    Separate the given DataFrame into two separate DataFrames based on a specified pathway.

    :param df: The input DataFrame.
    :param pathway: The pathway to separate the DataFrame by.
    :return: A tuple containing two separate DataFrames: positives and negatives.
    """
    df = df.dropna(subset=['NPclassif_pathway_results'])
    pathway_cols = get_npclassifier_pathway_columns_in_df(df)

    positives = pd.DataFrame()
    for col in pathway_cols:
        positives = pd.concat([positives, df[df[col] == pathway]])
    negatives = df[~df.index.isin(positives.index)]
    assert len(positives) + len(negatives) == len(df)
    problems = positives[~positives['NPclassif_pathway_results'].str.contains(pathway)]
    assert len(problems) == 0
    problems = negatives[negatives['NPclassif_pathway_results'].str.contains(pathway)]
    assert len(problems) == 0

    return positives, negatives


def add_pathway_information_columns(df: pd.DataFrame, compound_id_col: str) -> pd.DataFrame:
    """
    Add pathway information columns to a DataFrame.

    :param df: The DataFrame to add pathway information columns to.
    :return: The DataFrame with pathway information columns added.
    """
    df = df.dropna(subset=['NPclassif_pathway_results'])
    original_length = len(df)

    for pathway in NP_PATHWAYS:
        relevant_paths, other_paths = separate_into_pathway(df, pathway)
        relevant_paths[pathway] = 1
        other_paths[pathway] = 0

        pathway_df = pd.concat([relevant_paths, other_paths])[[compound_id_col, pathway]]
        pathway_df = pathway_df.drop_duplicates(subset=[compound_id_col, pathway])
        amibiguous_duplicates = pathway_df[pathway_df[compound_id_col].duplicated(keep=False)]
        if len(amibiguous_duplicates) > 0:  # A check that comps with same ID have been assigned same class
            print(
                f'WARNING: Some ambiguity for pathway: {pathway}. This is likely due to differing smiles strings for same given {compound_id_col}.')
            print(amibiguous_duplicates)

            raise ValueError

        df = df.merge(pathway_df, how='left', on=compound_id_col)
    assert len(df) == original_length  # Check no additions from merge

    return df


def get_pathway_based_diversity_measures(df: pd.DataFrame, taxon_grouping: str, compound_id_col: str) -> pd.DataFrame:
    """

    This method calculates various diversity measures for pathways based on a given DataFrame.

    This could be refactored to be more user-friendly.

    Parameters:
    - taxon_grouping: the column used to group compounds
    - compound_id_col: the column id column, to determine compoun uniqueness

    Returns:
    - measure_df: The updated pandas DataFrame containing the calculated diversity measures.

    Note:
    - The following diversity measures are calculated:
        - Shannon index (H)
        - Bias corrected Shannon index (Hbc)
        - Simpson index (G)
        - Pielou index (J)

    """

    group_compound_data_with_ohe_pathways = add_pathway_information_columns(df, compound_id_col)

    measure_df = get_group_level_version_for_all_pathways(group_compound_data_with_ohe_pathways, taxon_grouping=taxon_grouping,
                                                          use_distinct=True)

    ### Begin with Shannon index
    measure_df['H'] = 0
    for pathway in NP_PATHWAYS:
        measure_df[f'ln_mean_identified_as_{pathway}'] = np.log(measure_df[f'mean_identified_as_{pathway}']).replace(-np.inf, 0)
        measure_df['H'] = measure_df['H'] + measure_df[f'mean_identified_as_{pathway}'] * measure_df[
            f'ln_mean_identified_as_{pathway}']

    measure_df['H'] = -measure_df['H']
    ## Bias corrected shannon
    # From chao_nonparametric_2003, following beck_comparing_2010.
    # Note that there are updated metrics for calculating coverage e.g. chao_coveragebased_2012
    measure_df['number_singletons'] = 0
    for pathway in NP_PATHWAYS:
        measure_df.loc[measure_df[f'identified_{pathway}_count'] == 1, 'number_singletons'] += 1
    measure_df['sample_coverage'] = 1 - (measure_df['number_singletons'] / measure_df['identified_compounds_count'])

    measure_df['Hbc'] = 0
    for pathway in NP_PATHWAYS:
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
    for pathway in NP_PATHWAYS:
        measure_df['G'] = measure_df['G'] + (
                measure_df[f'mean_identified_as_{pathway}'] * measure_df[f'mean_identified_as_{pathway}'])
    measure_df['G'] = 1 - measure_df['G']

    measure_df['number_of_apparent_categories'] = 0
    for pathway in NP_PATHWAYS:
        measure_df[f'binary_identified_as_{pathway}'] = 0
        measure_df.loc[measure_df[f'identified_{pathway}_count'] > 0, f'binary_identified_as_{pathway}'] = 1
        measure_df['number_of_apparent_categories'] = measure_df['number_of_apparent_categories'] + measure_df[
            f'binary_identified_as_{pathway}']

    ## Pielou index measures evenness of the classes
    ## It normalises H by the 'richness' for the given genus, i.e. the number of different pathways present
    ## as discussed in corre_evaluation_2023.
    ## Where number_of_apparent_categories =1, this is left undefined
    measure_df['J'] = measure_df['H'] / (np.log(measure_df['number_of_apparent_categories']))

    measure_df = measure_df[[taxon_grouping, 'H', 'Hbc', 'G', 'J']]

    return measure_df
