import math

import pandas as pd


def get_pathway_version_resolved_at_taxon_level(df: pd.DataFrame, pathway: str, compound_grouping_col: str):
    """
    :param df: A pandas DataFrame containing the data
    :param pathway: A string specifying the column name of the class
    :param compound_grouping_col: A string specifying the taxon grouping column name (default is 'Genus')
    :return: A pandas DataFrame containing genus level version for pathway

    Given a dataset of organism-compound pairs with a binary column 'pathway', resolves the data to a higher taxonomic level based on categories
    (e.g. genus names) given in the 'compound_grouping_col' column.

    Returns a dataframe with mean values of the pathway for the new taxonomic groups, as well as normalised means and other related metrics.

    The normalisation of means follows:
    Daniele Micci-Barreca, ‘A Preprocessing Scheme for High-Cardinality Categorical Attributes in Classification and Prediction Problems’,
    ACM SIGKDD Explorations Newsletter 3, no. 1 (July 2001): 27–32, https://doi.org/10.1145/507533.507538.

    The means for each class are highly unreliable for small counts so we use a blend of posterior and prior probabilties to improve this i.e. low evidence examples are corrected towards the population mean.
    In below, k (as in original paper) determines half of the sample size for which we trust the mean estimate
    f denotes the smoothing effect to balance categorical average vs prior. Higher value means stronger regularization.
    Here we use the defaults used in the target encoder library

    """
    expected_mean = df[pathway].mean()

    genera_df = df.copy()

    num_genera_tested = len(genera_df[compound_grouping_col].unique().tolist())

    # count labelled species
    counts = genera_df.value_counts(compound_grouping_col)
    counts.name = 'identified_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on=compound_grouping_col, right_index=True)
    N_class_col = f'identified_{pathway}_count'
    genera_df[N_class_col] = genera_df.groupby([compound_grouping_col])[pathway].transform('sum')
    mean_col = f'mean_identified_as_{pathway}'
    genera_df[mean_col] = genera_df.groupby([compound_grouping_col])[pathway].transform('mean')
    expected_mean_col = f'expected_total_mean_for_{pathway}'
    genera_df[expected_mean_col] = expected_mean

    # Normalised mean for some analysis following:
    # Daniele Micci-Barreca, ‘A Preprocessing Scheme for High-Cardinality Categorical Attributes in Classification and Prediction Problems’,
    # ACM SIGKDD Explorations Newsletter 3, no. 1 (July 2001): 27–32, https://doi.org/10.1145/507533.507538.
    # The means for each class are highly unreliable for small counts
    # We can use a blend of posterior and prior probabilties to improve this i.e. low evidence examples are corrected towards the population mean.
    # In below, k (as in original paper) determines half of the sample size for which we trust the mean estimate
    # f denotes the smoothing effect to balance categorical average vs prior. Higher value means stronger regularization.
    # Here we use the defaults used in the target encoder library
    def weighting_factor(given_val: float, k: int = 20, f: float = 10) -> float:
        denom = 1 + (math.e ** ((k - given_val) / f))
        return 1 / denom

    factor_col = f'weigthing_factor_for_{pathway}'
    genera_df[factor_col] = genera_df['identified_compounds_count'].apply(weighting_factor)

    norm_mean_col = f'norm_mean_identified_as_{pathway}'
    genera_df[norm_mean_col] = (genera_df[factor_col] * genera_df[mean_col]) + (
            1 - genera_df[factor_col]) * expected_mean

    genera_df = genera_df[
        [compound_grouping_col, 'identified_compounds_count', N_class_col, mean_col, expected_mean_col, factor_col,
         norm_mean_col]]
    genera_df = genera_df.reset_index(drop=True)

    genera_df = genera_df.drop_duplicates(subset=[compound_grouping_col])
    genera_df.reset_index(drop=True, inplace=True)
    assert len(genera_df) == num_genera_tested
    return genera_df
