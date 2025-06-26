import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

from phytochempy.chemical_diversity_metrics import calculate_FAD_measures, get_pathway_based_diversity_measures


def rarefy_diversity_for_group(df: pd.DataFrame, compound_grouping: str, group: str, target_size: int, metrics: list, method_to_calculate,
                               iterations=1000,
                               compound_id_col: str = None):
    diversity_calc_means = {}
    diversity_calc_stds = {}
    group_df = df[df[compound_grouping] == group]
    bootstrap_results = []
    for _ in range(iterations):
        # Randomly subsample to the target size
        subsample = group_df.sample(target_size, replace=False)

        # Calculate pairwise distances for the subsample
        if compound_id_col is None:
            distances = method_to_calculate(subsample, compound_grouping)
        else:
            distances = method_to_calculate(subsample, compound_grouping, compound_id_col)

        if compound_id_col is None:
            assert len(distances) == 1
        else:
            assert len(distances) <= 1

        bootstrap_results.append(distances)
    bootstrap_df = pd.concat(bootstrap_results)

    if compound_id_col is None:
        assert len(bootstrap_df) == iterations

    for m in metrics:
        mean_distance = bootstrap_df[m].mean()

        diversity_calc_means[f'{m}_Rare'] = mean_distance
        diversity_calc_stds[f'{m}_Rare'] = bootstrap_df[m].std()
    diversity_calc_means[compound_grouping] = group
    diversity_calc_stds[compound_grouping] = group
    # Return the mean and standard deviation of diversity calcs over iterations
    return diversity_calc_means, diversity_calc_stds

def _process_group(args):
    (df, compound_grouping, group, target_size, iterations, compound_id_col) = args
    fad_means, fad_stds = rarefy_diversity_for_group(df, compound_grouping, group, target_size, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures,
                                                     iterations=iterations)
    group_df = pd.DataFrame(fad_means, index=[group])

    pathway_means, pathway_stds = rarefy_diversity_for_group(df, compound_grouping, group, target_size, ['H', 'Hbc', 'G', 'J'],
                                                             get_pathway_based_diversity_measures,
                                                             iterations=iterations, compound_id_col=compound_id_col)
    pathway_group_df = pd.DataFrame(pathway_means, index=[group])
    return (group_df, pathway_group_df)


def compile_rarified_calculations(df: pd.DataFrame, compound_grouping: str, target_size: int, compound_id_col: str, iterations=1000):
    fad_all_groups = []
    pathway_all_groups = []
    groups = df[compound_grouping].unique()
    num_cpus = max(multiprocessing.cpu_count() - 1, 1)
    args = [
        (df, compound_grouping, group, target_size, iterations, compound_id_col)
        for group in groups
    ]

    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(_process_group, a) for a in args]
        for future in as_completed(futures):
            group_df, pathway_group_df = future.result()
            fad_all_groups.append(group_df)
            pathway_all_groups.append(pathway_group_df)
    fad_all_groups = pd.concat(fad_all_groups)
    pathway_all_groups = pd.concat(pathway_all_groups)
    return fad_all_groups, pathway_all_groups
