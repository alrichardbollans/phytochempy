import pandas as pd


def rarefy_diversity_for_group(df: pd.DataFrame, compound_grouping: str, group: str, target_size: int, metrics: list, method_to_calculate,
                               iterations=1000,
                               compound_id_col: str = None):
    diversity_calc_means = {}
    diversity_calc_stds = {}
    group_df = df[df[compound_grouping] == group]
    bootstrap_df = pd.DataFrame()
    for _ in range(iterations):
        # Randomly subsample to the target size
        subsample = group_df.sample(target_size, replace=False)

        # Calculate pairwise distances for the subsample
        if compound_id_col is None:
            distances = method_to_calculate(subsample, compound_grouping)
        else:
            distances = method_to_calculate(subsample, compound_grouping, compound_id_col)  # pdist returns condensed distance matrix

        assert len(distances) == 1
        bootstrap_df = pd.concat([bootstrap_df, distances])

    assert len(bootstrap_df) == iterations

    for m in metrics:
        mean_distance = bootstrap_df[m].mean()

        diversity_calc_means[f'{m}_Rare'] = mean_distance
        diversity_calc_stds[f'{m}_Rare'] = bootstrap_df[m].std()
    diversity_calc_means[compound_grouping] = group
    diversity_calc_stds[compound_grouping] = group
    # Return the mean and standard deviation of diversity calcs over iterations
    return diversity_calc_means, diversity_calc_stds
