import os
import unittest
from operator import index

import pandas as pd

from phytochempy.chemical_diversity_metrics import add_pathway_information_columns, get_group_level_version_for_all_pathways, \
    get_pathway_based_diversity_measures, NP_PATHWAYS, rarefy_diversity_for_group
from phytochempy.chemical_diversity_metrics.compound_distance_metrics import _get_pairwise_distances_from_data, calculate_FAD_measures
from phytochempy.compound_properties import get_npclassifier_classes_from_df


class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO', 'C(C)N', 'C', 'CO'], 'Groups': ['A', 'A', 'A', 'B', 'B']})

    def test_multiple(self):
        rarified_A = rarefy_diversity_for_group(self.df, 'Groups', 'A', 2, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures)

        FADs = calculate_FAD_measures(self.df, 'Groups')
        AFADS = FADs[FADs['Groups'] == 'A']
        BFADS = FADs[FADs['Groups'] == 'B']

        assert round(rarified_A[0]['APWD'], 2) == round(AFADS['APWD'].iloc[0], 2)  # APWD is similar when target size is similar to group size

        rarified_B = rarefy_diversity_for_group(self.df, 'Groups', 'B', 2, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures)
        assert round(rarified_B[0]['APWD'], 2) == round(BFADS['APWD'].iloc[0], 2)  # APWD is similar when target size is similar to group size
        assert round(rarified_B[0]['FAD'], 2) == round(BFADS['FAD'].iloc[0], 2)  # FAD is similar when target size is same as group size

        print(AFADS)
    def test_compiling(self):
        all_groups = pd.DataFrame()
        for group in self.df['Groups'].unique():
            fad_means, fad_stds = rarefy_diversity_for_group(self.df,  'Groups', group, 2, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures)
            group_df = pd.DataFrame(fad_means, index=[group])
            all_groups = pd.concat([all_groups, group_df])
            print(fad_means)

if __name__ == '__main__':
    unittest.main()
