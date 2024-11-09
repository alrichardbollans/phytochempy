import os
import unittest
from operator import index

import pandas as pd

from phytochempy.chemical_diversity_metrics import add_pathway_information_columns, get_group_level_version_for_all_pathways, \
    get_pathway_based_diversity_measures, NP_PATHWAYS, rarefy_diversity_for_group, compile_rarified_calculations
from phytochempy.chemical_diversity_metrics.compound_distance_metrics import _get_pairwise_distances_from_data, calculate_FAD_measures
from phytochempy.compound_properties import get_npclassifier_classes_from_df


class MyTestCase(unittest.TestCase):

    def setUp(self):
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO', 'C(C)N', 'C', 'CO'], 'Groups': ['A', 'A', 'A', 'B', 'B']})
        self.df = get_npclassifier_classes_from_df(df, 'Standard_SMILES', 'inputs')

    def test_multiple(self):
        rarified_A = rarefy_diversity_for_group(self.df, 'Groups', 'A', 2, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures)

        FADs = calculate_FAD_measures(self.df, 'Groups')
        AFADS = FADs[FADs['Groups'] == 'A']
        BFADS = FADs[FADs['Groups'] == 'B']

        assert round(rarified_A[0]['APWD_Rare'], 2) == round(AFADS['APWD'].iloc[0], 2)  # APWD is similar when target size is similar to group size

        rarified_B = rarefy_diversity_for_group(self.df, 'Groups', 'B', 2, ['FAD', 'MFAD', 'APWD'], calculate_FAD_measures)
        assert round(rarified_B[0]['APWD_Rare'], 2) == round(BFADS['APWD'].iloc[0], 2)  # APWD is similar when target size is similar to group size
        assert round(rarified_B[0]['FAD_Rare'], 2) == round(BFADS['FAD'].iloc[0], 2)  # FAD is similar when target size is same as group size

        print(AFADS)

    def test_pways(self):

        rarified_A = rarefy_diversity_for_group(self.df, 'Groups', 'A', 3, ['H', 'Hbc', 'G', 'J'], get_pathway_based_diversity_measures,
                                                compound_id_col='Standard_SMILES', iterations=10)

        FADs = get_pathway_based_diversity_measures(self.df, 'Groups','Standard_SMILES')
        AFADS = FADs[FADs['Groups'] == 'A']
        BFADS = FADs[FADs['Groups'] == 'B']

        assert round(rarified_A[0]['H_Rare'], 2) == round(AFADS['H'].iloc[0], 2)


        # H is similar when target size is similar to group size
        try:
            rarified_B = rarefy_diversity_for_group(self.df, 'Groups', 'B', 3, ['H', 'Hbc', 'G', 'J'], get_pathway_based_diversity_measures,
                                                    compound_id_col='Standard_SMILES', iterations=10)
        except ValueError:
            print('ok')
        else:
            raise ValueError

    def test_compiling(self):
        fad_all_groups, pathway_all_groups = compile_rarified_calculations(self.df, 'Groups',2, 'Standard_SMILES', 10)

        assert len(fad_all_groups) == 2
        assert len(pathway_all_groups) == 2
        print(fad_all_groups)


if __name__ == '__main__':
    unittest.main()
