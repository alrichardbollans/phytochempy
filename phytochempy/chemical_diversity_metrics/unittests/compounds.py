import os
import unittest
import pandas as pd

from phytochempy.chemical_diversity_metrics import add_pathway_information_columns, get_group_level_version_for_all_pathways, \
    get_pathway_based_diversity_measures, NP_PATHWAYS
from phytochempy.chemical_diversity_metrics.compound_distance_metrics import _get_pairwise_distances_from_data, calculate_FAD_measures
from phytochempy.compound_properties import get_npclassifier_classes_from_df


class MyTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_something(self):
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'C(C)O']})
        distmat = _get_pairwise_distances_from_data(df)
        assert distmat == [0]
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO'], 'Groups':['A','A']})
        distmat = _get_pairwise_distances_from_data(df)
        assert distmat[0] > 0

        FADs = calculate_FAD_measures(df, 'Groups')
        assert FADs['FAD'].iloc[0]==distmat*2
        assert FADs['MFAD'].iloc[0]==(distmat*2)/2
        assert FADs['APWD'].iloc[0]==(distmat*2)/4

        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO', 'C(C)N'], 'Groups': ['A', 'A', 'A']})
        distmat = _get_pairwise_distances_from_data(df)
        assert len(distmat) == 3
        assert distmat[0] > 0
        assert distmat[1] > 0
        assert distmat[2] > 0

        FADs = calculate_FAD_measures(df, 'Groups')
        assert FADs['FAD'].iloc[0] == distmat.sum() * 2
        assert FADs['MFAD'].iloc[0] == (distmat.sum() * 2) / 3
        assert FADs['APWD'].iloc[0] == (distmat.sum() * 2) / 9

if __name__ == '__main__':
    unittest.main()
