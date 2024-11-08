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

    def test_zero(self):
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'C(C)O']})
        distmat = _get_pairwise_distances_from_data(df)
        assert distmat == [0]

    def test_one(self):

        df1 = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO'], 'Groups':['A','A']})
        distmat = _get_pairwise_distances_from_data(df1)
        assert distmat[0] > 0

        FADs = calculate_FAD_measures(df1, 'Groups')
        assert FADs['FAD'].iloc[0]==distmat*2
        assert FADs['MFAD'].iloc[0]==(distmat*2)/2
        assert FADs['APWD'].iloc[0]==(distmat*2)/2
    def test_two(self):
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO', 'C(C)N'], 'Groups': ['A', 'A', 'A']})
        distmat = _get_pairwise_distances_from_data(df)
        assert len(distmat) == 3
        assert distmat[0] > 0
        assert distmat[1] > 0
        assert distmat[2] > 0

        FADs = calculate_FAD_measures(df, 'Groups')
        assert FADs['FAD'].iloc[0] == distmat.sum() * 2
        assert FADs['MFAD'].iloc[0] == (distmat.sum() * 2) / 3
        assert FADs['APWD'].iloc[0] == (distmat.sum() * 2) / 6

    def test_multiple(self):
        df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'CO', 'C(C)N', 'C', 'CO'], 'Groups': ['A', 'A', 'A', 'B', 'B']})
        distmatA = _get_pairwise_distances_from_data(df[df['Groups'] == 'A'])
        distmatB = _get_pairwise_distances_from_data(df[df['Groups'] == 'B'])
        assert len(distmatA) == 3
        assert distmatA[0] > 0
        assert distmatA[1] > 0
        assert distmatA[2] > 0

        assert len(distmatB) == 1
        assert distmatB[0] > 0

        FADs = calculate_FAD_measures(df, 'Groups')
        AFADS = FADs[FADs['Groups'] == 'A']
        BFADS = FADs[FADs['Groups'] == 'B']
        assert AFADS['FAD'].iloc[0] == distmatA.sum() * 2
        assert BFADS['FAD'].iloc[0] == distmatB.sum() * 2


        assert AFADS['MFAD'].iloc[0] == (distmatA.sum() * 2) / 3
        assert BFADS['MFAD'].iloc[0] == (distmatB.sum() * 2) / 2


        assert AFADS['APWD'].iloc[0] == (distmatA.sum() * 2) / 6
        assert BFADS['APWD'].iloc[0] == (distmatB.sum() * 2) / 2

    def test_examples(self):
        # from Comastoma tenellum
        df = pd.DataFrame({'Standard_SMILES': ['COc1cc(O)c2c(=O)c3c(O)c(OC)ccc3oc2c1', 'COc1cc(O)c2c(=O)c3c(OC)c(OC)ccc3oc2c1'], 'Groups': ['Comastoma tenellum', 'Comastoma tenellum']})
        distmat = _get_pairwise_distances_from_data(df)

        FADs = calculate_FAD_measures(df, 'Groups')
        print(FADs)
if __name__ == '__main__':
    unittest.main()
