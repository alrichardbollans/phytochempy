import os
import unittest
import pandas as pd

from phytochempy.chemical_diversity_metrics import add_pathway_information_columns, get_group_level_version_for_all_pathways, \
    get_pathway_based_diversity_measures, NP_PATHWAYS
from phytochempy.compound_properties import get_npclassifier_classes_from_df


class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.df = pd.read_csv(os.path.join('inputs', 'Catalpa_bignonioides_deduplicated_data.csv'), index_col=0)
        self.COMPOUND_ID_COL = 'Standard_SMILES'

        self.TAXON_GROUPING = 'accepted_species'

        # self.data_with_npclass_classes = get_npclassifier_classes_from_df(self.df, 'Standard_SMILES',
        #                                                              'temp_outputs')
        # self.data_with_npclass_classes.to_csv(os.path.join('temp_outputs', 'data_with_npclass_classes.csv'))
        self.data_with_npclass_classes = pd.read_csv(os.path.join('temp_outputs', 'data_with_npclass_classes.csv'), index_col=0)

    def test_something(self):
        group_compound_data_with_ohe_pathways = add_pathway_information_columns(self.data_with_npclass_classes, self.COMPOUND_ID_COL)

        group_distinct_pathway_data = get_group_level_version_for_all_pathways(group_compound_data_with_ohe_pathways,
                                                                               taxon_grouping=self.TAXON_GROUPING,
                                                                               use_distinct=True)

        abundance_diversity = get_pathway_based_diversity_measures(self.data_with_npclass_classes, self.TAXON_GROUPING, self.COMPOUND_ID_COL)

        assert len(abundance_diversity) == len(self.df[self.TAXON_GROUPING].unique().tolist())
        assert len(group_distinct_pathway_data) == len(abundance_diversity)
        assert len(group_compound_data_with_ohe_pathways) == len(self.df)


if __name__ == '__main__':
    unittest.main()
