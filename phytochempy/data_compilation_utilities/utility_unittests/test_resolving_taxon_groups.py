import unittest
import pandas as pd
import numpy as np
from phytochempy.data_compilation_utilities.resolve_taxon_groupings import get_pathway_version_resolved_at_taxon_level

class TestGetPathwayVersionResolvedAtTaxonLevel(unittest.TestCase):

    def setUp(self):
        self.df = pd.DataFrame({
            'Genus': ['A', 'B', 'C', 'A', 'B', 'A', 'A', 'C'],
            'path': [0, 1, 0, 1, 0, 1, 1, 0],
        })
        self.pathway = 'path'
        self.compound_grouping_col = 'Genus'

    def test_function_runs_without_error(self):
        get_pathway_version_resolved_at_taxon_level(self.df, self.pathway, self.compound_grouping_col)

    def test_return_type(self):
        result = get_pathway_version_resolved_at_taxon_level(self.df, self.pathway, self.compound_grouping_col)
        self.assertIsInstance(result, pd.DataFrame)

    def test_resulting_df_shape(self):
        result = get_pathway_version_resolved_at_taxon_level(self.df, self.pathway, self.compound_grouping_col)
        self.assertEqual(len(result[self.compound_grouping_col].unique()), result.shape[0])

    def test_input_data_unchanged(self):
        df_copy = self.df.copy()
        get_pathway_version_resolved_at_taxon_level(self.df, self.pathway, self.compound_grouping_col)
        pd.testing.assert_frame_equal(self.df, df_copy)

    def test_resulting_df_values(self):
        result = get_pathway_version_resolved_at_taxon_level(self.df, self.pathway, self.compound_grouping_col)
        for genus, group_df in self.df.groupby(self.compound_grouping_col):
            row = result[result[self.compound_grouping_col] == genus].iloc[0]
            self.assertEqual(row['identified_compounds_count'], group_df.shape[0])
            self.assertEqual(row[f'identified_{self.pathway}_count'], group_df[self.pathway].sum())
            self.assertTrue(np.isclose(row[f'mean_identified_as_{self.pathway}'], group_df[self.pathway].mean()))

if __name__ == "__main__":
    unittest.main()
