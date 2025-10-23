import os.path
import unittest

import pandas as pd
from phytochempy.wikidata_searches import tidy_wikidata_output, generate_wikidata_search_query, submit_query, get_wikidata_id_for_taxon


class TestTidyWikidataOutput(unittest.TestCase):

    def setUp(self):
        self.manual_test_query = os.path.join('inputs', 'manual_query.csv')
        self.correct_output = os.path.join('inputs', 'correct_output.csv')

    def test_tidy_wikidata_output(self):
        # Arrange

        output_csv = 'output.csv'

        # Act
        tidy_wikidata_output(self.manual_test_query, output_csv)

        output_df = pd.read_csv(output_csv, index_col=0)
        correct_df = pd.read_csv(self.correct_output, index_col=0)

        pd.testing.assert_frame_equal(output_df[correct_df.columns], correct_df)

    def test_querying(self):
        my_query = generate_wikidata_search_query('Q1073514', 10)
        submit_query(my_query, 'wikidata_search.csv', 10)

        output_df = pd.read_csv('wikidata_search.csv', index_col=0)
        correct_df = pd.read_csv(self.manual_test_query, index_col=0)
        output_df = output_df[correct_df.columns]
        pd.testing.assert_frame_equal(output_df, correct_df)


class TestGettingID(unittest.TestCase):
    def test_examples(self):
        examples = {'Google': [], 'Pandanaceae': ['Q736182'], 'Gentianales': ['Q21754'], 'Cactaceae':['Q14560']}
        for e in examples:
            self.assertEqual(get_wikidata_id_for_taxon(e), examples[e])

    def test_bad_examples(self):

        try:
            get_wikidata_id_for_taxon(None)
        except TypeError:
            print('passed')
        self.assertEqual(get_wikidata_id_for_taxon('None)'), [])


if __name__ == "__main__":
    unittest.main()
