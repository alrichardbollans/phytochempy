import os.path
import unittest

import numpy as np
import pandas as pd

from phytochempy.compound_properties import npclassify_smiles, get_npclassif_classes_from_smiles, get_npclassifier_classes_from_df, \
    read_manual_npclassifier_input


class TestNPClassifierMethods(unittest.TestCase):

    def setUp(self):
        self.smiles_list = {'CC(=C)CCO': {'NPclassif_class_results': 'Acyclic monoterpenoids:Fatty alcohols',
                                          'NPclassif_class_results_0': 'Acyclic monoterpenoids',
                                          'NPclassif_class_results_1': 'Fatty alcohols',
                                          'NPclassif_isglycoside': False,
                                          'NPclassif_pathway_results': 'Fatty acids:Terpenoids',
                                          'NPclassif_pathway_results_0': 'Fatty acids',
                                          'NPclassif_pathway_results_1': 'Terpenoids',
                                          'NPclassif_superclass_results': 'Fatty acyls:Monoterpenoids',
                                          'NPclassif_superclass_results_0': 'Fatty acyls',
                                          'NPclassif_superclass_results_1': 'Monoterpenoids',
                                          'SMILES': 'CC(=C)CCO'},
                            'CCCCCCCCCCCCCCCC(=O)O': {'NPclassif_class_results': 'Branched fatty acids:Unsaturated fatty acids',
                                                      'NPclassif_class_results_0': 'Branched fatty acids',
                                                      'NPclassif_class_results_1': 'Unsaturated fatty acids',
                                                      'NPclassif_isglycoside': False,
                                                      'NPclassif_pathway_results': 'Fatty acids',
                                                      'NPclassif_pathway_results_0': 'Fatty acids',
                                                      'NPclassif_superclass_results': 'Fatty Acids and Conjugates',
                                                      'NPclassif_superclass_results_0': 'Fatty Acids and Conjugates',
                                                      'SMILES': 'CCCCCCCCCCCCCCCC(=O)O'},
                            "CC(=O)OC1=CC=CC=C1C(=O)O": {'SMILES': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                                                         'NPclassif_class_results': 'Simple phenolic acids',
                                                         'NPclassif_class_results_0': 'Simple phenolic acids',
                                                         'NPclassif_isglycoside': False,
                                                         'NPclassif_pathway_results': 'Shikimates and Phenylpropanoids',
                                                         'NPclassif_pathway_results_0': 'Shikimates and Phenylpropanoids',
                                                         'NPclassif_superclass_results': 'Phenolic acids (C6-C1)',
                                                         'NPclassif_superclass_results_0': 'Phenolic acids (C6-C1)'},
                            "C1=CC=CC=C1": {'SMILES': 'C1=CC=CC=C1',
                                            'NPclassif_class_results': np.nan, 'NPclassif_isglycoside': False,
                                            'NPclassif_pathway_results': 'Shikimates and Phenylpropanoids',
                                            'NPclassif_pathway_results_0': 'Shikimates and Phenylpropanoids',
                                            'NPclassif_superclass_results': np.nan},
                            "C1CCCCC1": {'SMILES': 'C1CCCCC1',
                                         'NPclassif_class_results': 'Lactones',
                                         'NPclassif_class_results_0': 'Lactones',
                                         'NPclassif_isglycoside': False,
                                         'NPclassif_pathway_results': 'Fatty acids',
                                         'NPclassif_pathway_results_0': 'Fatty acids',
                                         'NPclassif_superclass_results': 'Fatty esters',
                                         'NPclassif_superclass_results_0': 'Fatty esters'}}

    def test_npclassify_smiles(self):
        # Classify all SMILES strings in the list
        for smiles in self.smiles_list:
            classification_result = npclassify_smiles(smiles)
            self.assertEqual(classification_result, self.smiles_list[smiles])

    def test_nulls(self):
        self.assertIsNone(npclassify_smiles(None))
        self.assertIsNone(npclassify_smiles(np.nan))

    def test_list(self):
        result = get_npclassif_classes_from_smiles(list(self.smiles_list.keys()), 'temp_outputs')
        classes = result.columns.tolist()
        classes.remove('SMILES')
        classes.remove('NPclassif_isglycoside')
        for sm in self.smiles_list:
            this_df = result[result['SMILES'] == sm]
            self.assertEqual(len(this_df), 1)
            for c in classes:
                if c in self.smiles_list[sm]:
                    print(sm)
                    print(c)
                    correct = self.smiles_list[sm][c]
                    test = this_df[c].values[0]
                    if correct == correct:
                        self.assertEqual(correct, test)
                    else:
                        self.assertTrue(np.isnan(test))

    def test_df(self):

        df = pd.DataFrame.from_dict(self.smiles_list, orient='index')

        df = df.reset_index(drop=True)

        result = get_npclassifier_classes_from_df(df[['SMILES']], 'SMILES', 'temp_outputs')

        df = df[result.columns.tolist()]

        pd.testing.assert_frame_equal(df, result)

    def test_reading_manual_data(self):
        result = read_manual_npclassifier_input(os.path.join('test_inputs', 'manual_test_results.tsv'))
        correct = pd.read_csv(os.path.join('test_inputs', 'test_manual_input_correct.csv'))

        result = result[correct.columns.tolist()]

        pd.testing.assert_frame_equal(result,correct)


if __name__ == '__main__':
    unittest.main()
