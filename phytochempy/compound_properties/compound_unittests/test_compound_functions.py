import os
import unittest

import numpy as np
import pandas as pd

from phytochempy.compound_properties import resolve_cas_to_smiles, simplify_inchi_key, resolve_cas_to_inchikey, get_smiles_and_inchi_from_cas_ids, \
    add_CAS_ID_translations_to_df, fill_match_ids, standardise_smiles_to_MAIP_smiles, standardise_SMILES


class CASresolve(unittest.TestCase):
    def setUp(self):
        self.aspirin = {'CAS ID': '50-78-2', 'SMILES': 'CC(=O)Oc1ccccc1C(O)=O', 'InChIKey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'}
        self.arte = {'CAS ID': '63968-64-9', 'SMILES': 'C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4',
                     'InChIKey': 'BLUAFEHZUWYNDE-NNWCWBAJSA-N'}
        self.artem = {'CAS ID': '71963-77-4', 'SMILES': 'CO[C@H]1O[C@@H]2O[C@@]3(C)CC[C@H]4[C@H](C)CC[C@@H]([C@H]1C)[C@@]24OO3',
                      'InChIKey': 'SXYIRMFQILZOAM-HVNFFKDJSA-N'}

    def test_resolve_cas_to_smiles(self):
        result = resolve_cas_to_smiles(self.artem['CAS ID'])

        self.assertEqual(result, self.artem['SMILES'])

        result = resolve_cas_to_smiles(self.arte['CAS ID'])

        self.assertEqual(result, self.arte['SMILES'])
        result = resolve_cas_to_smiles(self.aspirin['CAS ID'])

        self.assertEqual(result, self.aspirin['SMILES'])

    def test_resolve_cas_to_inchikey(self):
        result = resolve_cas_to_inchikey(self.artem['CAS ID'])

        self.assertEqual(result, self.artem['InChIKey'])

        result = resolve_cas_to_inchikey(self.aspirin['CAS ID'])

        self.assertEqual(result, self.aspirin['InChIKey'])

        result = resolve_cas_to_inchikey(self.arte['CAS ID'])

        self.assertEqual(result, self.arte['InChIKey'])

    def test_list_resolve(self):
        result = get_smiles_and_inchi_from_cas_ids([self.aspirin['CAS ID'], self.arte['CAS ID'], self.artem['CAS ID']], 'temp_outputs')
        correct_df = pd.DataFrame([self.aspirin, self.arte, self.artem])
        pd.testing.assert_frame_equal(result, correct_df)

    def test_df_resolve(self):
        correct_df = pd.DataFrame([self.aspirin, self.arte, self.artem, self.aspirin])
        result = add_CAS_ID_translations_to_df(correct_df[['CAS ID']], 'CAS ID', 'temp_outputs')
        pd.testing.assert_frame_equal(result, correct_df)


class Testkeys(unittest.TestCase):

    def test_simplify_inchi_key_with_valid_input(self):
        input_data = "INCHI_KEY_SAMPLE_STRING"
        expected_output = "INCHI_KEY_SAMP"
        result = simplify_inchi_key(input_data)
        self.assertEqual(result, expected_output)

    def test_simplify_inchi_key_with_less_than_14_characters(self):
        input_data = "INCHI_KEY"
        expected_output = "INCHI_KEY"
        result = simplify_inchi_key(input_data)
        self.assertEqual(result, expected_output)

    def test_simplify_inchi_key_with_no_input(self):
        input_data = ""
        expected_output = ""
        result = simplify_inchi_key(input_data)
        self.assertEqual(result, expected_output)

    def test_simplify_inchi_key_with_None_input(self):
        input_data = None
        result = simplify_inchi_key(input_data)
        if result is not None:
            raise ValueError


class TestFillIds(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            'SMILES': ['CCO', 'BBX', 'C#C', 'CC', np.nan, 'CCCC','CC'],
            'InChIKey': ['OnyigWHJXaIubT', 'OnyigWHJXaIubT', 'QFJUORZEWOZLNI', np.nan, 'BTJIUGUIPKHNIE', 'BBX','THIS'],
            'CAS ID': [np.nan, '89-26', '75-25-2', np.nan, '75-00-3', '75-00-3','DEAD'],
        })

    def test_fill_CAS_IDs(self):
        actual_df = fill_match_ids(df=self.df, given_col='CAS ID')
        expected_df = self.df.copy()
        expected_df['CAS ID'] = ['89-26', '89-26', '75-25-2', 'DEAD', '75-00-3', '75-00-3','DEAD']
        pd.testing.assert_frame_equal(actual_df, expected_df)

    def test_fill_SMILES(self):
        actual_df = fill_match_ids(df=self.df, given_col='SMILES')
        expected_df = self.df.copy()
        expected_df['SMILES'] = ['CCO', 'BBX', 'C#C', 'CC', 'CCCC', 'CCCC','CC']
        pd.testing.assert_frame_equal(actual_df, expected_df)

    def test_fill_InChIKey(self):
        actual_df = fill_match_ids(df=self.df, given_col='InChIKey')
        expected_df = self.df.copy()
        expected_df['InChIKey'] = ['OnyigWHJXaIubT', 'OnyigWHJXaIubT', 'QFJUORZEWOZLNI','THIS', 'BTJIUGUIPKHNIE', 'BBX','THIS']
        pd.testing.assert_frame_equal(actual_df, expected_df)

    def test_with_empty_dataframe(self):
        df_empty = pd.DataFrame(columns=['SMILES', 'InChIKey', 'CAS ID'])
        actual_df = fill_match_ids(df=df_empty, given_col='CAS ID')
        pd.testing.assert_frame_equal(actual_df, df_empty)
class TestStandardiseSmiles(unittest.TestCase):
    def test_standardise_smiles_valid(self):
        self.assertEqual(standardise_smiles_to_MAIP_smiles('CO'), 'CO')
        self.assertEqual(standardise_smiles_to_MAIP_smiles("[Na]OC(=O)Cc1ccc(C[NH3+])cc1.c1nnn[n-]1.O"), 'NCc1ccc(CC(=O)O)cc1')

    def test_standardise_smiles_invalid(self):
        self.assertEqual(standardise_smiles_to_MAIP_smiles('Not a SMILES'), None)
        self.assertEqual(standardise_smiles_to_MAIP_smiles('NotaSMILES'), None)
        self.assertEqual(standardise_smiles_to_MAIP_smiles(None), None)

class TestrdkitStandardiseSmiles(unittest.TestCase):
    def test_standardise_smiles_valid(self):
        self.assertEqual(standardise_SMILES('CO'), 'CO')
        self.assertEqual(standardise_SMILES("[Na]OC(=O)Cc1ccc(C[NH3+])cc1.c1nnn[n-]1.O"), '[NH3+]Cc1ccc(CC(=O)[O-])cc1')

    def test_standardise_smiles_invalid(self):
        self.assertEqual(standardise_SMILES('Not a SMILES'), None)
        self.assertEqual(standardise_SMILES('NotaSMILES'), None)
        self.assertEqual(standardise_SMILES(None), None)
if __name__ == '__main__':
    unittest.main()
