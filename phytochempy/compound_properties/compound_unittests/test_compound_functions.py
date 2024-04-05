import unittest

import pandas as pd

from phytochempy.compound_properties import resolve_cas_to_smiles, simplify_inchi_key, resolve_cas_to_inchikey, get_smiles_and_inchi_from_cas_ids, \
    add_CAS_ID_translations_to_df


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
        correct_df = pd.DataFrame([self.aspirin, self.arte, self.artem,self.aspirin])
        result = add_CAS_ID_translations_to_df(correct_df[['CAS ID']], 'CAS ID','temp_outputs')
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


if __name__ == '__main__':
    unittest.main()
