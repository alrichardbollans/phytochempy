import unittest

from phytochempy.compound_properties import get_compound_ids_from_CAS_ID_from_knapsack


class TestGetCompoundIdsFromCasId(unittest.TestCase):
    def test_valid_cas_id(self):
        # Replace 'VALID_CAS_ID' with a real valid CAS ID or a mock if testing live would be unnecessary
        cas_id = '130-95-0'
        inchikey, smiles = get_compound_ids_from_CAS_ID_from_knapsack(cas_id)
        self.assertIsNotNone(inchikey)
        self.assertIsNotNone(smiles)

        self.assertEqual(inchikey,'LOUPRKONTZGTKE-GZJAXPACNA-N')
        self.assertEqual(smiles,'C=C[C@H]1C[N@]2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12')

    def test_invalid_cas_id(self):
        # Use a CAS ID that is unlikely to exist
        cas_id = 'INVALID_CAS_ID'
        inchikey, smiles = get_compound_ids_from_CAS_ID_from_knapsack(cas_id)
        self.assertIsNone(inchikey)
        self.assertIsNone(smiles)

    def test_missing_cas_id(self):
        # Test with an empty string for the CAS ID
        cas_id = ''
        inchikey, smiles = get_compound_ids_from_CAS_ID_from_knapsack(cas_id)
        self.assertIsNone(inchikey)
        self.assertIsNone(smiles)

    def test_numeric_cas_id(self):
        # Test with a numeric CAS ID provided as a string
        cas_id = '371159-83-0'
        inchikey, smiles = get_compound_ids_from_CAS_ID_from_knapsack(cas_id)
        self.assertTrue((inchikey is None or isinstance(inchikey, str)))
        self.assertTrue((smiles is None or isinstance(smiles, str)))
        self.assertEqual(inchikey,'IARIWKRWSRIQSJ-UHFFFAOYSA-N')
        self.assertEqual(smiles,'COc1cccc2[nH]c3cc(O)c(C)cc3c12')

if __name__ == '__main__':
    unittest.main()