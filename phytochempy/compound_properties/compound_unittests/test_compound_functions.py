import unittest
import urllib
from unittest.mock import patch

import numpy as np
from requests.exceptions import HTTPError

from phytochempy.metabolite_properties import resolve_cas_to_smiles, simplify_inchi_key, resolve_cas_to_inchikey


class CASresolve(unittest.TestCase):
    def setUp(self):
        self.aspirin = {'cas_id': '50-78-2', 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'}
        self.quinine = {'cas_id': '130-95-0', 'smiles': 'COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O', 'inchikey': 'LOUPRKONTZGTKE-WZBLMQSHSA-N'}

    def test_resolve_cas_to_smiles(self):
        result = resolve_cas_to_smiles(self.quinine['cas_id'])

        self.assertEqual(result, self.quinine['smiles'])
        result = resolve_cas_to_smiles(self.aspirin['cas_id'])

        self.assertEqual(result, self.aspirin['smiles'])

    def test_resolve_cas_to_inchikey(self):
        result = resolve_cas_to_inchikey(self.aspirin['cas_id'])

        self.assertEqual(result, self.aspirin['inchikey'])

        result = resolve_cas_to_inchikey(self.quinine['cas_id'])

        self.assertEqual(result, self.quinine['inchikey'])


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
