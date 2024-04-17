import os.path
import unittest

import pandas as pd
from pkg_resources import resource_filename
from tqdm import tqdm
from wcvpy.wcvp_download import wcvp_accepted_columns, get_all_taxa
from wcvpy.wcvp_name_matching import get_accepted_info_from_names_in_column

from phytochempy.compound_properties import COMPOUND_NAME_COLUMN
from phytochempy.knapsack_searches import get_knapsack_formulas_for_compound, get_knapsack_compounds_for_taxon

input_test_dir = resource_filename(__name__, 'test_inputs')
test_output_dir = resource_filename(__name__, 'test_outputs')


class SpecificTaxa(unittest.TestCase):

    def taxon_test(self, name, known_vals):
        table = get_knapsack_compounds_for_taxon(name)

        metabolites = table[COMPOUND_NAME_COLUMN].values.tolist()

        self.assertEqual(metabolites, known_vals)

        upper = get_knapsack_compounds_for_taxon(name.upper())
        lower = get_knapsack_compounds_for_taxon(name.lower())
        pd.testing.assert_frame_equal(upper.drop(columns=['knapsack_search_term']), lower.drop(columns=['knapsack_search_term']))

    def test_no_info(self):
        apoc_no_info_searches = ['XXX', 'Mandevilla', 'Tromotriche', 'Orthanthera', 'Bustelma', 'Cystostemma',
                                 'Widgrenia', 'Dactylostelma', 'Gothofreda', 'Oxystelma', 'Rhombonema',
                                 'Pattalias', 'Macbridea', 'Vadulia', 'Pentacyphus', 'Pentatropis',
                                 'Acustelma']

        all_genera_df = pd.DataFrame()
        for i in tqdm(range(len(apoc_no_info_searches)), desc="Searching genera in Knapsack…", ascii=False,
                      ncols=80):
            genus = apoc_no_info_searches[i]
            genus_table = get_knapsack_compounds_for_taxon(genus)
            if len(genus_table.index) > 0:
                all_genera_df = pd.concat([all_genera_df, genus_table])
        if len(all_genera_df.index) > 0:
            acc_df = get_accepted_info_from_names_in_column(all_genera_df, 'Organism',
                                                            families_of_interest=['Apocynaceae'])
            acc_df.to_csv(os.path.join(test_output_dir, 'apoc_no_info.csv'))
            acc_df = acc_df.dropna(subset=[wcvp_accepted_columns['id']])

            self.assertEqual(len(acc_df.index), 0)

        rub_no_info = ['Thouarsiora', 'Bemsetia', 'Tsiangia', 'Schetti', 'Taligalea', 'Zuccarinia', 'Jackia',
                       'Janotia', 'Jovetia', 'Joosia', 'Gouldia', 'Neobaumannia', 'Dentillaria', 'Kohautia']

        all_genera_df = pd.DataFrame()
        for i in tqdm(range(len(rub_no_info)), desc="Searching genera in Knapsack…", ascii=False,
                      ncols=80):
            genus = rub_no_info[i]
            genus_table = get_knapsack_compounds_for_taxon(genus)
            if len(genus_table.index) > 0:
                all_genera_df = pd.concat([all_genera_df, genus_table])
        if len(all_genera_df.index) > 0:
            acc_df = get_accepted_info_from_names_in_column(all_genera_df, 'Organism',
                                                            families_of_interest=['Rubiaceae'])
            acc_df.to_csv(os.path.join(test_output_dir, 'rub_no_info.csv'))
            acc_df = acc_df.dropna(subset=[wcvp_accepted_columns['id']])

            self.assertEqual(len(acc_df.index), 0)

    def test_meta_known_taxa(self):
        known_vals = ['(+)-Aspidospermidine', '(+)-Fendlerine', 'Alalakine', 'Aspidolimidine', 'Limaspermine',
                      'N-Acetylaspidospermidine']
        self.taxon_test('Aspidosperma album', known_vals)

        # Test some synonyms with vals
        self.taxon_test('Curarea', ['Isochondodendrine', 'Isochondodendrine', 'Isochondodendrine', 'Limacine',
                                    "(-)-Limacine 2'-beta-N-oxide", '(-)-Curine', '(-)-Krukovine',
                                    "(-)-Limacine 2'-alpha-N-oxide", '(-)-Limacine, 2-beta-N-oxide',
                                    '(+)-Candicusine', 'Limacusine'])

        self.taxon_test('Strychnos lucida', ['Brucine', 'Strychnine', 'Acanthoside B', '3-O-Caffeoylquinic acid', 'Loganin',
                                             'Adenosine', 'Loganate', 'Cantleyoside', 'Sylvestroside I', 'Sweroside',
                                             'Tachioside', 'alpha-Colubrine', 'beta-Colubrine', 'Brucine N-oxide',
                                             'Diaboline', 'Normacusine B', 'Pseudobrucine', 'Pseudostrychnine',
                                             '11-Methoxydiaboline', '3,4-di-O-caffeoylquinic acid',
                                             'Ligustrinoside', 'Picconioside I', 'Secoxyloganin', 'Staunoside C',
                                             'Triplostoside A'])


class Families(unittest.TestCase):
    def setUp(self):
        self.family = 'Pandanaceae'  # Use a small test family
        self.temp_csv = os.path.join(test_output_dir, 'knapsack_poly.csv')
        self.final_csv = os.path.join(test_output_dir, 'knapsack_tidied_poly.csv')
        # get_knapsack_compounds_in_family(self.family, self.temp_csv)
        # tidy_knapsack_results(self.temp_csv, self.final_csv, self.family, add_smiles_and_inchi=False)

    def test_family_example(self):
        metas_df = pd.read_csv(self.final_csv, index_col=0)
        correct_df = pd.read_csv(os.path.join(input_test_dir, 'knapsack_tidied_poly.csv'), index_col=0)
        pd.testing.assert_frame_equal(metas_df, correct_df)

    def test_all_genera_searched_for(self):
        metas_df = pd.read_csv(self.final_csv, index_col=0)

        rub_genera = get_all_taxa(families_of_interest=[self.family], ranks=['Genus'], accepted=True)[
            wcvp_accepted_columns['name']].values
        genera_not_in_knapsack = ['Freycinetia', 'Sararanga', 'Martellidendron', 'Benstonea']
        for genus in rub_genera:
            if genus not in genera_not_in_knapsack:
                self.assertIn(genus, metas_df['accepted_parent'].values)


class Formulas(unittest.TestCase):
    def test_get_formulas_for_metabolite(self):
        formulas = get_knapsack_formulas_for_compound('quinine')
        correct = ['C20H24N2O2', 'C19H22N2O2', 'C19H24N2O2', 'C23H26NO7']

        self.assertEqual(sorted(formulas), sorted(correct))


if __name__ == '__main__':
    unittest.main()
