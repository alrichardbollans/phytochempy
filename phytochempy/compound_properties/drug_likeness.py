import pandas as pd


def get_bioavailability_rules(in_df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """
    :param in_df: A pandas DataFrame containing the input data.
    :param smiles_col: The name of the column in `in_df` that contains the SMILES strings.
    :return: A pandas DataFrame with the calculated bioavailability rules.

    This method calculates the bioavailability rules for a given dataset using RDKit.

    The input `in_df` should be a pandas DataFrame with at least one column containing the SMILES strings. Any rows in `in_df` with empty or NaN values in the `smiles_col` will be excluded
    * from the computation.

    The method first drops duplicate rows based on the `smiles_col`. Then, it adds a new column 'rdkit_mol' to the DataFrame using RDKit's `PandasTools.AddMoleculeColumnToFrame` function
    *.

    Next, the method calculates various properties for each molecule in the DataFrame: molecular weight (`mw`), octanol-water partition coefficient (`logP`), number of hydrogen bond accept
    *ors (`HBA`), and number of hydrogen bond donors (`HBD`). These properties are computed using RDKit's built-in functions: `MolWt`, `MolLogP`, `NumHAcceptors`, and `NumHDonors`.

    The method then applies Lipinski's rule of five to determine whether a molecule passes or fails the rule. Lipinski's rule of five states that a molecule is likely to have good oral bio
    *availability if it satisfies at least three of the following four criteria: molecular weight <= 500, hydrogen bond acceptors <= 10, hydrogen bond donors <= 5, and logP <= 5. The result
    * of this calculation is stored in a new column named 'lipinski_pass'.

    Next, the method calculates the number of rotatable bonds (`rotatable_bonds`) and the polar surface area (`polar_surface_area`) for each molecule using RDKit's `NumRotatableBonds` and
    * `TPSA` functions.

    Finally, the method applies Veber's rule to determine whether a molecule passes or fails the rule. Veber's rule states that a molecule is likely to have good oral bioavailability if
    * it satisfies the following two criteria: number of rotatable bonds <= 10 and polar surface area <= 140. The result of this calculation is stored in a new column named 'veber_pass'.

    The method selects only the columns 'veber_pass', 'lipinski_pass', and `smiles_col` from the DataFrame and returns the resulting DataFrame.
    """
    from rdkit.Chem import PandasTools
    from rdkit.Chem.Crippen import MolLogP
    from rdkit.Chem.Descriptors import MolWt
    from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors, NumRotatableBonds
    from rdkit.Chem.MolSurf import TPSA
    df = in_df[~in_df[smiles_col].isna()].drop_duplicates(subset=smiles_col, keep='first')
    PandasTools.AddMoleculeColumnToFrame(df, smiles_col, 'rdkit_mol')
    df['mw'] = df['rdkit_mol'].apply(lambda x: MolWt(x) if x is not None else None)
    df['logP'] = df['rdkit_mol'].apply(lambda x: MolLogP(x) if x is not None else None)
    df['HBD'] = df['rdkit_mol'].apply(lambda x: NumHDonors(x) if x is not None else None)
    df['HBA'] = df['rdkit_mol'].apply(lambda x: NumHAcceptors(x) if x is not None else None)

    def lipinski_check(row):
        if row['mw'] is not None:
            conditions = [row['mw'] <= 500, row['HBA'] <= 10, row['HBD'] <= 5, row['logP'] <= 5]
            if conditions.count(True) >= 3:
                return 1
            else:
                return 0
        else:
            return None

    df['lipinski_pass'] = df.apply(lambda x: lipinski_check(x), axis=1)

    df['rotatable_bonds'] = df['rdkit_mol'].apply(lambda x: NumRotatableBonds(x) if x is not None else None)
    df['polar_surface_area'] = df['rdkit_mol'].apply(lambda x: TPSA(x) if x is not None else None)

    def veber_check(row):
        if row['mw'] is not None:
            if row['rotatable_bonds'] <= 10 and row['polar_surface_area'] <= 140:
                return 1
            else:
                return 0
        else:
            return None

    df['veber_pass'] = df.apply(lambda x: veber_check(x), axis=1)
    df = df[['veber_pass', 'lipinski_pass', smiles_col]]

    return df


def add_bioavailability_rules_to_df(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """
    Add bioavailability rules to a DataFrame.

    :param df: The DataFrame to which the bioavailability rules should be added.
    :param smiles_col: The name of the column in the DataFrame that contains the SMILES strings.
    :return: The DataFrame with the added bioavailability rules.
    """
    bio_av = get_bioavailability_rules(df, smiles_col)
    bio_av = bio_av.dropna(subset=[smiles_col])

    all_metabolites_with_info = pd.merge(df, bio_av, how='left', on='SMILES')

    return all_metabolites_with_info
