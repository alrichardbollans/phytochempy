import pandas as pd
from rdkit.Chem import rdFingerprintGenerator

_mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _get_pairwise_distances_from_data(df: pd.DataFrame):
    """
    Calculate pairwise distances between molecules in a DataFrame, based on SMILES strings in 'Standard_SMILES' column.

    Assuming distances are symmetric, this calculates the lower triangle elements of the symmetric distance matrix.
    To calculate for all distances, multiply by 2.

    Using rdkit (Tool: RDKit: Open-source cheminformatics. https://www.rdkit.org), distances are calculated using the Tanimoto metric
    (Tanimoto TT (17 Nov 1958). "An Elementary Mathematical theory of Classification and Prediction". Internal IBM Technical Report. 1957)

    Parameters:
    df (pd.DataFrame): The DataFrame containing molecule information.

    Returns:
    np.ndarray: The pairwise distance matrix.

    Example:
    df = pd.DataFrame({'Standard_SMILES': ['C(C)O', 'C(C)N', 'C(C)C']})
    distmat = _get_pairwise_distances_from_data(df)
    """

    from rdkit.Chem import PandasTools
    from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
    from rdkit.DataManip.Metric import GetTanimotoDistMat

    PandasTools.AddMoleculeColumnToFrame(df, 'Standard_SMILES', 'Molecule', includeFingerprints=True)
    df = df.dropna(subset=['Molecule'])[['Molecule', 'Standard_SMILES']]


    # Produce a hashed Morgan fingerprint for each molecule
    df['morgan_fingerprint'] = df['Molecule'].apply(lambda x: _mfpgen.GetFingerprint(x))
    distmat = GetTanimotoDistMat(df['morgan_fingerprint'].values)

    return distmat


def remove_groups_with_single_compounds(working_data: pd.DataFrame, compound_grouping: str):
    counts = working_data.value_counts(compound_grouping)
    groups_with_single_compounds = pd.DataFrame({compound_grouping: counts.index, 'N': counts.values})
    groups_with_single_compounds = groups_with_single_compounds[groups_with_single_compounds['N'] < 2][compound_grouping].values.tolist()

    if len(groups_with_single_compounds) > 0:
        print('Following groups have been removed from calculations as measures for groups containing single compounds are poorly defined.')
        print(groups_with_single_compounds)

    out_df = working_data[~working_data[compound_grouping].isin(groups_with_single_compounds)]
    return out_df


def calculate_FAD_measures(df: pd.DataFrame, compound_grouping: str):
    """

    This method calculates FAD-related measures for each unique taxon in a given DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing a 'compound_grouping' column by which compounds are grouped and a Standard_SMILES column
    for compounds.
    - compound_grouping (str): The column name in the DataFrame that represents the how compounds in the dataframe should be grouped.

    Returns:
    - pd.DataFrame: A DataFrame containing the FAD measures for each group.

    """

    duplicate_issues = df[df.duplicated(subset=['Standard_SMILES', compound_grouping])]
    if len(duplicate_issues) > 0:
        print(
            f'WARNING: Standard_SMILES records are duplicated within the grouping of {compound_grouping}, duplicates will be removed for diversity calculations')

        df = df.drop_duplicates(
            subset=[compound_grouping, 'Standard_SMILES'],
            keep='first')

    df = remove_groups_with_single_compounds(df, compound_grouping)

    groups = df.groupby(compound_grouping)
    results = []

    for taxon, taxon_data in groups:
        N = len(taxon_data)
        if N < 2:
            raise ValueError(
                f'Measures for taxa with single compounds are poorly defined. '
                f'In this case for taxon "{taxon}", '
                f'either remove these from your data, or make a pull request :)'
            )
        distances = _get_pairwise_distances_from_data(taxon_data)
        FAD = distances.sum() * 2
        num_distances = len(distances) * 2
        assert num_distances == N ** 2 - N
        MFAD = FAD / N
        APWD = FAD / num_distances

        results.append({
            compound_grouping: taxon,
            'FAD': FAD,
            'MFAD': MFAD,
            'APWD': APWD,
            'GroupSize_FAD': N
        })

    out_df = pd.DataFrame(results)

    return out_df
