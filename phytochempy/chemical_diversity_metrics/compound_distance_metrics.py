import pandas as pd


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
    df['morgan_fingerprint'] = df['Molecule'].apply(lambda x: GetMorganFingerprintAsBitVect(x, 2))
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

    FAD_outputs = {}
    MFAD_outputs = {}
    APWD_outputs = {}
    N_outputs = {}
    for taxon in df[compound_grouping].unique():
        taxon_data = df[df[compound_grouping] == taxon]
        if len(taxon_data) > 1:
            distances = _get_pairwise_distances_from_data(taxon_data)
            FAD = distances.sum() * 2
            FAD_outputs[taxon] = FAD

            N = len(taxon_data)

            number_of_distances = len(distances)*2
            assert number_of_distances == N ** 2 - N

            MFAD_outputs[taxon] = FAD / N

            APWD_outputs[taxon] = FAD / number_of_distances
            N_outputs[taxon] = N
        else:
            # FAD_outputs[taxon] = 0
            # MFAD_outputs[taxon] = 0
            # APWD_outputs[taxon] = 0
            # N_outputs[taxon] = len(taxon_data)
            raise ValueError(f'Measures for taxa with single compounds are poorly defined. In this case for taxon "{taxon}",'
                             f'Either remove these from your data, or make a pull request :)')

    out_df = pd.DataFrame.from_dict(FAD_outputs, orient='index', columns=['FAD'])

    out_df['MFAD'] = MFAD_outputs.values()
    out_df['APWD'] = APWD_outputs.values()
    out_df['GroupSize_FAD'] = N_outputs.values()  # The size of the compound groups used in these calculations.

    out_df = out_df.reset_index(names=[compound_grouping])
    return out_df
