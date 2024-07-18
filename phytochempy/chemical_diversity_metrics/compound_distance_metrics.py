import pandas as pd



def _get_pairwise_distances_from_data(df: pd.DataFrame):
    """
    Calculate pairwise distances between molecules in a DataFrame, based on SMILES strings in 'Standard_SMILES' column.

    Using rdkit (Tool: RDKit: Open-source cheminformatics. https://www.rdkit.org), distances are calculated using the Tanimoto metric
    (Tanimoto TT (17 Nov 1958). "An Elementary Mathematical theory of Classification and Prediction". Internal IBM Technical Report. 1957)

    Parameters:
    df (pd.DataFrame): The DataFrame containing molecule information.

    Returns:
    np.ndarray: The pairwise distance matrix.

    Example:
    df = pd.DataFrame({'Molecule': ['CCO', 'CCN', 'CCC'],
                       'Standard_SMILES': ['C(C)O', 'C(C)N', 'C(C)C']})
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


def calculate_FAD_measures(df: pd.DataFrame, taxon_grouping: str = 'Genus'):
    """

    This method calculates FAD-related measures for each unique taxon in a given DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing taxon data.
    - taxon_grouping (str): The column name in the DataFrame that represents the taxonomic grouping. Default is 'Genus'.

    Returns:
    - pd.DataFrame: A DataFrame containing the FAD measures for each taxon.

    Raises:
    - ValueError: If there is a taxon with only one compound.

    """
    FAD_outputs = {}
    MFAD_outputs = {}
    APWD_outputs = {}
    N_outputs = {}
    for taxon in df[taxon_grouping].unique():
        taxon_data = df[df[taxon_grouping] == taxon]
        if len(taxon_data) > 1:
            distances = _get_pairwise_distances_from_data(taxon_data)
            FAD_outputs[taxon] = distances.sum()
            MFAD_outputs[taxon] = distances.sum() / len(taxon_data)
            APWD_outputs[taxon] = distances.sum() / len(distances)
            N_outputs[taxon] = len(taxon_data)
        else:
            # FAD_outputs[taxon] = 0
            # MFAD_outputs[taxon] = 0
            # APWD_outputs[taxon] = 0
            # N_outputs[taxon] = len(taxon_data)
            raise ValueError(f'Measures for taxa with single compounds are poorly defined. In this case for {taxon}'
                             f'Either remove these from your data, or make a pull request :)')

    out_df = pd.DataFrame.from_dict(FAD_outputs, orient='index', columns=['FAD'])

    out_df['MFAD'] = MFAD_outputs.values()
    out_df['APWD'] = APWD_outputs.values()
    out_df['N'] = N_outputs.values()

    from sklearn.preprocessing import MinMaxScaler
    for index in ['MFAD', 'APWD', 'FAD']:
        scaler = MinMaxScaler()
        out_df[index + '_minmax'] = scaler.fit_transform(out_df[[index]])
        print(index)

    out_df = out_df.reset_index(names=[taxon_grouping])
    return out_df
