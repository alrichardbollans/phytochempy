# PlantChemicalDiversity

A repository for analysing chemical diversity in vascular plants.

## Installation

With pip, run:

`pip install git+https://github.com/alrichardbollans/phytochempy.git@1.0`

When using this package, please cite the appropriate data sources. These are detailed in References below.

## Usage Examples

## References

### Compound-Organism Pairs

#### KNApSAcK

Farit Mochamad Afendi, Taketo Okada,, Mami Yamazaki, Aki-Hirai-Morita, Yukiko Nakamura,
Kensuke Nakamura, Shun Ikeda, Hiroki Takahashi, Md. Altaf-Ul-Amin, Latifah, Darusman, Kazuki
Saito, Shigehiko Kanaya, “KNApSAcK Family Databases: Integrated Metabolite-Plant Species
Databases for Multifaceted Plant Research,” Plant Cell Physiol., 53, e1(1-12), (2012). doi:
10.1093/pcp/pcr165.

#### WikiData

Denny Vrandečić and Markus Krötzsch, ‘Wikidata: A Free Collaborative Knowledgebase’, Communications of the ACM 57, no. 10 (2014): 78–85.

and related initiative:
Adriano Rutz, Maria Sorokina, Jakub Galgonek, Daniel Mietchen, Egon Willighagen, Arnaud Gaudry, James G Graham,
Ralf Stephan, Roderic Page, Jiří Vondrášek, Christoph Steinbeck, Guido F Pauli, Jean-Luc Wolfender, Jonathan Bisson,
Pierre-Marie Allard (2022)
The LOTUS initiative for open knowledge management in natural products research.
eLife 11:e70780. https://doi.org/10.7554/eLife.70780

### Compound Properties

##### NPClassifier

Hyun Woo Kim et al., ‘NPClassifier: A Deep Neural Network-Based Structural Classification Tool for Natural Products’,
Journal of Natural Products 84, no. 11 (26 November 2021): 2795–2807, https://doi.org/10.1021/acs.jnatprod.1c00399.

GNPS Workflow:
Wang, Mingxun, Jeremy J. Carver, Vanessa V. Phelan, Laura M. Sanchez, Neha Garg, Yao Peng, Don Duy Nguyen et al. 'Sharing and community curation of
mass spectrometry data with Global Natural Products Social Molecular Networking.' Nature biotechnology 34, no. 8 (2016): 828-837.

#### ClassyFire

Yannick Djoumbou Feunang et al., ‘ClassyFire: Automated Chemical Classification with a Comprehensive, Computable Taxonomy’,
Journal of Cheminformatics 8, no. 1 (December 2016): 61, https://doi.org/10.1186/s13321-016-0174-y.
http://classyfire.wishartlab.com

##### ChEMBL Assays

The data supplied with the package was accessed 13/12/2023.

ChEMBL Database:
Mendez D, Gaulton A, Bento AP, Chambers J, De Veij M, Félix E, Magariños MP, Mosquera JF, Mutowo P, Nowotka M, Gordillo-Marañón M, Hunter F, Junco L,
Mugumbate G, Rodriguez-Lopez M, Atkinson F, Bosc N, Radoux CJ, Segura-Cabrera A, Hersey A, Leach AR. (2019) 'ChEMBL: towards direct deposition of
bioassay data.' Nucleic Acids Res., 47(D1) D930-D940.

ChEMBL Web Services:
Davies M, Nowotka M, Papadatos G, Dedman N, Gaulton A, Atkinson F, Bellis L, Overington JP. (2015) 'ChEMBL web services: streamlining access to drug
discovery data and utilities.' Nucleic Acids Res., 43(W1) W612-W620.

#### MAIP

Nicolas Bosc et al., “MAIP: A Web Service for Predicting Blood‐stage Malaria Inhibitors,”
Journal of Cheminformatics 13, no. 1 (December 2021): 13, https://doi.org/10.1186/s13321-021-00487-2.

##### Resolving CAS IDs

CIRpy: https://github.com/mcs07/CIRpy was used to resolve given CAS Registry IDs:
CAS, ‘CAS Registry System’, Journal of Chemical Information and Computer Sciences 18, no. 1 (1978): 58–58, https://doi.org/10.1021/ci60013a609
using the NCI/CADD Chemical Identifier Resolver:https://cactus.nci.nih.gov/chemical/structure

#### Bioavailability Metrics

Veber:
Daniel F. Veber et al., ‘Molecular Properties That Influence the Oral Bioavailability of Drug Candidates’, Journal of Medicinal Chemistry 45, no. 12 (
1 June 2002): 2615–23, https://doi.org/10.1021/jm020017n.

Lipinski:
Christopher A Lipinski et al., ‘Experimental and Computational Approaches to Estimate Solubility and Permeability in Drug Discovery and Development
Settings’, Advanced Drug Delivery Reviews 23, no. 1–3 (1997): 3–25.

Tool:
RDKit: Open-source cheminformatics. https://www.rdkit.org
