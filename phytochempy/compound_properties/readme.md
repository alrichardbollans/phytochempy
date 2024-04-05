```python

import os
import pandas as pd
from phytochempy.data_compilation_utilities import get_wikidata,get_knapsack_data,merge_and_tidy_compound_datasets,get_manual_files_to_upload,\
add_npclassifier_info,add_chembl_data,add_bioavailability_info,add_manual_info_files,tidy_final_dataset
### Example workflow

# Define context
comp_id_column = 'InChIKey'  # Where appropriate, which column should be used to determine compound uniqueness. This is not applicable to some properties, e.g. where SMILES must be used to generate data
families = ['Pandanaceae']
wiki_data_id_for_order = 'Q736182'
temp_outputs_folder = 'temp'
tidied_outputs_folder = 'tidied'

## Get compound-taxa pair data
get_wikidata(wiki_data_id_for_order, os.path.join(temp_outputs_folder, 'wikidata.csv'), os.path.join(tidied_outputs_folder, 'wikidata.csv'))
get_knapsack_data(families, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'knapsack_data.csv'))

## Merge and tidy the data
tidy_wiki_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'wikidata.csv'), index_col=0)
tidy_knapsack_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'knapsack_data.csv'), index_col=0)
all_compounds_in_taxa = merge_and_tidy_compound_datasets([tidy_wiki_data, tidy_knapsack_data],
                                                         os.path.join(tidied_outputs_folder, 'merged_data.csv'))

get_manual_files_to_upload(all_compounds_in_taxa, temp_outputs_folder)

## Add extra information related to the compound properties

# These steps can be included/removed as needed
# For the longer processes, to avoid repeats you can simply read the associated temp_output if the step has already been run
with_npclass_classes = add_npclassifier_info(all_compounds_in_taxa, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'npclassifier.csv'))

with_chembl_data = add_chembl_data(with_npclass_classes, os.path.join(temp_outputs_folder, 'chembl.csv'), compound_id_column=comp_id_column)
# with_chembl_data = pd.read_csv(os.path.join(temp_outputs_folder, 'chembl.csv'), index_col=0)
with_bioavailibility = add_bioavailability_info(with_chembl_data, os.path.join(tidied_outputs_folder, 'bioavailibility.csv'))
# with_bioavailibility = pd.read_csv(os.path.join(tidied_outputs_folder, 'bioavailibility.csv'), index_col=0)
# Issues with classyfire servers
# with_classyfire_classes = add_classyfire_info(with_chembl_data, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'classyfire.csv'))

## This step requires some manual input

all_info = add_manual_info_files(with_bioavailibility,
                                 maip_output_file=os.path.join(tidied_outputs_folder, 'example_maip_file.csv'))

### Then tidy and output final dataset
tidy_final_dataset(all_info, tidied_outputs_folder, os.path.join('outputs', 'all_taxa_compound_data.csv'), comp_id_column)
```