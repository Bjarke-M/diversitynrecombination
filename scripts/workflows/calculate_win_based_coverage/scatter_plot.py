import pandas as pd
import os

# Define the input path
input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'
data_dir = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/'

# Load metadata into dictionaries
dictionary_of_inds = {}
dictionary_of_species = {}

with open(input_file, 'r') as file:
    header = next(file)
    for line in file:
        pdgp_id, genus, species, froh, sex, ref_assembly = line.strip().split(',')
        species = f"{genus}_{species}"
        dictionary_of_inds.setdefault(ref_assembly, []).append(pdgp_id)
        dictionary_of_species.setdefault(species, []).append(pdgp_id)

# Define the list of window sizes
window_list = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 2000000, 10000000]

# Create a DataFrame to store the results
result_list = []

# Iterate through individuals and window sizes
for ref_assembly, pd_ids in dictionary_of_inds.items():
    for pd_id in pd_ids:
        for window_size in window_list:
            # Generate the file path
            file_path = os.path.join(data_dir, ref_assembly, pd_id, f"{pd_id}_coverage_{window_size}.txt")
            # Check if the file exists
            if os.path.isfile(file_path):
                # Read the data from the file
                data = pd.read_csv(file_path, header=None, names=['chrA', 'startA', 'stopA', 'chrB', 'startB', 'stopB', 'n'], sep='\t')
                
                # Calculate n/window_size
                n_over_window = data['n'].sum() / window_size
                
                # Append the result directly to the list
                result_list.append([ref_assembly, pd_id, species, n_over_window, data['chrA'].iloc[0]])

# Create a DataFrame from the result list
result_df = pd.DataFrame(result_list, columns=['ref_assembly', 'pd_id', 'species', 'n/window_size', 'chrA'])

print(result_df)
