import os
# a program to make a file for each individual in the dataset to ease the process of making a mask for each individual
directory_for_inds='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/'

# Function to create a directory if it doesn't exist
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

# Read the input file and process the data
input_file = ''
with open(input_file, 'r') as file:
    header = next(file)  # Read and skip the header line
    for line in file:
        fields = line.strip().split(',')
        pdgp_id, genus, species, froh, sex, ref_assembly = fields

        # Create a folder for each unique "ref_assembly"
        output_directory = os.path.join(directory_for_inds, ref_assembly)
        create_directory(output_directory)

        # Create a text file for the individual
        output_file = os.path.join(output_directory, f'{pdgp_id}.txt')
        open(output_file, 'w').close()  # Create an empty file


print("Files have been sorted into folders based on 'ref_assembly'.")
