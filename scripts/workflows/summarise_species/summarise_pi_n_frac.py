from cgi import test
from json import load
import sys
from weakref import ref
from numpy import append
import pandas as pd
from typing import Dict



sample_overview = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'

def create_input_paths(input_file):
    dictionary_of_species_n_pdid = {}
    with open(input_file, 'r') as file:
        header = next(file)  # Read and skip the header line
        for line in file:
            fields = line.strip().split(',')
            pdgp_id, genus, species, froh, sex, ref_assembly = fields
            if not ref_assembly:  # Check if ref_assembly is empty
                ref_assembly = 'unknown'
            # Add the pdgp_id to the appropriate category in the dictionary
            if ref_assembly not in dictionary_of_species_n_pdid:
                dictionary_of_species_n_pdid[ref_assembly] = [species]
            else:
                dictionary_of_species_n_pdid[ref_assembly].append(species)
    return dictionary_of_species_n_pdid

def load_in_datasets(files: Dict[str, str]) -> pd.DataFrame: # This function loads in the datasets
    dfs = []
    for key, file in files.items(): # key is pdgp_id and file is the path to the file
        df = pd.read_csv(file, sep='\t')
        df['species'] = key
        dfs += [df]
    return pd.concat(dfs)

def gather_pi_files(ref_assembly,window_size, dictionary_of_species): # This function generates a dictionary of files
    dictionary_of_files = {}
    for species in dictionary_of_species[ref_assembly]:
        if species not in dictionary_of_files:
            dictionary_of_files[species] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi'
        else:
            'this is not working'
    return dictionary_of_files

def gather_callable_frac_files(ref_assembly,window_size, dictionary_of_species): # This function generates a dictionary of files
    dictionary_of_files = {}
    for species in dictionary_of_species[ref_assembly]:
        if species not in dictionary_of_files:
            dictionary_of_files[species] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/callable_fraction/{ref_assembly}/{species}/{window_size}_callable_fraction.csv'
        else:
            'this is not working'
    return dictionary_of_files



def write_to_file(df, output_file): # This function writes the dataframe to a file
    df.to_csv(output_file, sep='\t', index=False)




def wrapper(ref_assembly, window_size, dictionary_of_species): 
    pi_files = load_in_datasets(gather_pi_files(ref_assembly, window_size, dictionary_of_species))
    callable_frac_files = load_in_datasets(gather_callable_frac_files(ref_assembly, window_size, dictionary_of_species))
    pi_files = pi_files.rename(columns={'CHROM': 'chr', 'BIN_START': 'start', 'BIN_END': 'end'})
    merged = pd.merge(callable_frac_files, pi_files, on=['chr', 'start', 'end', 'species'])
    return merged

## LOAD IN DATA
ref_assembly = sys.argv[1]
window_size = sys.argv[2]
output_file = sys.argv[3]

## RUN THE SCRIPT
dictionary_of_species = create_input_paths(sample_overview)
merged = wrapper(ref_assembly, window_size, dictionary_of_species)
write_to_file(merged, output_file)

