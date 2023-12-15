from cgi import test
from json import load
from operator import ge
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
                dictionary_of_species_n_pdid[ref_assembly] = {genus: [species]}
            if genus not in dictionary_of_species_n_pdid[ref_assembly]:
                dictionary_of_species_n_pdid[ref_assembly][genus] = [species]
            if species not in dictionary_of_species_n_pdid[ref_assembly][genus]:
                dictionary_of_species_n_pdid[ref_assembly][genus].append(species)
            # else:
            #     dictionary_of_species_n_pdid[ref_assembly][genus].append(species)
    return dictionary_of_species_n_pdid

def load_in_datasets(files) -> pd.DataFrame: # This function loads in the datasets
    dfs = []
    for genus in files:
        for species, file in  files[genus].items():
            df = pd.read_csv(file, sep='\t')
            df['species'] = species
            df['genus'] = genus
            dfs += [df]
    return pd.concat(dfs)

def gather_pi_files(ref_assembly,window_size, dictionary_of_species): # This function generates a dictionary of files
    dictionary_of_files = {}
    for genus in dictionary_of_species[ref_assembly]:
        for species in dictionary_of_species[ref_assembly][genus]:
            print(genus, species)
            if genus not in dictionary_of_files:
                dictionary_of_files[genus] = {species: f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi'}
            elif species not in dictionary_of_files[genus]:
                dictionary_of_files[genus][species] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi'
    return dictionary_of_files

def gather_callable_frac_files(ref_assembly,window_size, dictionary_of_species): # This function generates a dictionary of files
    dictionary_of_files = {}
    for genus in dictionary_of_species[ref_assembly]:
        for species in dictionary_of_species[ref_assembly][genus]:
            if genus not in dictionary_of_files:
                dictionary_of_files[genus] = {}
            if species not in dictionary_of_files[genus]:
                dictionary_of_files[genus][species] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/callable_fraction/{ref_assembly}/{species}/{window_size}_callable_fraction.csv'
    return dictionary_of_files

def gather_recombination_files(ref_assembly,window_size, dictionary_of_species): # This function generates a dictionary of files
    dictionary_of_files = {}
    for genus in dictionary_of_species[ref_assembly]:
        for species in dictionary_of_species[ref_assembly][genus]:
            if genus not in dictionary_of_files:
                dictionary_of_files[genus] = {}
            if species not in dictionary_of_files[genus]:
                dictionary_of_files[genus][species] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/recombinationrates/{ref_assembly}/{species}/nonpar/{species}_{window_size}_recomb.bed'
    return dictionary_of_files


def write_to_file(df, output_file): # This function writes the dataframe to a file
    df.to_csv(output_file, sep='\t', index=False)




def wrapper(ref_assembly, window_size, dictionary_of_species): 
    recombinationrates_files = load_in_datasets(gather_recombination_files(ref_assembly, window_size, dictionary_of_species))
    callable_frac_files = load_in_datasets(gather_callable_frac_files(ref_assembly, window_size, dictionary_of_species))
    frac_n_recombination = pd.merge(callable_frac_files, recombinationrates_files, on=['chr', 'start', 'end', 'genus','species'])
    pi_files = load_in_datasets(gather_pi_files(ref_assembly, window_size, dictionary_of_species))
    pi_files = pi_files.rename(columns={'CHROM': 'chr', 'BIN_START': 'start', 'BIN_END': 'end'})
    merged = pd.merge(frac_n_recombination, pi_files, on=['chr', 'start', 'end', 'genus','species'])
    return merged

## LOAD IN DATA
ref_assembly = sys.argv[1]
window_size = sys.argv[2]
output_file = sys.argv[3]

# ref_assembly = 'Pithecia_pithecia'
# window_size = 100000
# output_file = 'test.csv'


## RUN THE SCRIPT
dictionary_of_species = create_input_paths(sample_overview)
# print(dictionary_of_species['Pithecia_pithecia'])
# print(load_in_datasets(gather_pi_files(ref_assembly, window_size, dictionary_of_species)))
merged = wrapper(ref_assembly, window_size, dictionary_of_species)
write_to_file(merged, output_file)

