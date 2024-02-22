# This script takes in a species name and a window size and outputs a file with the callable fraction for each window
from asyncio import gather
import sys
from unittest import result
from numpy import median
import pandas as pd
from typing import Dict

def generatedictionary(input_file): # This function generates a dictionary of species and pdgp_ids
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
                dictionary_of_species_n_pdid[ref_assembly] = {species: [pdgp_id]}
            elif species not in dictionary_of_species_n_pdid[ref_assembly]:
                dictionary_of_species_n_pdid[ref_assembly][species] = [pdgp_id]
            else:
                dictionary_of_species_n_pdid[ref_assembly][species].append(pdgp_id)
    return dictionary_of_species_n_pdid 


def gather_files(ref_assembly, species, window_size, dictionary_of_pds): # This function generates a dictionary of files
    dictionary_of_files = {}
    for pd in dictionary_of_pds[ref_assembly][species]:
        if pd not in dictionary_of_files:
            dictionary_of_files[pd] = f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/tajimasd/mask/{ref_assembly}/{species}/{pd}/nonpar/csv/{pd}_masked_{window_size}.csv'
        else:
            'this is not working'
    return dictionary_of_files

def load_in_datasets(files: Dict[str, str]) -> pd.DataFrame: # This function loads in the datasets
    dfs = []
    for key, file in files.items(): # key is pdgp_id and file is the path to the file
        df = pd.read_csv(file)
        df['species'] = key
        dfs += [df]
    return pd.concat(dfs)

def summarise_df(df): # This function summarises the dataframe
    grouped_df = df.groupby(["chr", "start", "end"])
    result_df = pd.DataFrame({
    'freq_mean': grouped_df['freq'].mean(), # This is to get the mean of the freq column
    'freq_min': grouped_df['freq'].min(), # This is to get the min of the freq column
    'freq_max': grouped_df['freq'].max(), # This is to get the max of the freq column
    'freq_median': grouped_df['freq'].median(), # This is to get the median of the freq column
    'window_size_mean': grouped_df['window_size'].mean(), # This is to get the mean of the window_size column
    'sum_n_mean': grouped_df['sum_n'].mean() # This is to get the mean of the sum_n column
    }).reset_index() # This is to reset the index
    return result_df

def write_to_file(df, output_file): # This function writes the dataframe to a file
    df.to_csv(output_file, sep='\t', index=False)


def wrapper_function(ref_assembly, species, window_size, dictionary_of_pds): # This function is a wrapper function
    files = gather_files(ref_assembly, species, window_size, dictionary_of_pds) # This is to gather the files
    df = load_in_datasets(files) # This is to load in the datasets
    result_df = summarise_df(df) # This is to summarise the dataframe
    return result_df # This is to return the result dataframe




# ref_assembly = 'Cebus_albifrons'
# species = 'albifrons'
# window_size = '100000'
# output_file = 'test_albifrons_100000.csv'

# files = gather_files(ref_assembly, species, window_size, dictionary_of_pds)
# df = load_in_datasets(files)
# result_df = summarise_df(df)

ref_assembly = sys.argv[1]
species = sys.argv[2]
window_size = sys.argv[3]
output_file = sys.argv[4]

dictionary_of_pds = generatedictionary('/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt')
result_df = wrapper_function(ref_assembly, species, window_size, dictionary_of_pds)
write_to_file(result_df, output_file)