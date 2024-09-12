import os
import sys
import pandas as pd


def load_in_datasets(files) -> pd.DataFrame: # This function loads in the datasets
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        dfs += [df]
    return pd.concat(dfs)

def write_to_file(df, output_file): # This function writes the dataframe to a file
    df.to_csv(output_file, sep='\t', index=False)



ref_list = ['Daubentonia_madagascariensis','Saguinus_midas','Chlorocebus_aethiops','Rhinopithecus_roxellana',
                'Mandrillus_sphinx','Macaca_mulatta','Loris_tardigradus','Callithrix_jacchus',
                'Pithecia_pithecia','Lemur_catta-Thomas','Gorilla_gorilla_gorilla','Aotus_nancymaae','Sapajus_apella',
                'Pongo_abelii','Cercocebus_atys','Galago_moholi','Nomascus_leucogenys','Colobus_guereza',
                'Cebus_albifrons','Cercopithecus_mitis','Pongo_pygmaeus','Pan_troglodytes','Erythrocebus_patas',
                'Microcebus_murinus','Atele_fusciceps','Otolemur_garnettii','Nycticebus_pygmaeus']

def check_all_files_exits(list_of_species): #check if all files exist
    input_paths = []
    for ref_assembly in list_of_species:
        input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/tajimasd/combined/{ref_assembly}/{ref_assembly}_100000_combined.csv')
    return input_paths
list_of = check_all_files_exits(ref_list)

# fileone = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/tajimasd/combined/Pithecia_pithecia/Pithecia_pithecia_100000_combined.csv'
# filetwo = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/tajimasd/combined/Lemur_catta-Thomas/Lemur_catta-Thomas_100000_combined.csv'

#list_of_files = sys.argv[1:]

df = load_in_datasets(list_of)
output=sys.argv[1]

write_to_file(df, output)