import pandas as pd
import os

# Define the input path
input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'
data_dir = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/'


#/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/Saguinus_midas/PD_0017/PD_0017_coverage_2000000.txt

def generate_input_files(ind_coverage_file,window_size):
    with open(ind_coverage_file, 'r') as f:
        next(f)
        list_of_files = []
        for line in f:
            PDGP_ID,Genus,Species,FROH,Sex,ref_assembly = line.strip().split(',')
            #if ref_assembly != '':
            if ref_assembly == 'Saguinus_midas': 
                file_path = os.path.join(data_dir, ref_assembly, PDGP_ID, f'{PDGP_ID}_coverage_{window_size}.txt')
                list_of_files.append(file_path)
    return list_of_files

#def get_ind_name_from_path():


def extract_data_from_files(list_of_files,window_size):
    result_df = pd.DataFrame(columns=['chra', 'starta', 'coverage', 'frequency'])
    for file in list_of_files:
        data = pd.read_csv(file, sep='\t', header=None, names=['chra', 'starta', 'enda', 'chrb', 'startb', 'stopb' ,'coverage'])
        # Group by 'chra' and 'starta', then sum the 'coverage' for each group
        grouped_data = data.groupby(['chra', 'starta'])['coverage'].sum().reset_index()
        # Calculate 'sum_of_cov / window_size'
        grouped_data['frequency'] = grouped_data['coverage'] / window_size
        # Append the grouped data to the result DataFrame
        result_df = pd.concat([result_df, grouped_data])
        print(result_df)




extract_data_from_files(['/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/Saguinus_midas/PD_0017/PD_0017_coverage_2000000.txt'])