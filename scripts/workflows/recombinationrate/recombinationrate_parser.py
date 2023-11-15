import pandas as pd
import numpy as np
import sys
#get a window file and an input file with recombination rates and calculate the average recombination rate in each window
# file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/Cercopithecus_mitis/ascanius/nonpar/ascanius_100000.bed'
#recomb_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/Recombination_map_decode_2019'

# # get files from command line
file = sys.argv[1]
recomb_file = sys.argv[2]
output_file = sys.argv[3]

#parse files to pandas dataframes
recomb_file = pd.read_table(recomb_file, sep='\t', names=['chr', 'start', 'end', 'cmpermb', 'cm'])
df = pd.read_table(file, sep='\t', names=['chr', 'start', 'end'])

#find the previous end in the recomb file
def find_previous_end(chr, start):
    return recomb_file.loc[(recomb_file['chr'] == chr) & (recomb_file['end'] <= start), 'end'].max()

#find the next end in the recomb file
def find_next_end(chr, start):
    return recomb_file.loc[(recomb_file['chr'] == chr) & (recomb_file['end'] >= start), 'end'].min()

#returns nan if the either the previous or next end is nan, which means we only look at windows that are within the recomb map from DECODE
# to avoid extrapolating across somethings we do not understand (e.g. the recombination rate in the telomeres)
def calculate_relative_position(start, chr):
    relevant_rows = recomb_file[(recomb_file['chr'] == chr) & (recomb_file['end'] <= start)] # get all rows where the chromosome is the same as the window and the end is smaller than the start of the window
    
    previous_end = relevant_rows['end'].max() # get the largest end value
    cm_previous_end = relevant_rows.loc[relevant_rows['end'] == previous_end, 'cm'].max() # get the cm value for the previous end
    
    relevant_rows = recomb_file[(recomb_file['chr'] == chr) & (recomb_file['end'] >= start)] # get all rows where the chromosome is the same as the window and the end is larger than the start of the window
    
    next_end = relevant_rows['end'].min() # get the smallest end value
    cm_next_end = relevant_rows.loc[relevant_rows['end'] == next_end, 'cm'].max() # get the cm value for the next end

    cm_per_bp = ((cm_next_end - cm_previous_end) / (next_end - previous_end)) # calculate the cm per bp
    bp_from_end_to_start = (start - previous_end) # calculate the number of bp from the previous end to the start of the window
    
    cm_start = cm_previous_end + (cm_per_bp * bp_from_end_to_start)
    
    if chr == 'chrX':
        cm_start *= (2/3)
    return cm_start




def calculate_cm_per_mb_per_window(start, end, chr):
    cm_start = calculate_relative_position(start, chr)
    cm_end = calculate_relative_position(end, chr)
    cm_per_mb = ((cm_end - cm_start) / (end - start)) * 1000000
    return cm_start, cm_end, cm_per_mb

with open(output_file, 'w') as output:
    output.write('chr\tstart\tend\tcm_start\tcm_end\tcm_per_mb\n')
    for chr, start, end in df.itertuples(index=False):
        cm_start, cm_end, cm_per_mb = calculate_cm_per_mb_per_window(start, end, chr)
        output_line = f'{chr}\t{start}\t{end}\t{cm_start}\t{cm_end}\t{cm_per_mb}\n'
        output.write(output_line)
    output.close()
