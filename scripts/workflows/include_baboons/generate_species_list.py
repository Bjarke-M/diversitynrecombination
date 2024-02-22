import os

input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'
dictionary_of_inds = {}
with open(input_file, 'r') as file:
    header = next(file)  # Read and skip the header line
    for line in file:
        fields = line.strip().split(',')
        pdgp_id, genus, species, froh, sex, ref_assembly = fields
        if not ref_assembly:  # Check if ref_assembly is empty
            ref_assembly = 'unknown'
        # Add the pdgp_id to the appropriate category in the dictionary
        if ref_assembly not in dictionary_of_inds:
            dictionary_of_inds[ref_assembly] = [pdgp_id]
        else:
            dictionary_of_inds[ref_assembly].append(pdgp_id)

def real_paths(ids:list):
    real_paths = []
    for pd_id in ids:
        trace_file = f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/Baboons/{pd_id}/{pd_id}_concat.vcf.gz"
        real_paths.append(trace_file)
    # Create the directory if it does not exist
    output_directory = f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/Baboons/"
    os.makedirs(output_directory, exist_ok=True)
    # Create the file and write the file paths
    output_file = f"{output_directory}/all_inds_papio_files.txt"
    with open(output_file, 'w') as file:
        file.write('\n'.join(real_paths))
    return output_file



real_paths(dictionary_of_inds['Papio_anubis'])

