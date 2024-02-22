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

def real_paths(pd_id):
    real_paths = []
    for chr in chromosomes:
        trace_file = f"/home/bjarkemp/primatediversity/data/gVCFs_baboons_16_03_2021/{pd_id}/output.raw.snps.indels.{chr}.genotyped.g.vcf.gz"
        if os.path.exists(trace_file):
            real_paths.append(trace_file)
    
    # Create the directory if it does not exist
    output_directory = f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/Baboons/{pd_id}/"
    os.makedirs(output_directory, exist_ok=True)
    
    # Create the file and write the file paths
    output_file = f"{output_directory}/{pd_id}_files.txt"
    with open(output_file, 'w') as file:
        file.write('\n'.join(real_paths))
    
    return output_file

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
               'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chrX']


for ind in dictionary_of_inds['Papio_anubis']:
    real_paths(ind)

