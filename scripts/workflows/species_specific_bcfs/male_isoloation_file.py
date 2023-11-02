import sys

input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'

def generatedictionary(input):
    dictionary_of_species_n_pdid = {}
    with open(input_file, 'r') as file:
        header = next(file)  # Read and skip the header line
        for line in file:
            fields = line.strip().split(',')
            pdgp_id, genus, species, froh, sex, ref_assembly = fields
            if not ref_assembly:  # Check if ref_assembly is empty
                ref_assembly = 'unknown'
            # Add the pdgp_id to the appropriate category in the dictionary
            elif ref_assembly not in dictionary_of_species_n_pdid:
                dictionary_of_species_n_pdid[ref_assembly]={species:{sex:[pdgp_id]}}
            elif species not in dictionary_of_species_n_pdid[ref_assembly]:
                dictionary_of_species_n_pdid[ref_assembly][species] = {sex: [pdgp_id]}
            elif sex not in dictionary_of_species_n_pdid[ref_assembly][species]:
                dictionary_of_species_n_pdid[ref_assembly][species][sex] = [pdgp_id]
            else:
                dictionary_of_species_n_pdid[ref_assembly][species][sex].append(pdgp_id)
    return dictionary_of_species_n_pdid 


dictionary_of_inds = generatedictionary(input_file)

#print(dictionary_of_inds)
def get_pd_ids_from_species(dictionary_of_species_n_pdid, ref_assembly, species, outfile):
    try:
        result = '\n'.join(dictionary_of_species_n_pdid[ref_assembly][species]['M'])
    except KeyError:
        result = ''
    with open(outfile, 'w') as file:
        file.write(result)
    return result  



ref = sys.argv[1]
species = sys.argv[2]
outfile = sys.argv[3]

get_pd_ids_from_species(dictionary_of_inds, ref, species, outfile)