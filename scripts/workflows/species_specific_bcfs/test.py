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
            dictionary_of_inds[ref_assembly] = {species: [pdgp_id]}
        elif species not in dictionary_of_inds[ref_assembly]:
            dictionary_of_inds[ref_assembly][species] = [pdgp_id]
        else:
            dictionary_of_inds[ref_assembly][species].append(pdgp_id)

# def get_output_paths(dictionary_of_inds, ref_assembly, species):
#     return [
#         f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}_only.bcf.gz"
#         for species in dictionary_of_inds[ref_assembly]
#     ]
def get_output_paths(dictionary_of_inds):
    out_paths = []
    for ref_assembly in dictionary_of_inds:
        #print(ref_assembly)
        for species in dictionary_of_inds[ref_assembly]:
            #print(ref_assembly,species)
            out_paths.extend([
                    f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}_only.bcf.gz"
                    for species in dictionary_of_inds[ref_assembly]
                ])
    return out_paths


def get_pd_ids_from_species(dictionary_of_inds, ref_assembly, species):
    return ','.join(dictionary_of_inds[ref_assembly][species])

print(get_pd_ids_from_species(dictionary_of_inds, 'Saguinus_midas', 'mystax'))
