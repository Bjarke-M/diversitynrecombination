#Snakefile for calculating proportions of snps per window for each species

################################################## Function Bonanza ##################################################
def generatedictionary(input): #make a dictionary of species and pdgp_ids, to iterate over
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








vcftools --bcf data/species_specific_bcfs/Daubentonia_madagascariensis/madagascariensis/nonpar/merged_non_male_X/madagascariensis.bcf.gz --singletons --out test_single_ind