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
            if ref_assembly not in dictionary_of_species_n_pdid:
                dictionary_of_species_n_pdid[ref_assembly] = {species: [pdgp_id]}
            elif species not in dictionary_of_species_n_pdid[ref_assembly]:
                dictionary_of_species_n_pdid[ref_assembly][species] = [pdgp_id]
            else:
                dictionary_of_species_n_pdid[ref_assembly][species].append(pdgp_id)
    return dictionary_of_species_n_pdid 

def get_output_paths(dictionary_of_species_n_pdid, window_list, command):
    out_paths = []
    for ref_assembly in dictionary_of_species_n_pdid:
        if ref_assembly != 'unknown' and ref_assembly in ref_list:
            for species in dictionary_of_species_n_pdid[ref_assembly]:
                if command=='windowpi':
                    for window_size in window_list:
                        out_paths.append(f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/windowed_pi/{ref_assembly}/{species}/{species}_{window_size}.windowed.pi")
                elif command=='windowbed':
                    for window_size in window_list:
                        out_paths.append(f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/{ref_assembly}/{species}/{species}_{window_size}.bed")
                elif command=='maskcov':
                    for window_size in window_list:
                        for pd_id in dictionary_of_species_n_pdid[ref_assembly][species]:
                            out_paths.append(f"/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{ref_assembly}/{species}/{pd_id}/txt/{pd_id}_coverage_{window_size}.txt")
    return out_paths


## LOAD DICTIONARY AND WINDOWLIST AND REFS 
# ref_list = ['Daubentonia_madagascariensis','Saguinus_midas','Chlorocebus_aethiops','Rhinopithecus_roxellana',
#                 'Mandrillus_sphinx','Macaca_mulatta','Loris_tardigradus','Callithrix_jacchus',
#                 'Pithecia_pithecia','Lemur_catta-Thomas','Gorilla_gorilla_gorilla','Aotus_nancymaae','Sapajus_apella',
#                 'Pongo_abelii','Cercocebus_atys','Galago_moholi','Nomascus_leucogenys','Colobus_guereza',
#                 'Cebus_albifrons','Cercopithecus_mitis','Pongo_pygmaeus','Pan_troglodytes','Erythrocebus_patas',
#                 'Microcebus_murinus','Atele_fusciceps','Otolemur_garnettii','Nycticebus_pygmaeus']
ref_list = ['Microcebus_murinus']
#Missing g.vcf.gz species: Theropithecus_gelada
#Missing chainfiles: Carlito_syrichta Propithecus_coquereli'
input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'
list_of_windows = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 2000000, 10000000]
dictionary_of_inds=generatedictionary(input_file)

#print(dictionary_of_inds)
print(get_output_paths(dictionary_of_inds,list_of_windows,'maskcov'))