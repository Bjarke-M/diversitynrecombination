# this workflow is for calculating the pi away from genes for each species

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
            if ref_assembly not in dictionary_of_species_n_pdid:
                dictionary_of_species_n_pdid[ref_assembly] = {species: [pdgp_id]}
            elif species not in dictionary_of_species_n_pdid[ref_assembly]:
                dictionary_of_species_n_pdid[ref_assembly][species] = [pdgp_id]
            else:
                dictionary_of_species_n_pdid[ref_assembly][species].append(pdgp_id)
    return dictionary_of_species_n_pdid 


def generate_output_paths(species_list, window_list, command):
    input_paths = []
    for ref_assembly in species_list:
        for species in dictionary_of_inds[ref_assembly]:
            for window_size in window_list:
                if command=='distance_from_genes':
                    input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/distance_to_genes/{ref_assembly}/{species}_{window_size}_distance_to_genes.bed')
    return input_paths


def check_pi_files_for_combine_species(wildcards):
    input_paths = []
    ref_assembly = wildcards.ref_assembly
    window_size = wildcards.window_size
    for species in dictionary_of_inds[ref_assembly]:
            input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi')
    return input_paths
def check_callable_fracs_files_for_combine_species(wildcards):
    input_paths = []
    ref_assembly = wildcards.ref_assembly
    window_size = wildcards.window_size
    for species in dictionary_of_inds[ref_assembly]:
            input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/callable_fraction/{ref_assembly}/{species}/{window_size}_callable_fraction.csv')
    return input_paths
def check_distance_files_for_combine_species(wildcards):
    input_paths = []
    ref_assembly = wildcards.ref_assembly
    window_size = wildcards.window_size
    for species in dictionary_of_inds[ref_assembly]:
            input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/distance_to_genes/{ref_assembly}/{species}_{window_size}_distance_to_genes.bed')
    return input_paths

def check_all_files_exits(list_of_species): #check if all files exist
    input_paths = []
    for ref_assembly in list_of_species:
        input_paths.append(f'/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/distance_to_genes/{ref_assembly}_10000_combined.csv')
    return input_paths

window_list = [10000]
# species_list = ['Macaca_mulatta']
species_list = ['Daubentonia_madagascariensis','Saguinus_midas','Chlorocebus_aethiops',
                'Rhinopithecus_roxellana','Mandrillus_sphinx','Macaca_mulatta','Loris_tardigradus','Callithrix_jacchus',
                'Pithecia_pithecia','Lemur_catta-Thomas','Gorilla_gorilla_gorilla','Aotus_nancymaae','Sapajus_apella',
                'Pongo_abelii','Cercocebus_atys','Galago_moholi','Nomascus_leucogenys','Colobus_guereza',
                'Cebus_albifrons','Cercopithecus_mitis','Pongo_pygmaeus','Pan_troglodytes','Erythrocebus_patas',
                'Microcebus_murinus','Atele_fusciceps','Otolemur_garnettii','Nycticebus_pygmaeus']
#Missing g.vcf.gz species: Theropithecus_gelada
#Missing chainfiles: Carlito_syrichta Propithecus_coquereli'
dictionary_of_inds= generatedictionary(species_list)
list_of_inputs_to_combine = check_all_files_exits(species_list)

rule all:
    input:
        expand(generate_output_paths(species_list, window_list=window_list, command='distance_from_genes')),
        expand('/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/distance_to_genes/{ref_assembly}_{window_size}_combined.csv',ref_assembly=species_list,window_size=window_list),
        '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/distance_to_genes/combined_all/all_distances_species_combined.csv'

rule distance_from_nearest_gene:
    input:
        windows = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/{ref_assembly}/{species}/nonpar/{species}_{window_size}.bed',
        genes = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/distance_to_genes/sorted_hg38_protein_coding.bed'
    output:
        distance_to_genes = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/distance_to_genes/{ref_assembly}/{species}_{window_size}_distance_to_genes.bed'
    conda:
        'env/bedtools.yaml'
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:
        '''
        bedtools closest -a {input.windows} -b {input.genes} -d > {output.distance_to_genes}
        '''
    
rule combine_pi_frac_distance_to_genes:
    input:
        check_pi_files_for_combine_species,
        check_callable_fracs_files_for_combine_species,
        check_distance_files_for_combine_species
    output:
        combined = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/distance_to_genes/{ref_assembly}_{window_size}_combined.csv'
    params:
        ref_assembly='{ref_assembly}',
        window_size='{window_size}'
    conda:
        'env/pandas.yaml'
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:
        '''
        python /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/distance_to_genes/summaraise_pi_callable_distance.py {params.ref_assembly} {params.window_size} {output}
        '''



rule combine_all_species:
    input:
        list_of_inputs_to_combine
    output:
        file='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/distance_to_genes/combined_all/all_distances_species_combined.csv'
    conda:
        'env/pandas.yaml'
    resources:
        mem_mb=32000,
        runtime=30
    shell:
        '''
        python /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/distance_to_genes/combine_all_dist_and_species.py {output.file}
        '''