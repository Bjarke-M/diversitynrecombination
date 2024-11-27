




rule all:
    input:
        expand(generatedictionary(input_file))

rule species_specific_bcfs: #creates species specific bcfs, by extracting the species in the multi species bcf
    input:
        bcfs = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/Baboons/Papio_anubis/Papio_anubis_filtered_monomorphic_lifted_sorted.vcf.gz'
    output:
        species_specific = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/original/{species}_only.bcf.gz",
        isolation_file = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/{species}_isolation.txt",
        csi = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/original/{species}_only.bcf.gz.csi"
    params:
        ref = "{ref_assembly}",
        speciesname = lambda wildcards: wildcards.species
    resources:
        runtime= 600, #in minutes
        mem_mb= 16000 #in megabytes
    conda:
        'envs/bcftools.yaml'
    shell:
        """
        python3.10 /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/species_specific_bcfs/isolation_file.py {params.ref} {params.speciesname} {output.isolation_file}
        bcftools view -S {output.isolation_file} {input.bcfs} -O b -o {output.species_specific}
        bcftools index {output.species_specific}
        """


# rule bcf_without_par:  #removes the PAR region from the species specific bcfs 
#     input:
#         bcfs = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/original/{species}_only.bcf.gz"
#     output:
#         nonPAR_region = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/{species}_nonpar.bcf.gz",
#     params:
#         chromosome = 'chrX',
#         PAR_start = 2816500 #END Of XG gene in hg38
#     conda:
#         'envs/bcftools.yaml'
#     resources:
#         runtime= 600, #in minutes
#         mem_mb= 16000 #in megabytes
#     shell:
#         """
#         bcftools view -t ^{params.chromosome}:1-{params.PAR_start} -O b {input.bcfs} > {output.nonPAR_region}
#         bcftools index {output.nonPAR_region}
#         """


# rule bcf_with_par: # makes a bcf with only the PAR region of all individuals 
#     input:
#         bcfs = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/original/{species}_only.bcf.gz"
#     output:
#         PAR_region = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/par/{species}_par.bcf.gz",
#     params:
#         chromosome = 'chrX',
#         PAR_start = 2816500 #END Of XG gene in hg38
#     conda:
#         'envs/bcftools.yaml'
#     resources:
#         runtime= 600, #in minutes
#         mem_mb= 16000 #in megabytes
#     shell:
#         """
#         bcftools view -t {params.chromosome}:1-{params.PAR_start} -O b {input.bcfs} > {output.PAR_region}
#         bcftools index {output.PAR_region}
#         """

# rule extract_males: #to remove any false calls of heterozygosity from the male X chr to ensure the pi estimation to be correct, this creates a list of the males in the species
#     output:
#         sex_isolation_file = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/males_no_x/male_{species}_isolation.txt"
#     params:
#         ref = lambda wildcards: wildcards.ref_assembly,
#         speciesname = lambda wildcards: wildcards.species
#     resources:
#         runtime= 10, #in minutes
#         mem_mb= 1000 #in megabytes
#     shell:
#         """
#         python3.10 /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/species_specific_bcfs/male_isoloation_file.py {params.ref} {params.speciesname} {output.sex_isolation_file}
#         """

# rule split_sexes_remove_male_X: # uses the txt file from above to remove the males from the species specific bcfs 
#     input:
#         bcfs = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/{species}_nonpar.bcf.gz",
#         sex_isolation_file = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/males_no_x/male_{species}_isolation.txt"
#     output:
#         females = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/females/{species}.bcf.gz",
#         males = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/males_no_x/{species}.bcf.gz"
#     conda:
#         'envs/bcftools.yaml'
#     resources:
#         runtime= 600, #in minutes
#         mem_mb= 16000 #in megabytes
#     shell:
#         """
#         bcftools view -S ^{input.sex_isolation_file} {input.bcfs} -O b > {output.females}
#         bcftools index {output.females}
#         bcftools view -S {input.sex_isolation_file} {input.bcfs} -Ou | bcftools view -t ^chrX -O b > {output.males}
#         bcftools index {output.males}
#         """


# rule merge_sexes: #merge the females with the males without the X chr into one file, which is to be used in estimating pi
#     input:
#         females = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/females/{species}.bcf.gz",
#         males = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/males_no_x/{species}.bcf.gz"
#     output:
#         merged="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/merged_non_male_X/{species}.bcf.gz"
#     conda: 
#         'envs/bcftools.yaml'
#     resources:
#         runtime= 800, #in minutes
#         mem_mb= 32000 #in megabytes
#     shell:
#         """
#         bcftools merge {input.females} {input.males} -O b > {output.merged}
#         bcftools index {output.merged}
#         """

# rule windowed_pi_for_each_species:
#     input:
#         nonpar = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/nonpar/merged_non_male_X/{species}.bcf.gz",
#         par = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/species_specific_bcfs/{ref_assembly}/{species}/par/{species}_par.bcf.gz"
#     output:
#         windowed_pi = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi",
#         windowed_pi_par = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/par/{species}_{window_size}_par.windowed.pi"
#     params:
#         window = "{window_size}",
#         prefix = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}",
#         prefix_par = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/par/{species}_{window_size}_par"
#     log:
#         nonpar='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/pi_estimation/windowed_pi/log/{ref_assembly}/{species}_{window_size}.log',
#         par='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/pi_estimation/windowed_pi/log/{ref_assembly}/{species}_{window_size}_par.log'
#     conda:
#         'envs/vcftools.yaml'
#     resources:
#         mem_mb=16000, #memory in megabytes
#         runtime=300 #runtime in minutes
#     shell:
#         """
#         vcftools --bcf {input.nonpar} --window-pi {params.window} --out {params.prefix} 2> {log.nonpar}
#         vcftools --bcf {input.par} --window-pi {params.window} --out {params.prefix_par} 2> {log.par}
#         """

# rule make_window_file:
#     input:
#         window_pi = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/nonpar/{species}_{window_size}.windowed.pi",
#         window_par = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/windowed_pi/{ref_assembly}/{species}/par/{species}_{window_size}_par.windowed.pi"
#     output:
#         window_bed = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/{ref_assembly}/{species}/nonpar/{species}_{window_size}.bed",
#         window_bed_par = "/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/{ref_assembly}/{species}/par/{species}_{window_size}_par.bed"
#     resources:
#         mem_mb=16000, #memory in megabytes
#         runtime=200 #runtime in minutes
#     shell:
#         """
#         #get the start end and chr and make bed file maybe
#         awk 'NR > 1 {{print $1 "\t" $2 "\t" $3}}' {input.window_pi} > {output.window_bed}
#         awk 'NR > 1 {{print $1 "\t" $2 "\t" $3}}' {input.window_par} > {output.window_bed_par}
#         """








# # Rule to generate individual BED files based on specified criteria
# rule generate_ind_beds:
#     input:
#         vcf="/home/bjarkemp/primatediversity/data/het_data_11_04_2022/{pdgp_id}/{pdgp_id}_concat.vcf.gz",
#     output:
#         output_file="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}.bed",
#         modcov_file="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}_modcov.txt"
#     params:
#         MIN_HET_AD=3,  # Adjust this value as needed
#         GQ=30  # Adjust this value as needed
#     shell:
#         """
#         modcov=$(bcftools stats -d 2,500,1 {input.vcf} | grep 'DP' | \
#         grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk '{{print $3}}')
#         echo $modcov > {output.modcov_file}
#         min_cov=$((modcov / 2))
#         max_cov=$((modcov * 2))
#         # Define the output file path for the current individual
#         # Apply filters and process variants, output to the defined file path
#         bcftools view {input.vcf} | \
#             bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < {params.MIN_HET_AD} ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= {params.GQ}" | \
#             grep -v '#' | \
#             awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
#             bedtools merge | \
#             sort -k1,1 -k2,2n | \
#             bedtools merge > {output.output_file}
#         """

# #lift bed files to the human reference genome
# rule liftover_generate_bed:
#     input:
#         chain_file="/home/bjarkemp/primatediversity/data/chain_files_15_03_2022/{species}_To_hg38.liftOver.gz",
#         input_bed= '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}.bed'
#     output:
#         bed="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}_bed_lifted.bed",
#         unmapbed="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}_bed_lifted.bed.unmap"
#     shell:
#         """
#         # Run CrossMap.py to perform liftover and generate VCF files
#         CrossMap.py bed --chromid l {input.chain_file} {input.input_bed} {output.bed} --unmap-file {output.unmapbed}
#         """



# rule sort_lifted:
#     input:
#         input_bed= '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}_bed_lifted.bed'
#     output:
#         bed="/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{species}/{pdgp_id}/{pdgp_id}_lifted_sorted.bed"
#     shell:
#         """
#         sort -k1,1 -k2,2n {input.input_bed} > {output.bed}
#         """





# ## intersect callable and liftable parts ##
# rule intersect_windows:
#     input:
#         window_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/window_bed/{ref_assembly}/{species}/nonpar/{species}_{window_size}.bed',        
#         callable_vcf = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{ref_assembly}/{pd_id}/{pd_id}_lifted_sorted.bed'
#     output:
#         outfile = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{ref_assembly}/{species}/{pd_id}/nonpar/txt/{pd_id}_coverage_{window_size}.txt'
#     conda:
#         'envs/bedtools.yaml'
#     resources:
#         mem_mb=64000,
#         runtime=30
#     shell:
#         '''
#         bedtools intersect -a {input.window_file} -b {input.callable_vcf} -wo > {output.outfile}
#         '''


# ## calculate coverage per individual ##
# rule maskfile:
#     input:
#         coverage_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{ref_assembly}/{species}/{pd_id}/nonpar/txt/{pd_id}_coverage_{window_size}.txt'
#     output:
#         masked_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/{ref_assembly}/{species}/{pd_id}/nonpar/csv/{pd_id}_masked_{window_size}.csv'
#     params:
#         window_size = "{window_size}",
#         pd_id = "{pd_id}",
#         sex = lambda wildcards: get_sex(wildcards.pd_id) #NA for coverrage if sex = male and chr = X as we dont want to make the mistake of including the male x coverage when we summarise for species
#     conda:
#         'envs/maskfile.yaml'
#     resources:
#         mem_mb=64000,
#         runtime=30
#     script: '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/calculate_win_based_coverage/maskfile.R'


# # rule summarise_callable_fraction:
# #     input:
# #         expand(generate_input_paths(species_list,window_list)) #checks wheter all files exist not just the specific window_size it will generate the output for in each iteration
# #     output:
# #         output='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/callable_fraction/{ref_assembly}/{species}/{window_size}_callable_fraction.csv'
# #     params:
# #         ref_assembly='{ref_assembly}',
# #         species='{species}',
# #         window_size='{window_size}'
# #     resources:
# #         mem_mb=10000,
# #         runtime=30
# #     conda:
# #         'envs/pandas.yaml'
# #     shell:
# #         '''
# #         python /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/summarise_species/summarise_species_callable_fraction.py {params.ref_assembly} {params.species} {params.window_size} {output}
# #         ''' 



# rule combine_species:
#     input:
#         check_pi_files_for_combine_species,
#         check_callable_fracs_files_for_combine_species
#     output:
#         output='/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/results/combined/{ref_assembly}/{ref_assembly}_{window_size}_combined.csv'
#     params:
#         ref_assembly='{ref_assembly}',
#         window_size='{window_size}'
#     conda:
#         'envs/pandas.yaml'
#     resources:
#         mem_mb=10000,
#         runtime=30
#     shell:
#         '''
#         python /home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/scripts/workflows/summarise_species/summarise_pi_n_frac.py {params.ref_assembly} {params.window_size} {output}
#         '''
    
