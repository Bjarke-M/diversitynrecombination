rule windowed_pi:
    input:
        vcf = '../data1/liftover/{method}/{ref}/{species}_liftover.vcf.gz'
    output:
        windowed_pi= '../results1/windowed_stats/{method}/{ref}/{species}.windowed.pi'
    params:
        window_size= 10000
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        'envs/vcftools'
    shell:
    '''
    vcftools --gzvcf {input.vcf} --window-pi {params.window_size} --out {output.windowed_pi} 
    '''


# rule filter_synonymous:
#     input:
#         vcf=species_specific_lifted_filtered_files,
#     output:
#         pi_s=
    
    
# rule window_ps:
#     input:
#         windows=,
#         vcf_s =
#     output:
#         windowed_pis=
#     params:
#         window_size= lambda wildcards: 
#     resources:
#         runtime=10,
#         mem_mb=1000
#     conda:
#         envs/vcftools 
#     shell:
#     '''
#     vcftools -windowbased pi 
#     '''



# rule filter_nonsynonymous:
#     input:
#         vcf=species_specific_lifted_filtered_files,
#     output:
#         pi_s=
    
    
# rule window_pn:
#     input:
#         windows=,
#         vcf_s )
#     output:
#         windowed_pis=
#     params:
#         window_size= lambda wildcards: 
#     resources:
#         runtime=10,
#         mem_mb=1000
#     conda:
#         envs/vcftools 
#     shell:
#     '''
#     vcftools -windowbased pi 
#     '''

# rule window_recombination:
