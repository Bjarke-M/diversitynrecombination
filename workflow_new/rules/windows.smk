rule windowed_pi:
    input:
        species_specific_lifted_filtered_files=
    output:
        windowed_pi=
    params:
        window_size= lambda wildcards: 
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        envs/vcftools 
    shell:
    '''
    vcftools -windowbased pi 
    '''

rule window_singletons:
    input:
        windows=path to bed file,
        species_specific_lifted_filtered_files=
    output:
        singletons =
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        envs/base.yaml
    script:
    'scripts/singletons.py'



rule filter_synonymous:
    input:
        vcf=species_specific_lifted_filtered_files,
    output:
        pi_s=
    
    
rule window_ps:
    input:
        windows=,
        vcf_s =
    output:
        windowed_pis=
    params:
        window_size= lambda wildcards: 
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        envs/vcftools 
    shell:
    '''
    vcftools -windowbased pi 
    '''



rule filter_nonsynonymous:
    input:
        vcf=species_specific_lifted_filtered_files,
    output:
        pi_s=
    
    
rule window_pn:
    input:
        windows=,
        vcf_s )
    output:
        windowed_pis=
    params:
        window_size= lambda wildcards: 
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        envs/vcftools 
    shell:
    '''
    vcftools -windowbased pi 
    '''

rule window_recombination:
