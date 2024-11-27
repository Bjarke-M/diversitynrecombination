rule filter_vcf:
    input:
        vcf=
    output:
        clean_vcf =
        logfile = 
    resources:
        runtime=10,
        mem_mb=1000
    conda:
        envs/bcftools
    shell:
    '''
    bcftools filter 
    '''

    
