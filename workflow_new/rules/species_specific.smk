'results/vcf/hg38/ingroup/{ref}.vcf.gz'

rule get_pdids:
    input: 
        metadata = #path to file 
    output:
        files = expand(path/to/isolationfiles_{species}.txt, species = species_ids)
    conda:
        envs/base.yaml
    resources:
        runtime=5,
        mem_mb=1000
    script:
        "scripts/get_pdids.py"

rule species_specific_files:
    input:
        pdids = path/to/isolationfiles_{species}.txt
        reference = path/to/{reference_genome}.vcf.gz
    output:
        {reference_genome}/{species}.bcf.gz
    conda:
        envs/bcftools.yaml
    resources:
        runtime=10,
        mem_mb=8000
    shell:
    '''
    bcftools view -Bz -S {input.pdids} {input.reference} 
    '''