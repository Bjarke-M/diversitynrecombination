rule filter_monomorphic: # Rule for filtering monomorphic sites from the gVCFs
    input:
        gvcf =  '../data1/{species}/{species}_batch_{batch}_gt.gvcf.gz'# Input GVCF file
    output:
        vcf = "../data1/{species}/{species}_batch_{batch}_gt.vcf.gz"       # Output VCF file
    conda:
        "envs/bcftools.yaml"  # Environment with bcftools installed
    shell:
        """
        bcftools view \
            -v snps,indels \
            -Oz \
            {input.gvcf} \
            -o {output.vcf}
            
        # Index the output VCF
        bcftools index {output.vcf}
        """


rule combine_batches_per_species: # Rule for combining the differenct batches before liftover
    input:
        vcf = expand("../data1/{species}/{species}_batch_{batch}_gt.vcf.gz", batch=batch,
                                                                            species = lambda wildcards: species) 
    output:
        vcf = "../data1/{species}/{species}_gt.vcf.gz" # Output VCF file
    conda:
        "envs/bcftools.yaml"  # Environment with bcftools installed
    shell:
        """
        bcftools concat \
            -Oz \
            {input.vcf} \
            -o {output.vcf}
            
        # Index the output VCF
        bcftools index {output.vcf}
        """

