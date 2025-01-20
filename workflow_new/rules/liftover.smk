rule lift_species_gatk:
    input:
        vcf= '../data1/{species}/{species}_gt.vcf.gz',
        chain='resources/chain/{ref}.liftOver.gz',
        ref="resources/ref/hg38.fasta",
    output:
        vcf='../data1/liftover/gatk/{ref}/{species}_liftover.vcf.gz',
        rejected='../data1/liftover/gatk/{ref}/{species}_liftover.vcf.gz.rejected'
    conda:
        'envs/gatk.yaml'
    resources:
        runtime=60 * 12,
        mem_mb=16000
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk LiftoverVcf --java-options "-Xmx15g" -I {input.vcf} -O {output.vcf} -C \
        {input.chain} -R {input.ref} --REJECT {output.rejected} --RECOVER_SWAPPED_REF_ALT true \
        --WRITE_ORIGINAL_ALLELES true --WRITE_ORIGINAL_POSITION true
        """


rule lift_species_crossmap:
    input:
        vcf= '../data1/{species}/{species}_gt.vcf.gz',
        chain='resources/chain/{ref}.liftOver.gz',
        ref="resources/ref/hg38.fasta",
    output:
        vcf='../data1/liftover/crossmap/{ref}/{species}_liftover.vcf.gz',
        rejected='../data1/liftover/crossmap/{ref}/{species}_liftover.vcf.gz.rejected'
    conda:
        'envs/crossmap.yaml'
    resources:
        runtime=60 * 12,
        mem_mb=16000
    shell:
        """
        # Run CrossMap.py to perform liftover and generate VCF files
        CrossMap.py vcf --chromid l --no-comp-alleles {input.chain_file} {input.input_vcf} {input.reference_fasta} {output.vcf}
        """
