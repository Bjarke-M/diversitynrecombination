checkpoint split_ingroup_vcf:
    input:
        "resources/vcf/ingroup/{ref}.snps.vcf.gz",
    output:
        directory("results/vcf/ingroup/{ref}/chunks"),
        touch("results/vcf/ingroup/{ref}/chunks.done")
    params:
        size=1000000
    conda:
        "envs/base.yaml"
    script:
        "scripts/split_vcf.py"

rule lift_ingroup:
    input:
        chain='resources/chain/{ref}.liftOver.gz',
        touch="results/vcf/ingroup/{ref}/chunks.done",
        ref="resources/ref/hg38.fasta",
        dict="resources/ref/hg38.dict"
    output:
        vcf='results/vcf/hg38/ingroup/{ref}/chunks/old_ref/{i}.vcf.gz',
        rejected='results/vcf/hg38/ingroup/{ref}/chunks/old_ref/{i}.rejected.vcf.gz'
    conda:
        'envs/gatk.yaml'
    resources:
        runtime=60 * 12,
        mem_mb=16000
    params:
        vcf='results/vcf/ingroup/{ref}/chunks/{i}.vcf.gz'
    shell:
        """
        gatk IndexFeatureFile -I {params.vcf}
        gatk LiftoverVcf --java-options "-Xmx15g" -I {params.vcf} -O {output.vcf} -C \
        {input.chain} -R {input.ref} --REJECT {output.rejected} --RECOVER_SWAPPED_REF_ALT true \
        --WRITE_ORIGINAL_ALLELES true --WRITE_ORIGINAL_POSITION true
        """


# merge lifted chunks
rule merge_lifted_chunks_ingroup:
    input:
        touch="results/vcf/ingroup/{ref}/chunks.done",# taking directory as input caused ChildIOException
        vcfs=lambda w: expand('results/vcf/hg38/ingroup/{ref}/chunks/old_ref/{i}.vcf.gz',i=glob_wildcards(f'results/vcf/ingroup/{w.ref}/chunks/{{i}}.vcf.gz').i,ref=w.ref),
        index=lambda w: expand('results/vcf/hg38/ingroup/{ref}/chunks/old_ref/{i}.vcf.gz.tbi',i=glob_wildcards(f'results/vcf/ingroup/{w.ref}/chunks/{{i}}.vcf.gz').i,ref=w.ref)
    output:
        'results/vcf/hg38/ingroup/{ref}.vcf.gz'
    conda:
        'envs/gatk.yaml'
    params:
        input_flags=lambda wildcards, input: ' '.join(f'-I {vcf}' for vcf in input.vcfs)
    shell:
        """
        gatk MergeVcfs {params.input_flags} -O {output}
        """