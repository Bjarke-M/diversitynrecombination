rule filter_gvcf:


rule liftover_gvcf:

rule merge_liftovergvcf:

rule calculate_windowbased_callability:
    input: 
        windows = 
        merged_lifted_gvcf =
    output:
        callability_file = 
    resources:
        runtime=20,
        mem_mb=1000
    conda:
        envs/base.yaml
    script:
    "scripts/get_callability.py"


    


