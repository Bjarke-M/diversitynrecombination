'''
------------------------------------------------------------------------------------------------------------------------
Preparing the gvcfs, by  generating a callability mask for every individual basedon the gvcf.
Then, lifting it to human coordinates, both for a jointcalled vcf and the callability masks.
------------------------------------------------------------------------------------------------------------------------
'''

import os
import glob
import pandas as pd
from gwf import Workflow

# Absolute path to the input data (gvcfs)
vcfs_path = "/faststorage/project/primatediversity/data/gVCF_jointCalling_17_05_2021/"
gvcf_path = "/faststorage/project/primatediversity/data/gVCF_transfer_5.5.2021/"
outpath_het = "/faststorage/project/primatediversity/data/het_data_11_04_2022/"
outpath_vcf_lift = outpath_het + "lifted_GtVar/"
outpath_zarr_lift = outpath_het + "lifted_Zarr/"
outpath_zarr_ref = outpath_het + "ref_Zarr/"
human_fa = "/faststorage/project/primatediversity/data/chain_files_15_03_2022/hg38/hg38.fasta"
chain_path = "/faststorage/project/primatediversity/data/chain_files_15_03_2022"
# bed_sizes = "/home/eriks/primatediversity/people/erik/Primate_Het_X_Autosome/data/liftover_beds/"

gwf = Workflow(defaults={"account": "primatediversity"})


########################################################################################################################
############################################### ---- VCF CONCAT ---- ###################################################
########################################################################################################################


def get_ID_concat(idx, target):
    ID = target[1][0].split("/")[6]
    return '{}_concat'.format(ID)


def gvcf_concat(file_list, ID, outpath):
    """Concatenating all the gvcf chunks and then tabix indexing"""
    o = outpath + "/" + ID + "_concat.vcf.gz"
    inputs = file_list
    outputs = [o]
    options = {'cores': 1, 'memory': "8g", 'walltime': "07:00:00"}
    pd.Series(file_list).to_csv(outpath + "/file_list.txt", index=False, header=False)
    spec = """
    bcftools concat -O z -n -f {file_list} -o {outfile}_temp
    mv {outfile}_temp {outfile}
    """.format(file_list=outpath + "/file_list.txt", outfile=o)
    return (inputs, outputs, options, spec)


########################################################################################################################
############################################# ---- CALLABILITY MASK ---- ###############################################
########################################################################################################################


def get_ID_callmask(idx, target):
    ID = target[1][0].split("/")[6][:14]
    return '{}_callmask'.format(ID)


def callMask(file_list, ID, outpath, min_het, gq):
    """Reference callability."""
    i = outpath + "/" + ID + "_concat.vcf.gz"
    inputs = [i]
    outfile = outpath + "/" + ID 
    outputs = [outfile + ".bed"]
    options = {'cores': 1, 'memory': "8g", 'walltime': "20:00:00", "account": 'primatediversity'}
    spec = """
    tabix {gvcf} -f

    bcftools stats -d 2,500,1 {gvcf} | grep 'DP' |\
    grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk {awk} > {stat_path}_modcov.txt

    modcov=$(<{stat_path}_modcov.txt)
    min_cov=$((modcov/2))
    max_cov=$((modcov*2))
    echo $max_cov
    echo $min_cov

    bcftools view -Ou {gvcf} | bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < {MIN_HET_AD} ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= {GQ} " {gvcf} |\
    grep -v '#' | \
    awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
    bedtools merge | \
    sort -k1,1 -k2,2n | \
    bedtools merge > {outfile}_temp

    mv {outfile}_temp {outfile}.bed
    """.format(gvcf=i, awk="'{print $3}'",
               outfile=outfile, stat_path=outfile,
               MIN_HET_AD=min_het, GQ=gq)
    return (inputs, outputs, options, spec)


########################################################################################################################
############################################# ---- CALLABILITY MASK ---- ###############################################
########################################################################################################################


def get_ID_vcf_lift(idx, target):
    ID = target[1][0].split("/")[-1][:10]
    return '{}_vcf_lift'.format(ID)


def vcf_lift(ref, outpath, jointcall_path, chain_path, human_fa):
    """Lifting to human coordinates"""
    ref_name = ref.split("/")[-1][:-20]
    ID_path = outpath  + ref_name
    inputs = [ref]
    outputs = [ID_path + "_lifted.vcf.gz"]
    options = {'cores': 2, 'memory': "16g", 'walltime': "60:00:00"}
    spec = """
    less {ref} | head -n 500 |  grep 'PGDP/REFERENCES_FROZEN' | \
    awk '{{print $11}}' | awk -F/ '{{print $8}}' | awk -F. '{{print $1}}' | xargs > {ID_path}_ref.txt
    ref_genome=$(<{ID_path}_ref.txt)
    echo {chain_path}/${{ref_genome}}_To_hg38.liftOver.gz

    CrossMap.py vcf {chain_path}/${{ref_genome}}_To_hg38.liftOver.gz {ref} {human_fa} {ID_path}_lifted.vcf --no-comp-alleles --chromid l
    bgzip {ID_path}_lifted.vcf
    bgzip {ID_path}_lifted.vcf.unmap
    """.format(ref=ref, ID_path=ID_path, jointcall_path=jointcall_path, chain_path=chain_path, human_fa=human_fa)
    return (inputs, outputs, options, spec)


def get_ID_vcf_to_zarr(idx, target):
    ID = target[0].split("/")[-1][:10]
    return '{}_vcf_to_zarr'.format(ID)


def vcf_to_zarr(ref, outpath):
    ref_name = ref.split("/")[-1][:-20]
    output = outpath + ref_name
    short_name = ref[:-6]
    inputs = ref
    outputs = output + "/calldata/"
    options = {'memory': '50g',
               'walltime': '60:00:00'}
    spec = """
    bcftools sort {ref} -o {ref}
    tabix {ref}
    python scripts/vcf_to_zarr.py {ref} {output}
    """.format(ref=ref,  short_name=short_name, output=output)
    print(spec)
    return (inputs, outputs, options, spec)


def get_ID_lifted_vcf_to_zarr(idx, target):
    ID = target[0].split("/")[-1][:10]
    return '{}_lifted_vcf_to_zarr'.format(ID)


def lifted_vcf_to_zarr(ref, outpath):
    ref_name = ref.split("/")[-1][:-14]
    plain_ref = ref[:-3]
    output = outpath + ref_name
    inputs = ref
    outputs = output + "/calldata/"
    options = {'memory': '15g',
               'walltime': '24:00:00'}
    # bcftools sort {ref} -o {ref}
    # tabix {ref}
    spec = """
    python scripts/vcf_to_zarr.py {ref} {output}
    """.format(ref=ref, plain_ref=plain_ref, output=output)
    return (inputs, outputs, options, spec)


def get_ID_compress(idx, target):
    ID = target[1][0].split("/")[6]
    return '{}_compress'.format(ID)


def vcf_compress(file_list, ID, outpath):
    """Compressing vcfs"""
    i = outpath + "/" + ID + "_concat.vcf.gz"
    ID_path = outpath + "/" + ID
    inputs = [ID_path + "_lifted.vcf"]
    outputs = [ID_path + "_lifted.vcf.gz"]
    options = {'cores': 1, 'memory': "8g", 'walltime': "10:00:00"}
    spec = """
    gzip {ID_path}_lifted.vcf
    gzip {ID_path}_lifted.vcf.unmap
    """.format(ID_path=ID_path)
    print(spec)
    return (inputs, outputs, options, spec)


def get_ID_bed_lift(idx, target):
    ID = target[1][0].split("/")[6]
    return '{}_bed_lift'.format(ID)


def bed_lift(file_list, ID, outpath, chain_path):
    """Lifting to human coordinates"""
    i = outpath + "/" + ID + "_concat.vcf.gz"
    ID_out = outpath + "/" + ID
    inputs = [i, outpath + "/" + ID + ".bed"]
    ID_out = outpath + "/" + ID
    outputs = [ID_out + "_mapped.bed"]
    options = {'cores': 2, 'memory': "16g", 'walltime': "24:00:00"}
    spec = """
    less {ID}_concat.vcf.gz | head -n 500 |  grep 'PGDP/REFERENCES_FROZEN' | \
    awk '{{print $11}}' | awk -F/ '{{print $8}}' | awk -F. '{{print $1}}' | xargs > {ID}_ref.txt
    ref_genome=$(<{ID}_ref.txt)
    echo $ref_genome

    bedtools intersect -a {ID}.bed -b {chain_path}/${{ref_genome}}_To_hg38.liftOver.bed -sorted > {ID}_intersect.bed

    CrossMap.py bed {chain_path}/${{ref_genome}}_To_hg38.liftOver.gz {ID}_intersect.bed \
    {ID}_mapped.bed --unmap-file {ID}_unmapped.bed
    """.format(ID=ID_out, chain_path=chain_path)
    return (inputs, outputs, options, spec)


def bed_sort(file_list, ID, outpath):
    """Sorting the lifted bed files"""
    ID_path = outpath + "/" + ID
    inputs = [ID_path + "_mapped.bed"]
    outputs = [ID_path + "_mapped_sorted.bed"]
    options = {'cores': 1, 'memory': "8g", 'walltime': "2:00:00"}
    spec = """
    sort -k1,1 -k2,2n {file_name}_mapped.bed | \
    bedtools merge > {file_name}_mapped_sorted.bed 
    """.format(file_name = ID_path)
    return (inputs, outputs, options, spec)


########################################################################################################################
################################################ ---- RUN PIPELINE ---- ################################################
########################################################################################################################


### get individual callability with Lukas protocol


os.makedirs("steps/depth_stats", exist_ok=True)
os.makedirs(outpath_het, exist_ok=True)
os.makedirs(outpath_vcf_lift, exist_ok=True)
os.makedirs(outpath_zarr_lift, exist_ok=True)

### Four possible gvcf_paths.

gvcf_path = "/faststorage/project/primatediversity/data/gVCF_transfer_5.5.2021/" # Primary
# gvcf_path = "/faststorage/project/primatediversity/data/gVCFs_02_03_2021/" #Second
# gvcf_path = "/faststorage/project/primatediversity/data/gVCF_transfer_np_remap_28.5.2021/" #Third

ID_paths = glob.glob(gvcf_path+"*")
ref_ID_inputs = []
ref_ID_callmask = []
for PD_path in ID_paths:
    PD_ID = PD_path.split("/")[-1]
    d_call = {}
    i, x, path_list = 0, True, []
    ref_path = "{}/{}/{}_ref.txt".format(outpath_het, PD_ID, PD_ID)
    ref_name = "?"
    if os.path.exists(ref_path):
        ref_name = pd.read_csv(ref_path, names=["ref"]).iat[0, 0]
    if ref_name == "Propithecus_coquereli" or ref_name == "Carlito_syrichta":
        continue
    while x == True:
        p  = gvcf_path + "{}/{}_{}.raw.snps.indels.g.vcf.gz".format(PD_ID, PD_ID, i)
        if os.path.exists(p):
            path_list.append(p)
            i += 1
        else:
            x = False
    #Manual removal of two species in which I do not have a liftover file.
    if len(path_list) == 0:
        pass
        # If you want a warning for missing inds - they are in the remap folder ir 02_03
        # print("There is no data for ", PD_ID)
    else:
        d_call["file_list"] = path_list
        d_call["ID"] = PD_ID
        d_call["outpath"] = outpath_het + PD_ID
        ref_ID_callmask.append(d_call)
    os.makedirs(outpath_het + "/{}".format(PD_ID), exist_ok=True)

ref_ID_callmask = ref_ID_callmask

### To generate the callmasks, I am going to merge all the chunks for each raw individual.

### Section 1: concat gvcf and generate callmask bed

# concat_gwf = gwf.map(gvcf_concat, ref_ID_callmask, name=get_ID_concat)

# callmask_gwf = gwf.map(callMask, ref_ID_callmask, name=get_ID_callmask, extra={"min_het": 3, "gq": 30})

### Section 2: lift variant vcf and convert to zarr, as well as the original refs.
vcfs_path = "/home/eriks/primatediversity/data/gVCF_jointCalling_17_05_2021/Atele_fusciceps/" #Alternative path due to permission problem
GtVar_l = glob.glob(vcfs_path + "*ConcatGtVar.g.vcf.gz")  
print(GtVar_l)

vcf_to_zarr_output = gwf.map(vcf_to_zarr, GtVar_l, name=get_ID_vcf_to_zarr,
                             extra={"outpath": outpath_zarr_ref})

# vcf_lift_gwf = gwf.map(vcf_lift, GtVar_l, name=get_ID_vcf_lift,
#                         extra={"outpath": outpath_vcf_lift, "jointcall_path": vcfs_path,
#                                 "chain_path": chain_path, "human_fa": human_fa})

# lifted_vcf_to_zarr_output = gwf.map(lifted_vcf_to_zarr, vcf_lift_gwf.outputs, name=get_ID_lifted_vcf_to_zarr,
#                              extra={"outpath": outpath_zarr_lift})



## The compress function was only ran because I had made a mistake in which the compress step was
## not included in the vcf_lift job
## vcf_compress_gwf = gwf.map(vcf_compress, ref_ID_callmask[0:100], name=get_ID_compress)

### Section 3: Lift bed and sort it.

# print(ref_ID_callmask)

# bed_lift_gwf = gwf.map(bed_lift, ref_ID_callmask, name=get_ID_bed_lift,
#                        extra={"chain_path": chain_path})

# bed_sort_gwf = gwf.map(bed_sort, ref_ID_callmask)
