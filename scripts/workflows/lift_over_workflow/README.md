Here is a workflow that should process the vcf files and perform liftover to hg38 using CrossMap, from all the primate species.

NOTE it worked and im quite satisfied

the command run on the cluster 
#snakemake --cluster "sbatch -A primatediversity -t 24:00:00 -N 1 --mem 16g" -j 3 