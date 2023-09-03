Here is a workflow that should process the vcf files and perform liftover to hg38 using CrossMap, from all the primate species.

NOTE it has not been run yet

#snakemake --cluster "sbatch -A primatediversity -t 10:00:00 -N 1 --mem 16g" -j 2