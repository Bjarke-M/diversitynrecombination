The folder includes the workflow used to make coverage plots for a lot of different window sizes across species and individuals.

the command used to run the plotting of coverage 
snakemake --use-conda --rerun-incomplete --cluster "sbatch -A primatediversity -t 00:30:00 -N 1 --mem 16g" -j 150