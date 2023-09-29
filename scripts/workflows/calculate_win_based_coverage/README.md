the command used to run the plotting of coverage 
snakemake --use-conda --rerun-incomplete --cluster "sbatch -A primatediversity -t 00:30:00 -N 1 --mem 16g" -j 150