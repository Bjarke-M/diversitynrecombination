# Workflow for estimating pi across windows using vcftools

This workflow uses vcftools to estimate pi across windows. The input file is a VCF file containing SNP data for a set of samples. The output is a file containing pi estimates for each window.

## Requirements

- vcftools
- a VCF file containing SNP data

## Usage

1. Install vcftools if it is not already installed.
2. Prepare a VCF file containing SNP data.
3. Modify the `vcftools_pi.sh` script to specify the input file and window size.
4. Run the `vcftools_pi.sh` script.

## Output

The output is a file named `pi.txt` containing pi estimates for each window.

## Script details

The `pi estimation snakefile` script performs the following steps:

1. Calculates pi estimates for each window using vcftools.
2. Writes the pi estimates to a file named `.windowed.pi`.
