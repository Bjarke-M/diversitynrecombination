from cyvcf2 import VCF

def count_singletons_in_windows(vcf_path, windows):
    """
    Count the number of singletons within specified genomic windows in a VCF file.

    Parameters:
        vcf_path (str): Path to the VCF file.
        windows (list of tuples): List of windows, each specified as (chromosome, start, end).
    
    Returns:
        dict: A dictionary where keys are windows (chr:start-end) and values are singleton counts.
    """
    vcf = VCF(vcf_path)
    singleton_counts = {f"{chrom}:{start}-{end}": 0 for chrom, start, end in windows}

    for variant in vcf:
        # Get variant position and chromosome.
        chrom = variant.CHROM
        pos = variant.POS
        
        # Get the genotype array for the variant.
        genotypes = variant.genotypes  # List of (GT, phased, [ploidy])
        
        # Count the number of heterozygous or homozygous alternate alleles.
        alt_counts = sum(1 for gt in genotypes if sum(gt[:2]) == 1)  # Heterozygous
        alt_counts += sum(1 for gt in genotypes if sum(gt[:2]) == 2)  # Homozygous alt
        
        # If exactly one individual has the alternate allele, it's a singleton.
        if alt_counts == 1:
            # Check if the variant falls within any of the specified windows.
            for chrom_window, start, end in windows:
                if chrom == chrom_window and start <= pos <= end:
                    key = f"{chrom}:{start}-{end}"
                    singleton_counts[key] += 1

    for window, count in singleton_counts.items():
        print(f"{window}: {count} singletons")
    return singleton_counts

# Example usage:
vcf_file = "your_file.vcf"
windows = [
    ("chr1", 100000, 200000),
    ("chr1", 300000, 400000),
    ("chr2", 500000, 600000)
]

count_singletons_in_windows(vcf_file, windows)
