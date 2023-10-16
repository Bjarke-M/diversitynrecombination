import sgkit as sg
import numpy as np


# Load genotype data
genotypes = sg.load_dataset("1000G_phase3_chr22")

# Define window size and step size
window_size = 100000
step_size = 50000

# Calculate window-based heterozygosity
heterozygosity = sg.stats.window_heterozygosity(
    genotypes=genotypes,
    size=window_size,
    step=step_size,
    backend="numpy"
)

# Print the results
print(heterozygosity)
