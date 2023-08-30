import os
import pandas as pd

# Define the path where the species folders are located
species_path = "primatediversity/data/gVCF_jointCalling_17_05_2021"

# Get a list of species folders within the specified path
species_folders = [species for species in os.listdir(species_path) if os.path.isdir(os.path.join(species_path, species))]

# Create a DataFrame to store the species names
species_df = pd.DataFrame({"Species Name": species_folders})

# Save the DataFrame as a CSV file
output_csv_path = "species_names.csv"
species_df.to_csv(output_csv_path, index=False)

print("Species names saved to:", output_csv_path)
