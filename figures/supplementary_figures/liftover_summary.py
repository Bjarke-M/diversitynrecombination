import os
import re

# Define the directory where the log files are stored
log_directory = "data/liftover/logs/"



ref_list = ['Daubentonia_madagascariensis','Saguinus_midas','Chlorocebus_aethiops','Rhinopithecus_roxellana',
                'Mandrillus_sphinx','Macaca_mulatta','Loris_tardigradus','Callithrix_jacchus',
                'Pithecia_pithecia','Lemur_catta-Thomas','Gorilla_gorilla_gorilla','Aotus_nancymaae','Sapajus_apella',
                'Pongo_abelii','Cercocebus_atys','Galago_moholi','Nomascus_leucogenys','Colobus_guereza',
                'Cebus_albifrons','Cercopithecus_mitis','Pongo_pygmaeus','Pan_troglodytes','Erythrocebus_patas',
                'Microcebus_murinus','Atele_fusciceps','Otolemur_garnettii','Nycticebus_pygmaeus']

# Function to extract total entries and failed to map from log file
def extract_counts_from_log(log_file):
    total_entries = 0
    failed_to_map = 0
    with open(log_file, 'r') as file:
        for line in file:
            total_match = re.search(r'Total entries: (\d+)', line)
            failed_match = re.search(r'Failed to map: (\d+)', line)
            if total_match:
                total_entries = int(total_match.group(1))
            elif failed_match:
                failed_to_map = int(failed_match.group(1))
    return total_entries, failed_to_map

# Function to collect counts for all species
def collect_counts_for_species():
    species_counts = {}
    for species_name in os.listdir(log_directory):
        species_dir = os.path.join(log_directory, species_name)
        if os.path.isdir(species_dir):
            log_file_path = os.path.join(species_dir, f"{species_name}_lifted.log")
            if os.path.exists(log_file_path):
                total, failed = extract_counts_from_log(log_file_path)
                species_counts[species_name] = (total, failed)
    return species_counts

# Function to print the table
def print_table(species_counts):
    print("Reference_genome\tTotal_Entries\tFailed_to_Map")
    for species, counts in species_counts.items():
        total, failed = counts
        print(f"{species}\t{total}\t{failed}")

# Main function
def main():
    species_counts = collect_counts_for_species()
    print_table(species_counts)

if __name__ == "__main__":
    main()
