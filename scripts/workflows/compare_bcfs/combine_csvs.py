import os
import csv
import sys

def combine_csvs(input_dir, output_file):
    # Header for the combined CSV file
    header = ['species', 'unique_erik', 'unique_me', 'intersect', 'unique_erik_unmap', 'unique_me_unmap', 'intersect_unmap']

    # Initialize a list to store the data from each CSV file
    combined_data = [header]

    # Iterate over the CSV files in the directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".csv"):
            file_path = os.path.join(input_dir, filename)
            
            # Extract the species name from the filename (assuming a specific naming convention)
            species = filename.split('.')[0]
            
            # Read the CSV file and append its data to the combined_data list
            with open(file_path, 'r') as csv_file:
                csv_reader = csv.reader(csv_file)
                next(csv_reader)  # Skip the header row
                for row in csv_reader:
                    combined_data.append([species] + row)

    # Write the combined data to the output CSV file
    with open(output_file, 'w', newline='') as csv_output:
        csv_writer = csv.writer(csv_output)
        csv_writer.writerows(combined_data)

    print("Combined CSV file created:", output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_csvs.py <input_directory> <output_file>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_csv = sys.argv[2]

    combine_csvs(input_directory, output_csv)
