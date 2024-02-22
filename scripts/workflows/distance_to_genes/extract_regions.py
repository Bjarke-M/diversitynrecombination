# take a bed file and extract chr, and exon start and end positions for each mRNA coding gene
#output a bed file with chr, start, end, gene name 
#usage: python extract_regions.py <bed_file> <output_file>

import sys

# Open the input and output files
def extract_regions(bed_file, output_file):
    with open(bed_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
        regions = {}  # Store regions as a dictionary
        # Read the input file line by line
        for line in infile:
            # Split the line into an array
            array = line.split('\t')
            # Get the gene name
            gene = array[3].split('.')[0]
            # Check if the gene starts with NM_
            if gene.startswith('NM_'):
                # Create a key with the chromosome, start and end positions
                key = (array[0], array[1], array[2])
                # If this key is new, initialize it with an empty list
                if key not in regions:
                    regions[key] = []
                # Append the gene name to the list of this key
                regions[key].append(array[3])
                print(regions)
        # Write the regions to the output file
    #     for key, gene_names in regions.items():
    #         gene_list_str = ",".join(gene_names)
    #         outfile.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{gene_list_str}\n")
    # return outfile

# Get input and output file names from command line arguments
infile = sys.argv[1]
outfile = sys.argv[2]

# Call the function
extract_regions(infile, outfile)
