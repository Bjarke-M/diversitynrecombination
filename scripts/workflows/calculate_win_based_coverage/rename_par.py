import sys
def rename_par_region(input_file, output_file, par_end):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
                parts = line.strip().split('\t')  # Split the line into parts using a tab character as the delimiter
                chromosome, start, end = parts
                start = int(start)
                end = int(end)  
                par_end = int(par_end)
                if chromosome=='chrX' and int(end) <= int(par_end):
                    chromosome = 'PAR'
                outfile.write(f"{chromosome}\t{start}\t{end}\n")

rename_par_region(sys.argv[1], sys.argv[2], sys.argv[3])