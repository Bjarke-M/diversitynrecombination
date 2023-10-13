import sys
def process_file(input_file, output_file, window_size, par_end, chr_end):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')  # Split the line into parts using a tab character as the delimiter
            chromosome, start, end = parts # Convert relevant variables to integers
            start = int(start)
            end = int(end)  
            par_end = int(par_end)
            chr_end = int(chr_end)
            window_size = int(window_size)
            if chromosome == 'chrX' and start < par_end and end <= par_end:
                 # If the chromosome is 'chrX' and the start is less than par_end, and the end is less than or equal to par_end,
                # set the chromosome to 'PAR' and update the 'end' and 'next_start' values
                chromosome = 'PAR'
                end = int(end)
                next_start = int(par_end)

            elif chromosome == 'chrX' and start < par_end and end > par_end:
                # If the chromosome is 'chrX', the start is less than par_end, and the end is greater than par_end,
                # set the chromosome to 'PAR', set 'end' to par_end, and update 'next_start'
                chromosome = 'PAR'
                end = int(par_end)
                next_start = int(par_end)

            elif chromosome == 'chrX' and int(start) > par_end and int(start)+window_size <= chr_end:
                # If the chromosome is 'chrX', the start is greater than par_end, and (start + window_size) is less than or equal to chr_end,
                # set the chromosome to 'chrX', update 'start', 'end', and 'next_start'
                chromosome = 'chrX'
                start = next_start
                end = int(next_start)+window_size
                next_start = int(end)

            elif chromosome == 'chrX' and next_start+window_size >= chr_end:
                # If the chromosome is 'chrX', and (next_start + window_size) is greater than or equal to chr_end,
                # update 'start' and 'end' values
                start = str(int(next_start))
                end = chr_end

                
            # Write the updated values to the output file, separated by tabs
            outfile.write(f"{chromosome}\t{start}\t{end}\n")

# Call the function with your input file, output file, window size, and initial PAR end



process_file(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
