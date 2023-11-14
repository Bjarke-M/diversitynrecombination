library(tidyverse)
generate_df <- function(inputfile, window_size, sex, outfile) {
  # Read the coverage data from the input file
  coverage <- read.table(inputfile, col.names = c('chrA', 'startA', 'stopA', 'chrB', 'startB', 'stopB', 'n')) %>%
    select(chrA, startA, stopA, n) %>%
    group_by(chrA, startA, stopA) %>%
    summarize(sum_n = sum(n)) %>%
    mutate(
      window_size = stopA - startA,
      freq = ifelse(chrA == 'chrX' & sex == 'M', NA, sum_n / window_size) # make sure that we calculate the X chromosome correctly for the males in the data set.
    ) %>%
    select(chrA, startA, stopA, window_size, sum_n, freq)
  
  write.csv(coverage, file = outfile, row.names = FALSE)
}

generate_df(snakemake@input[[1]], snakemake@params[['window_size']], snakemake@params[['sex']], snakemake@output[[1]])
    