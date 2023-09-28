library(tidyverse)

distribution_plot <- function(input_file, species_name, window_size, pd_id, outfile) {
  # Read the coverage data from the input file
  coverage <- read_table(input_file, col_names = c('chrA', 'startA', 'stopA', 'chrB', 'startB', 'stopB', 'n'))
  
  # Calculate the count of n/window_size == 1 per facet
  count_equal_to_1 <- coverage %>%
    group_by(chrA) %>%
    filter(n / window_size == 1) %>%
    summarise(count_equal_to_1 = n())
  
  # Create the histogram plot
  hist_plot <- coverage %>% 
    ggplot(aes(x = n/window_size)) +
    geom_histogram(bins = 100) +
    geom_text(data = count_equal_to_1, aes(x = 0.5, y = 0.5, label = paste("#100%:", count_equal_to_1)), vjust = -1) +
    facet_wrap(~ chrA)+
    labs(
      title = sprintf("%s_%s_%d_coverage",
                      species_name, pd_id, window_size),
      x = '% coverage',
      y = 'counts'
    )
  
  ggsave(filename = outfile, plot = hist_plot, width = 8, height = 4)
}

distribution_plot(snakemake@input[[1]], snakemake@params[['species']], snakemake@params[['window_size']], snakemake@params[['pd_id']], snakemake@output[[1]])

