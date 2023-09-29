library(tidyverse)

distribution_plot <- function(input_file, species_name, window_size, pd_id, outfile) {
  # Read the coverage data from the input file
  coverage <- read_table(input_file, col_names = c('chrA', 'startA', 'stopA', 'chrB', 'startB', 'stopB', 'n'))  %>% 
  select(chrA, stopA, n) %>%
  group_by(chrA, stopA) %>%
  summarise(freq = sum(n)/window_size)
 
 median_freq_df <- coverage %>%
 group_by(chrA) %>%
 na.omit() %>%
 summarise(median_freq = round(median(freq),2))
  
  
  # Create the histogram plot
  hist_plot <- coverage %>% 
    ggplot(aes(x = freq)) +
    geom_histogram(bins = 30) +
    geom_text(data = median_freq_df, aes(x = 0.5, y = 0.5, label = paste("Median:", median_freq)), vjust = -1) +
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