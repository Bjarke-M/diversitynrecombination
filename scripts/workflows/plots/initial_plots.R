library(tidyverse)
library(gridExtra)
data_file <- snakemake@input[[1]]
out <- snakemake@output[[1]]

data <- read_delim(data_file, delim = '\t') %>%
  na.omit() %>% group_by(species) %>% mutate(recomb_bin = ntile(cm_per_mb,100))


plot1 <- data %>% filter(!(chr=='chrX')) %>% 
  filter(freq_mean > 0.7) %>% 
  group_by(species,recomb_bin) %>% 
  summarise(
    mean_pi = mean(PI),
    freq_mean = mean(freq_mean),
    mean_cm_per_mb = mean(cm_per_mb)
  ) %>%
  ggplot(aes(x = mean_cm_per_mb, y = mean_pi*freq_mean, col = species))+
  geom_point()+
  geom_smooth()+
  facet_wrap(species~.)+
  theme_minimal()+
  theme(legend.position = 'none')



plot2 <- data %>% filter(!(chr=='chrX')) %>% 
  filter(freq_mean > 0.7) %>% 
  group_by(species,recomb_bin) %>% 
  summarise(
    mean_pi = mean(PI),
    freq_mean = mean(freq_mean),
    mean_cm_per_mb = mean(cm_per_mb)
  ) %>%
  ggplot(aes(x = log10(mean_cm_per_mb), y = mean_pi*freq_mean, col = species))+
  geom_point()+
  geom_smooth()+
  facet_wrap(species~.)+
  theme_minimal()+
  theme(legend.position = 'none')

plot3 <- data %>% filter(!(chr=='chrX')) %>% 
  filter(freq_mean > 0.7) %>% 
  group_by(species,recomb_bin) %>% 
  summarise(
    mean_pi = mean(PI),
    freq_mean = mean(freq_mean),
    mean_cm_per_mb = mean(cm_per_mb)
  ) %>%
  ggplot(aes(x = recomb_bin, y = mean_pi*freq_mean, col= species))+
  geom_point()+
  geom_smooth()+
  facet_wrap(species~.)+
  theme_minimal()+
  theme(legend.position = 'none')

combined <- grid.arrange(plot1,plot2,plot3, ncol =1)

ggsave(out,combined, width = 10, height = 10)
