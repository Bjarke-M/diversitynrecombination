library(tidyverse)

# Read in data
# Pan_troglodytes_1kb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_1000.windowed.pi", delim = "\t")
# Pan_troglodytes_5kb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_5000.windowed.pi", delim = "\t")
# Pan_troglodytes_10kb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_10000.windowed.pi", delim = "\t")
# Pan_troglodytes_50kb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_50000.windowed.pi", delim = "\t")
# Pan_troglodytes_100kb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_100000.windowed.pi", delim = "\t")
# Pan_troglodytes_1mb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_1000000.windowed.pi", delim = "\t")
# Pan_troglodytes_2mb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_2000000.windowed.pi", delim = "\t")
# Pan_troglodytes_10mb <- read_delim("../results/windowed_pi/Pan_troglodytes/troglodytes/nonpar/troglodytes_10000000.windowed.pi", delim = "\t")

# list_of_dfs <- list(Pan_troglodytes_1kb,Pan_troglodytes_5kb,Pan_troglodytes_10kb,Pan_troglodytes_50kb,Pan_troglodytes_100kb,Pan_troglodytes_1mb,Pan_troglodytes_2mb,Pan_troglodytes_10mb)
# list_of_paths <- list('Pan/Pan_troglodytes_1kb','Pan/Pan_troglodytes_5kb','Pan/Pan_troglodytes_10kb','Pan/Pan_troglodytes_50kb','Pan/Pan_troglodytes_100kb','Pan/Pan_troglodytes_1mb','Pan/Pan_troglodytes_2mb','Pan/Pan_troglodytes_10mb')


# Gorilla_gorilla_1kb <- read_delim("../results/windowed_pi/Gorilla_gorilla_gorilla/beringei/nonpar/beringei_1000.windowed.pi", delim = "\t")
# Gorilla_gorilla_10kb <- read_delim("../results/windowed_pi/Gorilla_gorilla_gorilla/beringei/nonpar/beringei_10000.windowed.pi", delim = "\t")
# Gorilla_gorilla_100kb <- read_delim("../results/windowed_pi/Gorilla_gorilla_gorilla/beringei/nonpar/beringei_100000.windowed.pi", delim = "\t")
# Gorilla_gorilla_1mb <- read_delim("../results/windowed_pi/Gorilla_gorilla_gorilla/beringei/nonpar/beringei_1000000.windowed.pi", delim = "\t")
# Gorilla_gorilla_2mb <- read_delim("../results/windowed_pi/Gorilla_gorilla_gorilla/beringei/nonpar/beringei_2000000.windowed.pi", delim = "\t")

# list_of_dfs <- list(Gorilla_gorilla_1kb,Gorilla_gorilla_10kb,Gorilla_gorilla_100kb,Gorilla_gorilla_1mb,Gorilla_gorilla_2mb)
# list_of_paths <- list('Gorilla/Gorilla_gorilla_1kb','Gorilla/Gorilla_gorilla_10kb','Gorilla/Gorilla_gorilla_100kb','Gorilla/Gorilla_gorilla_1mb','Gorilla/Gorilla_gorilla_2mb')

Macaca_mulatta_1kb <- read_delim("../results/windowed_pi/Macaca_mulatta/mulatta/nonpar/mulatta_1000.windowed.pi", delim = "\t")
Macaca_mulatta_10kb <- read_delim("../results/windowed_pi/Macaca_mulatta/mulatta/nonpar/mulatta_10000.windowed.pi", delim = "\t")
Macaca_mulatta_100kb <- read_delim("../results/windowed_pi/Macaca_mulatta/mulatta/nonpar/mulatta_100000.windowed.pi", delim = "\t")

plot_pi_total_length <- function(df,title){
    plot <- df %>% group_by(CHROM) %>%
       summarise(mean_pi = mean(PI),
                median_pi = median(PI),
                length = max(BIN_END)) %>%
        ggplot(aes(x = length, y = median_pi, col = CHROM)) +
        geom_point() +
        ggtitle(title) # Add title here
    return(plot)
}

for (i in seq_along(list_of_dfs)) {
    df <- list_of_dfs[[i]]
    path <- paste0(list_of_paths[[i]], ".pdf")
    title <- gsub("/", "_", list_of_paths[[i]])
    plot <- plot_pi_total_length(df, title)
    ggsave(path, plot)
}

