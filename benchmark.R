library(rbenchmark)
library(colorRamps)
library(tidyverse)
library(magrittr)



GB_01 <- read_chrom('data/01_GB.cdf', 5L)
GB_02 <- read_chrom('data/02_GB.cdf', 5L)
# Load bechmark.RData
mpca <- m_prcomp(all_chrom)
bench <- lapply(c(1, seq(10, 100, by = 10)), function(x, chrom1,
                                                      chrom2,mpca_data,
                                                      mpca){
  bench_data <- benchmark('read' = read_chrom('data/01_GB.cdf', 5L), 
                          'plot_chrom' = plot(chrom1,
                                              color.palette = colorRamps::matlab.like),
                          'smooth1' = wsmooth(chrom1, penalty = 1, lambda = 1e1),
                          'smooth2' = wsmooth(chrom1, penalty = 2, lambda = 1e1),
                          'baseline' =  baseline_corr(chrom1, lambda = 1e3),
                          '2DCOW' = twod_cow(chrom1, chrom2, c(20, 40), c(2, 8)),
                          'MPCA' = m_prcomp(mpca_data),
                          'plot_loading' = plot_loading(mpca,
                                                        color.palette = colorRamps::matlab.like),
                          replications = x)
  bench_data %<>% mutate(Repetitions = factor(x))
  bench_data
}, chrom1 = GB_01, chrom2 = GB_02, mpca_data = all_chrom, mpca = mpca)

data <- bench %>% bind_rows() %>% 
  mutate(Repetitions = as.numeric(Repetitions)) 

save(data, file = 'bech_results.RData')
