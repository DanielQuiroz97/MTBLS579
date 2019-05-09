library(tidyverse)
library(magrittr)
metadata <- read_csv('metadata.csv')
raw_fl <- list.files(pattern = '.cdf')
extracted <- data.frame(raw_fl = raw_fl,
                        names_fl = tools::file_path_sans_ext(raw_fl))
metadata <- merge(metadata, extracted, by.x = "Names", by.y = 'names_fl')
rm(extracted, raw_fl)

