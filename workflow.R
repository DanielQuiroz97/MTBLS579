# Add needed libraries
library(tidyverse)
library(RGCxGC)
library(colorRamps)

#### Importing Raw Chromatograms ####

## Import metadata ##
metadata <- read_csv('metadata.csv')
##
GB_01 <- read_chrom(name = 'data/01_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_02 <- read_chrom(name = 'data/02_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_04 <- read_chrom(name = 'data/04_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_05 <- read_chrom(name = 'data/05_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_06 <- read_chrom(name = 'data/06_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_07 <- read_chrom(name = 'data/07_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_08 <- read_chrom(name = 'data/08_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_09 <- read_chrom(name = 'data/09_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_10 <- read_chrom(name = 'data/10_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_11 <- read_chrom(name = 'data/11_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_12 <- read_chrom(name = 'data/12_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_13 <- read_chrom(name = 'data/13_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_14 <- read_chrom(name = 'data/14_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_16 <- read_chrom(name = 'data/16_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_17 <- read_chrom(name = 'data/17_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_18 <- read_chrom(name = 'data/18_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_19 <- read_chrom(name = 'data/19_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_21 <- read_chrom(name = 'data/21_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_22 <- read_chrom(name = 'data/22_GB.cdf', 5L, 
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_23 <- read_chrom(name = 'data/23_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_24 <- read_chrom(name = 'data/24_GB.cdf', 5L, 
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_25 <- read_chrom(name = 'data/25_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_26 <- read_chrom(name = 'data/26_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_27 <- read_chrom(name = 'data/27_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_28 <- read_chrom(name = 'data/28_GB.cdf', 5L, 
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_29 <- read_chrom(name = 'data/29_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_30 <- read_chrom(name = 'data/30_GB.cdf', 5L, 
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_31 <- read_chrom(name = 'data/31_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_32 <- read_chrom(name = 'data/32_GB.cdf', 5L,
                    x_cut = c(8, 49), y_cut = c(1.3, 5))
GB_34 <- read_chrom(name = 'data/34_GB.cdf', 5L, 
                    x_cut = c(8, 49), y_cut = c(1.3, 5))

# Plot all raw chromatograms
lapply(ls(), function(x){
  jpeg(paste0('plots/raw/', x, '.jpg'), width = 8, height = 7, units = 'in',
       res = 350)
  plot(eval(as.name(x)),  nlevels = 150, color.palette = matlab.like)
  dev.off()
})


#### Preprocessing ####
## Smoothing ##
GB_01_smth <- wsmooth(GB_01, penalty = 2, lambda = 1e1)
GB_02_smth <- wsmooth(GB_02, penalty = 2, lambda = 1e1)
GB_04_smth <- wsmooth(GB_04, penalty = 2, lambda = 1e1)
GB_05_smth <- wsmooth(GB_05, penalty = 2, lambda = 1e1)
GB_06_smth <- wsmooth(GB_06, penalty = 2, lambda = 1e1)
GB_07_smth <- wsmooth(GB_07, penalty = 2, lambda = 1e1)
GB_08_smth <- wsmooth(GB_08, penalty = 2, lambda = 1e1)
GB_09_smth <- wsmooth(GB_09, penalty = 2, lambda = 1e1)
GB_10_smth <- wsmooth(GB_10, penalty = 2, lambda = 1e1)
GB_11_smth <- wsmooth(GB_11, penalty = 2, lambda = 1e1)
GB_12_smth <- wsmooth(GB_12, penalty = 2, lambda = 1e1)
GB_13_smth <- wsmooth(GB_13, penalty = 2, lambda = 1e1)
GB_14_smth <- wsmooth(GB_14, penalty = 2, lambda = 1e1)
GB_16_smth <- wsmooth(GB_16, penalty = 2, lambda = 1e1)
GB_17_smth <- wsmooth(GB_17, penalty = 2, lambda = 1e1)
GB_18_smth <- wsmooth(GB_18, penalty = 2, lambda = 1e1)
GB_19_smth <- wsmooth(GB_19, penalty = 2, lambda = 1e1)
GB_21_smth <- wsmooth(GB_21, penalty = 2, lambda = 1e1)
GB_22_smth <- wsmooth(GB_22, penalty = 2, lambda = 1e1)
GB_23_smth <- wsmooth(GB_23, penalty = 2, lambda = 1e1)
GB_24_smth <- wsmooth(GB_24, penalty = 2, lambda = 1e1)
GB_25_smth <- wsmooth(GB_25, penalty = 2, lambda = 1e1)
GB_26_smth <- wsmooth(GB_26, penalty = 2, lambda = 1e1)
GB_27_smth <- wsmooth(GB_27, penalty = 2, lambda = 1e1)
GB_28_smth <- wsmooth(GB_28, penalty = 2, lambda = 1e1)
GB_29_smth <- wsmooth(GB_29, penalty = 2, lambda = 1e1)
GB_30_smth <- wsmooth(GB_30, penalty = 2, lambda = 1e1)
GB_31_smth <- wsmooth(GB_31, penalty = 2, lambda = 1e1)
GB_32_smth <- wsmooth(GB_32, penalty = 2, lambda = 1e1)
GB_34_smth <- wsmooth(GB_34, penalty = 2, lambda = 1e1)

# Plot all smoothing chromatograms
smothing_chrom <- ls()[grepl(pattern = 'smth',ls())]
lapply(smothing_chrom, function(x){
    jpeg(paste0('plots/smoothing/', x, '.jpg'), width = 8, height = 7, units = 'in',
       res = 350)
  plot(eval(as.name(x)),  nlevels = 150, color.palette = matlab.like)
  dev.off()
})

## Baseline Correction ##
# labmda = 1e3
GB_01_bsl <- baseline_corr(GB_01_smth, lambda = 1e3)
GB_02_bsl <- baseline_corr(GB_02_smth, lambda = 1e3)
GB_04_bsl <- baseline_corr(GB_04_smth, lambda = 1e3)
GB_05_bsl <- baseline_corr(GB_05_smth, lambda = 1e3)
GB_06_bsl <- baseline_corr(GB_06_smth, lambda = 1e3)
GB_07_bsl <- baseline_corr(GB_07_smth, lambda = 1e3)
GB_08_bsl <- baseline_corr(GB_08_smth, lambda = 1e3)
GB_09_bsl <- baseline_corr(GB_09_smth, lambda = 1e3)
GB_10_bsl <- baseline_corr(GB_10_smth, lambda = 1e3)
GB_11_bsl <- baseline_corr(GB_11_smth, lambda = 1e3)
GB_12_bsl <- baseline_corr(GB_12_smth, lambda = 1e3)
GB_13_bsl <- baseline_corr(GB_13_smth, lambda = 1e3)
GB_14_bsl <- baseline_corr(GB_14_smth, lambda = 1e3)
GB_16_bsl <- baseline_corr(GB_16_smth, lambda = 1e3)
GB_17_bsl <- baseline_corr(GB_17_smth, lambda = 1e3)
GB_18_bsl <- baseline_corr(GB_18_smth, lambda = 1e3)
GB_19_bsl <- baseline_corr(GB_19_smth, lambda = 1e3)
GB_21_bsl <- baseline_corr(GB_21_smth, lambda = 1e3)
GB_22_bsl <- baseline_corr(GB_22_smth, lambda = 1e3)
GB_23_bsl <- baseline_corr(GB_23_smth, lambda = 1e3)
GB_24_bsl <- baseline_corr(GB_24_smth, lambda = 1e3)
GB_25_bsl <- baseline_corr(GB_25_smth, lambda = 1e3)
GB_26_bsl <- baseline_corr(GB_26_smth, lambda = 1e3)
GB_27_bsl <- baseline_corr(GB_27_smth, lambda = 1e3)
GB_28_bsl <- baseline_corr(GB_28_smth, lambda = 1e3)
GB_29_bsl <- baseline_corr(GB_29_smth, lambda = 1e3)
GB_30_bsl <- baseline_corr(GB_30_smth, lambda = 1e3)
GB_31_bsl <- baseline_corr(GB_31_smth, lambda = 1e3)
GB_32_bsl <- baseline_corr(GB_32_smth, lambda = 1e3)
GB_34_bsl <- baseline_corr(GB_34_smth, lambda = 1e3)

# Plot all smoothing chromatograms
bsl_chrom <- ls()[grepl(pattern = 'bsl$',ls())]
lapply(bsl_chrom, function(x){
  jpeg(paste0('plots/baseline/', x, '.jpg'), width = 8, height = 7, units = 'in',
       res = 350)
  plot(eval(as.name(x)),  nlevels = 150, color.palette = matlab.like)
  dev.off()
})


## Peak Alignment ##
# Non-carriage control
conC_con <- list(GB_19 = GB_19_bsl, GB_21 = GB_21_bsl, GB_22 = GB_22_bsl,
            GB_23 = GB_23_bsl, GB_24 = GB_24_bsl, GB_25 = GB_25_bsl,
            GB_26 = GB_26_bsl, GB_27 = GB_27_bsl, GB_28 = GB_28_bsl,
            GB_29 = GB_29_bsl, GB_30 = GB_30_bsl, GB_31 = GB_31_bsl,
            GB_32 = GB_32_bsl, GB_34 = GB_34_bsl)
ncar_con <- batch_2DCOW(GB_18_bsl, conC_con, c(20, 40), c(2, 8))
names(ncar_con@Batch_2DCOW)[1] <- 'GB_18'
# S. Paratyphi carriage
spC_chrom <- list(GB_05 = GB_05_bsl, GB_06 = GB_06_bsl, GB_07 = GB_07_bsl)
SPc_align <- batch_2DCOW(GB_04_bsl, spC_chrom,  c(20, 40), c(2, 8))
names(SPc_align@Batch_2DCOW)[1] <- 'GB_04'
# S. Typhi carriage
stC_chrom <- list(GB_02 = GB_02_bsl, GB_08 = GB_08_bsl, GB_09 = GB_09_bsl,
                  GB_10 = GB_10_bsl, GB_11 = GB_11_bsl, GB_12 = GB_12_bsl,
                  GB_13 = GB_13_bsl, GB_14 = GB_14_bsl,
                  GB_16 = GB_16_bsl, GB_17 = GB_17_bsl)
stc_align <- batch_2DCOW(GB_01_bsl, stC_chrom,  c(20, 40), c(2, 8))
names(stc_align@Batch_2DCOW)[1] <- 'GB_01'

all_chrom <- join_chromatograms(ncar_con, SPc_align, stc_align,
                                groups = metadata)

M579_pca <- m_prcomp(all_chrom)
save(M579_pca, file = 'MPCA.RData')