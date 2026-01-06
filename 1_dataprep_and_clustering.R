# ResAs run 2 after troubleshoot 1 #

library(TwoSampleMR)
library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(mrclust)
library(writexl)

setwd("/Users/marma900/Documents/RStudio/ResAs Project")


#### BMI data prep: meta analysis Mutie et al. 2025 combining FinnGen r9 and GIANT ####

# load Pascal's BMI data
BMI_GWAS <- readRDS("GT_FN_lifted_META_Clumped.rds")

# load Pascal's lifted (converted to GRCh37) r9 FinnGen data to get chromosome and position to add to meta dataset
FG_37 <- readRDS("finngen_lifted_toGRCh37.rds")
FG_37 <- select(FG_37, chr, pos, rsids)
BMI_ftd <- inner_join(BMI_GWAS, FG_37, by = c("SNP" = "rsids"))

BMI_ftd <-  unite(BMI_ftd, chr.pos.BMI, chr, pos, sep = ".", remove = FALSE)

# check SNPs p<5e-8, remove non-biallelic from effect and other allele, and duplicates
BMI_ftd <- filter(BMI_ftd, pval.exposure < 5e-8)

BMI_ftd <- BMI_ftd[nchar(BMI_ftd$effect_allele.exposure) == 1, ]

BMI_ftd <- BMI_ftd[nchar(BMI_ftd$other_allele.exposure) == 1, ]

BMI_ftd <- BMI_ftd[!duplicated(BMI_ftd$SNP), ]

# add sample size as reported from FinnGen+GIANT to run Steiger filtering
BMI_ftd <- mutate(BMI_ftd, sample_size = 588284)

# format to exposure dataset
bmi_exp <- format_data(BMI_ftd, type = "exposure", 
                       snp_col = "SNP", 
                       beta_col = "beta.exposure", 
                       se_col = "se.exposure", 
                       # eaf_col = "eaf.exposure",
                       samplesize_col = "sample_size",
                       
                       effect_allele_col = "effect_allele.exposure", 
                       other_allele_col = "other_allele.exposure",
                       pval_col = "pval.exposure", 
                       chr_col = "chr", 
                       pos_col = "pos" )


#### BC data prep ####
# load BC mixed GWAS data, keep meta elements only

BC_GWAS_mix <- vroom::vroom("/Users/marma900/Documents/RStudio/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")

BC_meta <- select(BC_GWAS_mix, var_name, Freq.Gwas, Effect.Meta, Baseline.Meta, Beta.meta, sdE.meta, p.meta)

BC_meta <- separate(BC_meta, var_name, into = c("chr", "pos", "NonEffAll", "EffAll"), sep = "_")

## both
# create mini BMI set with only rsIDs and chr.pos.BMI to align to BC set
# make chr.pos.BC column in BC set 
# join both by chr.pos columns -> will only keep SNPs of interest and have rsIDs

BMI_exp_mini <-select(BMI_ftd, chr.pos.BMI, SNP)

BC_meta <-  unite(BC_meta, chr.pos.BC, chr, pos, sep = ".", remove = FALSE)

BC_meta_ftd <- BC_meta %>% inner_join(BMI_exp_mini, by = c("chr.pos.BC" = "chr.pos.BMI"))

# add columns to outcome (BC) as found in BCAC publication: N.cases, N.controls, N.total

BC_meta_ftd <- mutate(BC_meta_ftd, 
                      N.cases = 133384,
                      N.controls = 113789,
                      N.total = 247173)

## prepare BC as outcome data 

bc_outc <- format_data(BC_meta_ftd, type = "outcome", 
                       snp_col = "SNP", 
                       beta_col = "Beta.meta", 
                       se_col = "sdE.meta", 
                       eaf_col = "Freq.Gwas",
                       effect_allele_col = "Effect.Meta", 
                       other_allele_col = "Baseline.Meta",
                       
                       ncase_col = "N.cases",
                       ncontrol_col = "N.controls",
                       samplesize_col = "N.total",
                       
                       
                       pval_col = "p.meta", 
                       chr_col = "chr", 
                       pos_col = "pos")


#### harmonisation - joining both ####

# harmonisation, and remove any palindromic SNPs
dat <- harmonise_data(
  exposure_dat = bmi_exp, 
  outcome_dat = bc_outc)

dat <- dat[dat$mr_keep == TRUE, ]



#### TSMR Steiger filtering ####
# remove all SNPs with BOTH significant Steiger p-value 
# AND stronger association of SNP with outcome than with exposure (i.e. steiger_dir = FALSE)
# i.e. remove all SNPs that are significantly more associated with the outcome than exposure
# or, keep all SNPs sign. more associated with exposure than outcome, 
# and SNPs who seem more associated with outcome than exposure, but not significantly

dat <- steiger_filtering(dat)

dat <- dat[!(dat$steiger_dir == FALSE & dat$steiger_pval < 0.05), ]

#### MR Clust ####

# load data if needed 
load("preClust_FN-GN-BCAC_libSt.Rda")

# take values from TwoSampleMR "dat" dataframe, name after vectors that will be needed for MRClust

bx <- dat$beta.exposure
by <- dat$beta.outcome
bxse <- dat$se.exposure
byse <- dat$se.outcome
snp_names <- dat$SNP

ratio_est <- by / bx
ratio_est_se <- byse / abs(bx)

# run MRClust with values needed for it
# important!!! set seed immediately before 

set.seed(5)
res_mrc_FgGnBCAC_libSt <- mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, 
                                      bx = bx, by = by, bxse = bxse, byse = byse, obs_names = snp_names)

# plot all results

plot_BMI_BC_best = res_mrc_FgGnBCAC_libSt$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with BMI") +
  ggplot2::ylab("Genetic association with BC") +
  ggplot2::ggtitle("All SNPs and their association with BMI and BC (FG, Giant, BCAC; lib St ft)")

plot_BMI_BC_best


#### save res and clust data ####

clustdata_FgGnBCAC_libSt <- res_mrc_FgGnBCAC_libSt$plots$two_stage$data
save(clustdata_FgGnBCAC_libSt, file = "clustdata_FgGnBCAC_run2.Rda")

save(res_mrc_FgGnBCAC_libSt, file = "clustresults_FgGnBCAC_run2.Rda")

#### TSMR per cluster ####

# if needed, load full data and cluster data

load("preClust_FN-GN-BCAC_libSt.Rda")
load("clustdata_FgGnBCAC_run2.Rda")

# create dataset for every cluster
clust1_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 1, ]
clust2_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 2, ]
clust3_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 3, ]
clust0_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 4, ]

# join clust and data prepared in TSMR for individual clusters, by SNP rsIDs

clust1_lib_forTSMR <- dat %>% inner_join(clust1_lib, by = c("SNP" = "observation"))
clust2_lib_forTSMR <- dat %>% inner_join(clust2_lib, by = c("SNP" = "observation"))
clust3_lib_forTSMR <- dat %>% inner_join(clust3_lib, by = c("SNP" = "observation"))
clust0_lib_forTSMR <- dat %>% inner_join(clust0_lib, by = c("SNP" = "observation"))

# run TSMR for each cluster

clust1_TSMR_res_lib <- mr(clust1_lib_forTSMR)
clust2_TSMR_res_lib <- mr(clust2_lib_forTSMR)
clust3_TSMR_res_lib <- mr(clust3_lib_forTSMR)
clust0_TSMR_res_lib <- mr(clust0_lib_forTSMR)

# join TSMR results
allclusts_tsmr <- bind_rows(clust1_TSMR_res_lib, clust2_TSMR_res_lib, 
                            clust3_TSMR_res_lib, clust0_TSMR_res_lib,
                          .id = "cluster")

allclusts_tsmr <- filter(allclusts_tsmr, method=='Inverse variance weighted')

#### prepare txt files for PGS calculations, using FinnGen+GIANT weights ####

load("clustdata_FgGnBCAC_run2.Rda")
BMI_GWAS <- readRDS("GT_FN_lifted_META_Clumped.rds")

# create dataset for every cluster
clust1_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 1, ]
clust2_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 2, ]
clust3_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 3, ]
clust0_lib<- clustdata_FgGnBCAC_libSt[clustdata_FgGnBCAC_libSt$cluster == 4, ]

# join clust and data prepared in TSMR for individual clusters, by SNP rsIDs

clust1_lib_forPGS <- BMI_GWAS %>% inner_join(clust1_lib, by = c("SNP" = "observation"))
clust2_lib_forPGS <- BMI_GWAS %>% inner_join(clust2_lib, by = c("SNP" = "observation"))
clust3_lib_forPGS <- BMI_GWAS %>% inner_join(clust3_lib, by = c("SNP" = "observation"))
clust0_lib_forPGS <- BMI_GWAS %>% inner_join(clust0_lib, by = c("SNP" = "observation"))

# rename columns
clust1_libSt_forPGS <- rename(clust1_libSt_forPGS, EA = effect_allele.exposure, BETA = beta.exposure)
clust2_libSt_forPGS <- rename(clust2_libSt_forPGS, EA = effect_allele.exposure, BETA = beta.exposure)
clust3_libSt_forPGS <- rename(clust3_libSt_forPGS, EA = effect_allele.exposure, BETA = beta.exposure)
clust0_libSt_forPGS <- rename(clust0_libSt_forPGS, EA = effect_allele.exposure, BETA = beta.exposure)

# keep only columns needed 
clust1_libSt_forPGS <- select(clust1_libSt_forPGS, SNP, EA, BETA, cluster)
clust2_libSt_forPGS <- select(clust2_libSt_forPGS, SNP, EA, BETA, cluster)
clust3_libSt_forPGS <- select(clust3_libSt_forPGS, SNP, EA, BETA, cluster)
clust0_libSt_forPGS <- select(clust0_libSt_forPGS, SNP, EA, BETA, cluster)

# export as txt files
write.table(clust1_libSt_forPGS, "clusters_PGS_TS1/clust1FgGt_TS1.txt", quote = FALSE, row.names = FALSE)
write.table(clust2_libSt_forPGS, "clusters_PGS_TS1/clust2FgGt_TS1.txt", quote = FALSE, row.names = FALSE)
write.table(clust3_libSt_forPGS, "clusters_PGS_TS1/clust3FgGt_TS1.txt", quote = FALSE, row.names = FALSE)
write.table(clust0_libSt_forPGS, "clusters_PGS_TS1/clust0FgGt_TS1.txt", quote = FALSE, row.names = FALSE)
