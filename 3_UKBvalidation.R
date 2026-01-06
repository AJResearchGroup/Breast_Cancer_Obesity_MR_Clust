### new script for UKB validation - bmi and bc linear models
### but with more covariates than before:  age, age_sq, centre, sex, PCs, gen array 

library(dplyr)
library(broom)
library(writexl)


## load BMI + covariates data, keep only fields needed
load("/proj/sens2017538/nobackup/marina/data/bmi_cov_phenos.Rda") 

bmi_phenos <- select(bmi_cov_phenos, 
                     f.eid, BMI, sexGEN, age, age_sq, centre, gen_array2, PC1:PC20)

# create column with rank-based inverse normal transformation for BMI
bmi_phenos <- rename(bmi_phenos, BMI.raw = BMI)

rntransf <- function(x) {
  tmp.1 <- (order(order(x[!is.na(x)]))-0.5)/length(x[!is.na(x)])
  tmp.2 <- rep(NA,length(x))
  tmp.2[!is.na(x)] <- tmp.1
  return(qnorm(tmp.2,0,1))
}

BMI <- rntransf(bmi_phenos$BMI.raw)

bmi_phenos$BMI <- BMI

## load bc pheno data
load("/proj/sens2017538/nobackup/marina/data/bc_phenos2.Rda")

# join according to IID
bc_bmi_phenos <- inner_join(bc_phenos, bmi_phenos, by = c("IID" = "f.eid"))


# load PGS for each cluster
clust1PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust1FgGt_TS1.profile", header = TRUE)
clust2PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust2FgGt_TS1.profile", header = TRUE)
clust3PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust3FgGt_TS1.profile", header = TRUE)
clust0PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust0FgGt_TS1.profile", header = TRUE)

# join sets
clust1prs_pheno <- inner_join(bc_bmi_phenos, clust1PRS, by = c("IID" = "IID"))
clust2prs_pheno <- inner_join(bc_bmi_phenos, clust2PRS, by = c("IID" = "IID"))
clust3prs_pheno <- inner_join(bc_bmi_phenos, clust3PRS, by = c("IID" = "IID"))
clust0prs_pheno <- inner_join(bc_bmi_phenos, clust0PRS, by = c("IID" = "IID"))

#### linear model for BMI ####

# run lm and save summaries with CI 
clust1_lm <- lm(BMI ~ SCORESUM 
                + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                  PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                  PC19 + PC20, 
                data = clust1prs_pheno)
clust1_lm <- tidy(clust1_lm, conf.int = TRUE)

clust2_lm <- lm(BMI ~ SCORESUM 
                + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                  PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                  PC19 + PC20, 
                data = clust2prs_pheno)
clust2_lm <- tidy(clust2_lm, conf.int = TRUE)

clust3_lm <- lm(BMI ~ SCORESUM 
                + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                  PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                  PC19 + PC20, 
                data = clust3prs_pheno)
clust3_lm <- tidy(clust3_lm, conf.int = TRUE)

clust0_lm <- lm(BMI ~ SCORESUM 
                + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                  PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                  PC19 + PC20, 
                data = clust0prs_pheno)
clust0_lm <- tidy(clust0_lm, conf.int = TRUE)


# add identifier column to lm dataframes
clust1_lm$clust <- "clust1"
clust2_lm$clust <- "clust2"
clust3_lm$clust <- "clust3"
clust0_lm$clust <- "clust0"

# join lm tibbles, keep only SCORESUM (PGS)
allclusts_lm <- bind_rows(clust1_lm, clust2_lm, clust3_lm, clust0_lm,
                          .id = "cluster")
allclusts_lm <- filter(allclusts_lm, term=='SCORESUM')


#### log model for BC ####

# first filter to keep only women ! 
clust1prs_BCpheno <- filter(clust1prs_pheno, sexGEN == 0)
clust2prs_BCpheno <- filter(clust2prs_pheno, sexGEN == 0)
clust3prs_BCpheno <- filter(clust3prs_pheno, sexGEN == 0)
clust0prs_BCpheno <- filter(clust0prs_pheno, sexGEN == 0)


# run logm and save summaries with CI
clust1_lm_bc <- glm(BC ~ SCORESUM 
                    + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                      PC19 + PC20, 
                    data = clust1prs_BCpheno, family = binomial)
clust1_lm_bc <- tidy(clust1_lm_bc, conf.int = TRUE)


clust2_lm_bc <- glm(BC ~ SCORESUM 
                    + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                      PC19 + PC20, 
                    data = clust2prs_BCpheno, family = binomial)
clust2_lm_bc <- tidy(clust2_lm_bc, conf.int = TRUE)

clust3_lm_bc <- glm(BC ~ SCORESUM
                    + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                      PC19 + PC20, 
                    data = clust3prs_BCpheno, family = binomial)
clust3_lm_bc <- tidy(clust3_lm_bc, conf.int = TRUE)

clust0_lm_bc <- glm(BC ~ SCORESUM 
                    + age + age_sq + centre + gen_array2 + PC1 + PC2 + PC3 +
                      PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                      PC19 + PC20, 
                    data = clust0prs_BCpheno, family = binomial)
clust0_lm_bc <- tidy(clust0_lm_bc, conf.int = TRUE)

# add identifier column to log model dataframes
clust1_lm_bc$clust <- "clust1"
clust2_lm_bc$clust <- "clust2"
clust3_lm_bc$clust <- "clust3"
clust0_lm_bc$clust <- "clust0"

# join lm tibbles, keep only SCORESUM (PGS)
allclusts_lm_bc <- bind_rows(clust1_lm_bc, clust2_lm_bc, clust3_lm_bc, clust0_lm_bc,
                             .id = "cluster")
allclusts_lm_bc <- filter(allclusts_lm_bc, term=='SCORESUM')

#### join BMI and BC lin and log models ####

# add identifiers
allclusts_lm$model <- "lm"
allclusts_lm_bc$model <- "glm"

# bind rows
all_models <- bind_rows(allclusts_lm, allclusts_lm_bc, .id = "model_n")

# keep only columns wanted
all_models <- select(all_models, 
                     cluster:model)

# export 
write_xlsx(all_models, "/proj/sens2017538/nobackup/marina/results/UKBvalidation_allmodels.xlsx")

## calculated CIs and Wald ratio MR estimates in excel sheet (Cluster CI calculations)
