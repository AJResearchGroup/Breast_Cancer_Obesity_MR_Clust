library(dplyr)
library(broom)
library(tidyr)


### protein PGS lm ####

# load results of lm on protein and PGS

load("/proj/sens2017538/nobackup/marina/results/pgs_lm_V2.Rda")

# keep proteins only, (no filter for significant ones yet), label as outcome data

proteins <- filter(results_df, grepl("c1|c2|c3|c0", term))

proteins <- rename(proteins,
                   beta_outc = estimate,
                   outc_SE = StdError,
                   outc_CI_lower = CI_Lower,
                   outc_CI_upper = CI_Upper,
                   outc_p = p_value)


### BMI PGS lm (but without women in the pgs_lm)####


# load phenotype data
load("/proj/sens2017538/nobackup/marina/data/bmi_cov_phenos.Rda")


# remove individuals used in lm protein PGS to avoid overlap, i.e. females from protein analysis
load("/proj/sens2017538/nobackup/marina/data/bigset_fem_forprotPGSlm_V2.Rda")
bmi_cov_phenos <- anti_join(bmi_cov_phenos, bigset_fem, join_by(f.eid == f.eid))


# create column with rank-based inverse normal transformation
bmi_cov_phenos <- rename(bmi_cov_phenos, BMI.raw = BMI)

rntransf <- function(x) {
  tmp.1 <- (order(order(x[!is.na(x)]))-0.5)/length(x[!is.na(x)])
  tmp.2 <- rep(NA,length(x))
  tmp.2[!is.na(x)] <- tmp.1
  return(qnorm(tmp.2,0,1))
}

BMI <- rntransf(bmi_cov_phenos$BMI.raw)

bmi_cov_phenos$BMI <- BMI

# prepare covariates for BMI linear model
bmi_covariates <- c("age", "age_sq", "sexSR", "centre", "gen_array2", "PC1", "PC2", "PC3",
                    "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11",
                    "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18",
                    "PC19", "PC20")


# load PGS for each cluster
clust1PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust1FgGt_TS1.profile", header = TRUE)
clust2PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust2FgGt_TS1.profile", header = TRUE)
clust3PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust3FgGt_TS1.profile", header = TRUE)
clust0PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmipgs_TS1/clust0FgGt_TS1.profile", header = TRUE)

# join sets
clust1prs_pheno <- inner_join(bmi_cov_phenos, clust1PRS, by = c("f.eid" = "IID"))
clust2prs_pheno <- inner_join(bmi_cov_phenos, clust2PRS, by = c("f.eid" = "IID"))
clust3prs_pheno <- inner_join(bmi_cov_phenos, clust3PRS, by = c("f.eid" = "IID"))
clust0prs_pheno <- inner_join(bmi_cov_phenos, clust0PRS, by = c("f.eid" = "IID"))

# run lm and save summaries with CI
clust1_lm <- lm(formula = as.formula(paste0("BMI ~ SCORESUM +", paste0(bmi_covariates, collapse = "+"))),
                data = clust1prs_pheno)
clust1_lm <- tidy(clust1_lm, conf.int = TRUE)

clust2_lm <- lm(formula = as.formula(paste0("BMI ~ SCORESUM +", paste0(bmi_covariates, collapse = "+"))),
                data = clust2prs_pheno)
clust2_lm <- tidy(clust2_lm, conf.int = TRUE)

clust3_lm <- lm(formula = as.formula(paste0("BMI ~ SCORESUM +", paste0(bmi_covariates, collapse = "+"))),
                data = clust3prs_pheno)
clust3_lm <- tidy(clust3_lm, conf.int = TRUE)

clust0_lm <- lm(formula = as.formula(paste0("BMI ~ SCORESUM +", paste0(bmi_covariates, collapse = "+"))),
                data = clust0prs_pheno)
clust0_lm <- tidy(clust0_lm, conf.int = TRUE)

# join lm tibbles, remove intercept, change cluster "4" identifier to 0
allclusts_lm <- bind_rows(clust1_lm, clust2_lm, clust3_lm, clust0_lm,
                          .id = "cluster")
allclusts_lm <- mutate(allclusts_lm, cluster = ifelse(cluster == "4", "0", cluster))
allclusts_lm <- filter(allclusts_lm, term=='SCORESUM')

# rename "estimate" to "beta_exp", keep only identifier + beta_exp, make integer
betaexp <- select(allclusts_lm, cluster, estimate)
betaexp <- rename(betaexp, "beta_exp"="estimate")
betaexp$cluster <- as.integer(betaexp$cluster)

### join beta_exp to the proteins ####
proteins <- inner_join(proteins, betaexp)

### calculate MR beta, SE and CI ####

proteins <- mutate(proteins, MR_beta = beta_outc / beta_exp)
proteins <- mutate(proteins, MR_se = outc_SE / beta_exp)
proteins <- mutate(proteins, MR_CIlow = MR_beta - (1.96 * MR_se))
proteins <- mutate(proteins, MR_CIhigh = MR_beta + (1.96 * MR_se))

### join protein names to file ####

# remove _p from protein letter code (careful with dashes and underscores)
named_proteins_withcalcMR <- proteins
named_proteins_withcalcMR$protein <- sub("_p$", "", proteins$protein)

# save 
save(named_proteins_withcalcMR, file = "/proj/sens2017538/nobackup/marina/data/named_prots_calcMR_V3.Rda")



