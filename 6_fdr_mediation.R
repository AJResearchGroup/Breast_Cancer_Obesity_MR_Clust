### looking for mediating proteins, but using FDR correction instead of Bonferroni ##
### first on p values from the lm protein ~ PGSbmi + covariates 
### then on p values from Azimi et al.(2025) TSMR proteins-to-BC 

library(dplyr)
library(readxl)
library(ggplot2)
library(ggvenn)
library(tidyr)
library(broom)


### FDR part ####


## part I: BMI-protein, Wald ratio method MR estimates

# get data
load("/proj/sens2017538/nobackup/marina/data/named_prots_calcMR_V3.Rda")

# keep only columns needed and rename
BMIprot_MR <- select(named_proteins_withcalcMR, 
                     c(protein, cluster, term, MR_beta, MR_se, outc_p, 
                       MR_CIlow, MR_CIhigh))

BMIprot_MR <- rename(BMIprot_MR, 
                     beta_BMIprotMR = MR_beta, 
                     se_BMIprotMR = MR_se,
                     p_BMIprotMR = outc_p, 
                     CIlow_BMIprotMR = MR_CIlow,
                     CIhigh_BMIprotMR = MR_CIhigh)

# apply fdr correction on p values
BMIprot_MR <- mutate(BMIprot_MR, 
                      p_adj_BMIprotMR = p.adjust(BMIprot_MR$p_BMIprotMR, method = "fdr"))

# filter by adjusted p < 0.05 
BMIprot_MR_ftd <- filter(BMIprot_MR, p_adj_BMIprotMR <0.05)



## part II: protein to BC, TSMR

# get proteins-to-BC TSMR results and adjust first rows to get column names
protBC_TSMR <- read_excel("/proj/sens2017538/nobackup/marina/data/prot_breastcancer_TSMR_atleast_1cis.xlsx", 
                              col_names = FALSE)

protBC_TSMR <- protBC_TSMR[-1,]
colnames(protBC_TSMR) <- c("Protein","id_exposure","id_outcome","outcome",
                               "exposure","method","nsnp","b","se","pval", 
                               "OR","CI_lower", "CI_upper")
protBC_TSMR <- protBC_TSMR[-1,]

# keep IVW values only
protBC_TSMR <- filter(protBC_TSMR, method == "Inverse variance weighted")

# make columns numeric if needed
protBC_TSMR$b <- as.numeric(protBC_TSMR$b)
protBC_TSMR$pval <- as.numeric(protBC_TSMR$pval)
protBC_TSMR$se <- as.numeric(protBC_TSMR$se)
protBC_TSMR$OR <- as.numeric(protBC_TSMR$OR)
protBC_TSMR$CI_lower <- as.numeric(protBC_TSMR$CI_lower)
protBC_TSMR$CI_upper <- as.numeric(protBC_TSMR$CI_upper)
protBC_TSMR$nsnp <- as.numeric(protBC_TSMR$nsnp)

# keep only columns needed and rename
protBC_TSMR <- select(protBC_TSMR, c(Protein, b, se, pval))

protBC_TSMR <- rename(protBC_TSMR, 
                      protein = Protein, 
                      beta_protBC = b,
                      se_protBC = se,
                      p_protBC = pval)

# calculate CIs 
protBC_TSMR <- mutate(protBC_TSMR, 
                      CIlow_protBC = beta_protBC - (1.96 * se_protBC))
protBC_TSMR <- mutate(protBC_TSMR, 
                      CIhigh_protBC = beta_protBC + (1.96 * se_protBC))

# apply fdr correction to p_protBC
protBC_TSMR <- mutate(protBC_TSMR, 
                          p_adj_protBC = p.adjust(protBC_TSMR$p_protBC, method = "fdr"))

# filter by adjusted values <0.05
protBC_TSMR_ftd <- filter(protBC_TSMR, p_adj_protBC < 0.05)


## part III : join to keep only proteins where p<0.05 in both steps

sign_proteins <- inner_join(BMIprot_MR_ftd, protBC_TSMR_ftd)



#### calculate mediating effects ####

## join both mediation steps and select proteins identified as significant in at least one of the clusters
## to calculate mediating effects+CIs for all (non-null) clusters using Monte Carlo (Rmediation package)

sign_prots_allclusts <- inner_join(BMIprot_MR, protBC_TSMR)
sign_prots_allclusts <- filter(sign_prots_allclusts, cluster > 0)
sign_prots_allclusts <- filter(sign_prots_allclusts, grepl("CPM|CST6|MET\\b", protein))

sign_prots_allclusts <- mutate(sign_prots_allclusts, 
                               med_eff = beta_BMIprotMR * beta_protBC)

prod_mediation_rslts <- select(sign_prots_allclusts, 
                               protein, cluster, term, beta_BMIprotMR, se_BMIprotMR, 
                               beta_protBC, se_protBC)

## getting mediation estimates and CIs 

CI_MET1 <- RMediation::medci(mu.x = prod_mediation_rslts[1,4], 
                             mu.y = prod_mediation_rslts[1,6], 
                             se.x = prod_mediation_rslts[1,5], 
                             se.y = prod_mediation_rslts[1,7],
                             alpha = 0.05, 
                             type = "MC" ) 
CI_MET2 <- RMediation::medci(mu.x = prod_mediation_rslts[2,4], 
                             mu.y = prod_mediation_rslts[2,6], 
                             se.x = prod_mediation_rslts[2,5], 
                             se.y = prod_mediation_rslts[2,7],
                             alpha = 0.05, 
                             type = "MC" )
CI_MET3 <- RMediation::medci(mu.x = prod_mediation_rslts[3,4], 
                             mu.y = prod_mediation_rslts[3,6], 
                             se.x = prod_mediation_rslts[3,5], 
                             se.y = prod_mediation_rslts[3,7],
                             alpha = 0.05, 
                             type = "MC" )

CI_CST6_1 <- RMediation::medci(mu.x = prod_mediation_rslts[4,4], 
                             mu.y = prod_mediation_rslts[4,6], 
                             se.x = prod_mediation_rslts[4,5], 
                             se.y = prod_mediation_rslts[4,7],
                             alpha = 0.05, 
                             type = "MC" ) 
CI_CST6_2 <- RMediation::medci(mu.x = prod_mediation_rslts[5,4], 
                             mu.y = prod_mediation_rslts[5,6], 
                             se.x = prod_mediation_rslts[5,5], 
                             se.y = prod_mediation_rslts[5,7],
                             alpha = 0.05, 
                             type = "MC" )
CI_CST6_3 <- RMediation::medci(mu.x = prod_mediation_rslts[6,4], 
                             mu.y = prod_mediation_rslts[6,6], 
                             se.x = prod_mediation_rslts[6,5], 
                             se.y = prod_mediation_rslts[6,7],
                             alpha = 0.05, 
                             type = "MC" )

CI_CPM1 <- RMediation::medci(mu.x = prod_mediation_rslts[7,4], 
                             mu.y = prod_mediation_rslts[7,6], 
                             se.x = prod_mediation_rslts[7,5], 
                             se.y = prod_mediation_rslts[7,7],
                             alpha = 0.05, 
                             type = "MC" )
CI_CPM2 <- RMediation::medci(mu.x = prod_mediation_rslts[8,4], 
                             mu.y = prod_mediation_rslts[8,6], 
                             se.x = prod_mediation_rslts[8,5], 
                             se.y = prod_mediation_rslts[8,7],
                             alpha = 0.05, 
                             type = "MC" )
CI_CPM3 <- RMediation::medci(mu.x = prod_mediation_rslts[9,4], 
                             mu.y = prod_mediation_rslts[9,6], 
                             se.x = prod_mediation_rslts[9,5], 
                             se.y = prod_mediation_rslts[9,7],
                             alpha = 0.05, 
                             type = "MC" )

## format as df and add identifier column

MET1_df <- as.data.frame(CI_MET1)
MET1_df <- mutate(MET1_df, ID="MET_1")
MET2_df <- as.data.frame(CI_MET2)
MET2_df <- mutate(MET2_df, ID="MET_2")
MET3_df <- as.data.frame(CI_MET3)
MET3_df <- mutate(MET3_df, ID="MET_3")

CST6_1_df <- as.data.frame(CI_CST6_1)
CST6_1_df <- mutate(CST6_1_df, ID="CST6_1")
CST6_2_df <- as.data.frame(CI_CST6_2)
CST6_2_df <- mutate(CST6_2_df, ID="CST6_2")
CST6_3_df <- as.data.frame(CI_CST6_3)
CST6_3_df <- mutate(CST6_3_df, ID="CST6_3")

CPM1_df <- as.data.frame(CI_CPM1)
CPM1_df <- mutate(CPM1_df, ID="CPM_1")
CPM2_df <- as.data.frame(CI_CPM2)
CPM2_df <- mutate(CPM2_df, ID="CPM_2")
CPM3_df <- as.data.frame(CI_CPM3)
CPM3_df <- mutate(CPM3_df, ID="CPM_3")

## rotate for more succinct format

all_med_CIs <- bind_rows(MET1_df, MET2_df, MET3_df, 
                         CST6_1_df, CST6_2_df, CST6_3_df, 
                         CPM1_df, CPM2_df, CPM3_df)
all_med_CIs <- mutate(all_med_CIs, CI = c("low", "high", "low", "high", "low", "high", "low", "high", "low", "high", 
                                          "low", "high", "low", "high", "low", "high", "low", "high"))

wide_medCIs <- pivot_wider(all_med_CIs, 
                           names_from = CI, 
                           names_prefix = "CI_", 
                           values_from = X95..CI)
wide_medCIs <- separate(wide_medCIs, col = ID, into = c("Protein", "Cluster"), sep = "_")

## export 

library(writexl)
write_xlsx(wide_medCIs,"/proj/sens2017538/nobackup/marina/results/mediation_CIs.xlsx")


