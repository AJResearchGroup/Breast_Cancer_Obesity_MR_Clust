
library(dplyr)
library(readxl)
library(ggplot2)
library(ggvenn)
library(tidyr)
library(broom)

#### Sobel mediation analysis, post-FDR ####

## (first, repetition of BMIprot_MR and protBC_TSMR preparation)

## BMIprot_MR data preparation ##

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

## protBC_TSMR data preparation ##

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


## joining and joint manipulations ## 
for_mediation <- inner_join(BMIprot_MR, protBC_TSMR)
for_mediation <- filter(for_mediation, cluster > 0)
for_mediation <- filter(for_mediation, grepl("CPM|CST6|MET\\b", protein))
for_mediation <- select(for_mediation, protein, cluster, term, 
                        beta_BMIprotMR, se_BMIprotMR, beta_protBC, se_protBC)


## calculate beta_diff and se_diff for all cluster pairs

# get protein names 
prot_list <- unique(for_mediation$protein) 

# initiate df for results
med_res <- data.frame()

# make "temp", df which is only for one protein
for (prot in prot_list) {
  temp <- for_mediation[which(for_mediation$protein == prot), ]
  
  # get all possible combinations of clusters (differences we want)
  combs <- combn(temp$cluster, 2)
  numb_comb <- ncol(combs)  # nr of columns of combs = the number of combinations
  
  # initialise "temp_diff", temporary df that will hold calculated differences
  temp_diff <- matrix(nrow = numb_comb, ncol = 5)
  
  # for all combination pairs
  for (i in 1:numb_comb){
    # take one pair ("i'th" column of combs)
    a_comb <- combs[,i]
    # assign beta from the row corresponding to the 1st term in the i'th pair to x1
    x1 <- temp$beta_BMIprotMR[a_comb[1]] 
    # assign beta from the row corresponding to the 2nd term in the i'th pair to x2
    x2 <- temp$beta_BMIprotMR[a_comb[2]]
    # calculate difference for that pair
    beta_diff = x1-x2
    
    # save info in temp_diff
    temp_diff[i, 1] <- prot
    temp_diff[i, 2] <- a_comb[1]
    temp_diff[i, 3] <- a_comb[2]
    temp_diff[i, 4] <- beta_diff
    
    # assign SEs for same pair to x3 and x4
    x3 <- temp$se_BMIprotMR[a_comb[1]]
    x4 <- temp$se_BMIprotMR[a_comb[2]]
    # calculate se_diff for that pair
    se_diff <- sqrt((x3)^2 + (x4)^2)
    
    # save se_diff in temp_diff
    temp_diff[i, 5] <- se_diff
    
  }
  # attach temp_diff to main results df
  med_res <- rbind(med_res, temp_diff)
}
# rename columns, make beta and se_diff numeric
colnames(med_res) <- c("protein", "clust.a", "clust.b", "beta_diff", "se_diff") 
med_res$beta_diff <- as.numeric(med_res$beta_diff)
med_res$se_diff <- as.numeric(med_res$se_diff)


## get beta_diff_mediation, se_mediation and p_mediation ## 

# add beta_PROT_bc and se_PROT_bc to med_res
med_res <- left_join(med_res, select(protBC_TSMR, protein, beta_protBC, se_protBC))

# calculate beta_diff_mediation
med_res <- mutate(med_res, 
                  beta_diff_mediation = beta_diff * beta_protBC)

# calculate se_mediation (se from the delta method in Sobel's test)
med_res <- mutate(med_res, 
                  se_mediation = sqrt(beta_diff^2 * se_protBC^2 + beta_protBC^2 * se_diff^2))

# calculate p_mediation
med_res <- mutate(med_res, 
                  p_mediation = 2*pnorm(-abs(beta_diff_mediation)/se_mediation))



### added Monte-Carlo (?) bootstrap part since p values not quite unequivocal ####

df_temp<-med_res[, 1:3]
df_temp<-cbind(df_temp,data.frame(L95=rep(NA,nrow(df_temp)),U95=rep(NA,nrow(df_temp))))
for (r in 1:nrow(med_res))
{
  set.seed(876)
  
  # assign values from specific row to values needed in calculation
  beta_diff <- med_res[r, 4]
  se_diff <- med_res[r, 5]
  beta_protBC <- med_res[r, 6]
  se_protBC <- med_res[r, 7]
  
  # perform calculation
  quantiles <- quantile(rnorm(100000, beta_diff, se_diff)
                        *rnorm(100000, beta_protBC, se_protBC), 
                        probs=c(0.025,0.975))
  
  # get quantile values
  df_temp$L95[r] <- unname(quantiles[1])
  df_temp$U95[r] <- unname(quantiles[2])
  
}

med_res <- cbind(med_res, df_temp[,c(4,5)])


## keep important results
final_mediation_results <- select(med_res, 
                                  protein, clust.a, clust.b, 
                                  beta_diff_mediation, se_mediation, p_mediation, 
                                  L95, U95)


#### export med_res ####
library(writexl)
write_xlsx(med_res, "/proj/sens2017538/nobackup/marina/results/mediation_results.xlsx")

