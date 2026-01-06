# for submission with bash

library(dplyr)
library(tidyr)

load("/proj/sens2017538/nobackup/marina/data/bigset_fem_forprotPGSlm_V2.Rda")

# loop building

# Initialise df to store results
results_df <- data.frame(
  protein = character(),
  cluster = integer(),
  term = character(),
  estimate = numeric(),
  StdError = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

protein_name <- colnames(bigset_fem)[37:2958]


covariates <- c("age", "age_sq", "centre", "gen_array2", "PC1", "PC2", "PC3",
                "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11",
                "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18",
                "PC19", "PC20", "inclusion_ppp", "Batch.x")


for (i in protein_name) {
  for (c in 0:3) {
    
    # Construct the cluster variable
    cluster_var <- paste0("c", c)
    
    # Regression model
    model <- lm(
      formula = as.formula(paste0("`", i, "` ~", cluster_var, "+", paste0(covariates, collapse = " + "))),
      data = bigset_fem)
    
    # Extract summary of the model
    model_summary <- summary(model)
    coeff <- coef(model_summary)
    conf_int <-confint(model, level = 0.95)
    
    
    # Save relevant results in a data frame
    for (term in rownames(coeff)) {
      results_df <- rbind(
        results_df,
        data.frame(
          protein = i,
          cluster = c,
          term = term,
          estimate = coeff[term, "Estimate"],
          StdError = coeff[term, "Std. Error"],
          CI_Lower = conf_int[term, 1],
          CI_Upper = conf_int[term, 2],
          p_value = coeff[term, "Pr(>|t|)"]
        )
      )
      
    }
  }
}

save(results_df, file = "/proj/sens2017538/nobackup/marina/results/pgs_lm_V2.Rda")
