################################################################################
# Concomitant lipid-lowering and anticoagulation therapy versus                #
# anticoagulation monotherapy for venous thromboembolism prevention:           #
# a systematic review and Bayesian meta-analysis                               #
# Code for Biomarker levels                                                    #
#         Created: 25.07.2025                                                  #
#         Last Updated: 13.08.2025                                             #
#         Coded by: Harsini R S                                                #
################################################################################

#------------------------------------------------------------------------------#
# Set up
#------------------------------------------------------------------------------#

# Load necessary libraries
library(rstanarm)
library(brms)
library(dplyr)
library(ggplot2)
library(ggridges)
library(readxl)
library(readr)
library(tidyr)
library(tidybayes)
library(forplo)
library(glue)
library(stringr)
library(forcats)
library(moments)
library(extrafont)
library(ggrepel)

##Loading fonts for the first time
font_import()

##load fonts
loadfonts(device = "win", quiet = TRUE)

#set working directory
getwd()
setwd("Z:/meta_ana")

#check Rtools for functioning of brms and stan
Sys.which("make")
Sys.setenv(PATH = paste("C:/rtools44/usr/bin", Sys.getenv("PATH"), sep=";"))

#------------------------------------------------------------------------------#
# Data Input
#------------------------------------------------------------------------------#

#read data sheet
bio_data <- read_xlsx("outcome_lltAC.xlsx", sheet="biomarkers")

##cal log values
bio_data <- bio_data %>%
  mutate(
    # Calculate Mean
    log_mean.e = log(mean.e),
    log_mean.c = log(mean.c),
    # Calculate Standard dev
    log_sd.e = log(sd.e),
    log_sd.c = log(sd.c)
  )

#cal log-md and log-smd and SE
bio_data <- bio_data %>%
  mutate(
    # Calculate Mean Difference (MD)
    MD = log_mean.e - log_mean.c,
    # Calculate Standard Error for MD (MD_SE)
    MD_SE = sqrt((log_sd.e^2 / n.e) + (log_sd.c^2 / n.c)),
    # Calculate Standardized Mean Difference (SMD) using pooled SD
    SMD = (log_mean.e - log_mean.c) / sqrt(((n.e - 1) * log_sd.e^2 + (n.c - 1) * log_sd.c^2) / (n.e + n.c - 2)),
    # Calculate Standard Error for SMD (SMD_SE)
    SMD_SE = sqrt((1 / n.e) + (1 / n.c) + ((SMD^2) / (2 * (n.e + n.c)))),
    SMD_var = SMD_SE^2
  )

ddimer <-subset(bio_data, Outcome == "D-dimer")
CRP <-subset(bio_data, Outcome == "CRP")


#------------------------------------------------------------------------------#
# Bayesian Meta Analysis Model
#------------------------------------------------------------------------------#

# Fit the Bayesian model using brms
model <- brm(
  #by default, the se terms replace the sigma, resulting in sigma 0 and NA for Rhat and ESS
  formula = SMD | se(SMD_SE) ~  1+(1 | Study), #by default, the se terms replace the sigma, resulting in sigma 0 and NA for Rhat and ESS
  data = ddimer, # prepared dataset - ddimer or CRP
  family = gaussian(),  # Gaussian distribution for the continuous log-transformed IRR
  prior = c(
    prior(normal(0,1), class = "Intercept"),  # Prior for the intercept
    prior(cauchy(0, 0.5), class = "sd")# Prior for the random effect SD
  ),
  chains = 4,  # Number of chains for MCMC sampling
  iter = 2000,  # Number of iterations per chain
  warmup = 1000,  # Number of warm-up iterations (burn-in)
  control = list(adapt_delta = 0.995), # Control parameter for better convergence
  save_pars = save_pars(all = TRUE)
)

##model diagnostics
pairs(model)

# Check the model summary
summary(model)

#extract group level estimates
ranef(model)

##posterior distribution
post.samples <- as_draws_df(model)

#rename for clarity
post.samples <- post.samples %>%
  rename(
    smd = b_Intercept,             # pooled effect size
    tau = sd_Study__Intercept     # between-study heterogeneity (τ)
  )

##check the probability of our pooled effect being smaller than 0
smd.ecdf <- ecdf(post.samples$smd)
smd.ecdf(0)

# For calculating I² from the posterior for VTE data
post.samples$I2 <- (post.samples$tau^2) / 
  (post.samples$tau^2 + mean(ddimer$SMD_SE^2)) * 100


# Summarize I²
quantile(post.samples$I2, probs = c(0.025, 0.5, 0.975))  # 95% credible interval for I²

##Skewness
skew<- skewness(post.samples$smd)
# Suppose 'post.samples' is a data frame with columns like:
# r_Study[Chang.2017,Intercept], r_Study[Del.Giorno.2024,Intercept], etc.
post.samples <- as.data.frame(post.samples)
# Extract study columns (adjust the pattern if needed)
study_cols <- grep("^r_Study\\[.*,Intercept\\]", names(post.samples), value=TRUE)
# Calculate skewness for each study
skewness_results <- sapply(post.samples[, study_cols], skewness)
# Make it a nice table
skewness_table <- data.frame(
  Study = sub("^r_Study\\[(.*),Intercept\\]", "\\1", names(skewness_results)),
  Skewness = skewness_results
)
# Remove the dot between name and year in the 'Study' column and add space
skewness_table$Study <- gsub("\\.", " ", skewness_table$Study)
skewness_table <- rbind(
  skewness_table,
  data.frame(Study = "Pooled Effect", Skewness = skew)
)

#------------------------------------------------------------------------------#
# Summarise results
#------------------------------------------------------------------------------#

#calculate the actual effect size of each study by adding the pooled effect size b_Intercept to the estimated deviation of each study
study.draws <- spread_draws(model, r_Study[Study,], b_Intercept) %>% 
  mutate(b_Intercept = r_Study + b_Intercept)

#generate the distribution of the pooled effect
pooled.effect.draws <- spread_draws(model, b_Intercept) %>% 
  mutate(Study = "Pooled Effect")

#bind study.draws and pooled.effect.draws together in one data frame.
#(1) clean the study labels (i.e. replace dots with spaces), and (2) reorder the study factor levels by effect size (high to low).
#The result is the data we need for plotting, which we save as forest.data.
forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(Study = str_replace_all(Study, "[.]", " ")) %>% 
  mutate(Study = reorder(Study, b_Intercept))

#the forest plot should also display the effect size (SMD and credible interval) of each study.
# group it by Study, and then use the mean_qi function to calculate these values. We save the output as forest.data.summary.
forest.data.summary<- group_by(forest.data, Study) %>% 
  median_qi(b_Intercept)

##sorting to ensure pooled effect is last
forest.data.summary<- forest.data.summary %>%
  arrange(Study == "Pooled Effect")#comes last

##add columns
forest.data.summary <- forest.data.summary%>%
  left_join(select(ddimer, Study, AC, LLT, Design), by = "Study")##choose datatset major, ICH, GIB 
forest.data.summary<- forest.data.summary%>%
  left_join(select(skewness_table, Study, Skewness), by = "Study")

##complete pooled effect data
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "AC"] <- "ACs"
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "LLT"] <- "Statins"
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "Design"] <- " "


##convert to dataframe for forplo
##convert forest.data.summary to specific outcome major, ICH, GIB forest plot data summary
##eg forest.data.summary.ich<-as.data.frame(forest.data.summary)
forest.data.summary.ddimer<-as.data.frame(forest.data.summary)
#or
forest.data.summary.crp<-as.data.frame(forest.data.summary)

###plotting###
loadfonts(device = "win", quiet = TRUE) 

#------------------------------------------------------------------------------#
# Forest plots
#------------------------------------------------------------------------------#

##forest plot single outcome
#save image
png(file = "forest_crp.png", width = 4000, height = 2000, res = 300)
forplo2(exp(forest.data.summary.crp[,c(2:4)]), 
        normal.pdf = 1:nrow(forest.data.summary.crp), ci.edge=F, row.labels = forest.data.summary.crp$Study,
        normal.skew = forest.data.summary.crp$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "CRP levels",
        add.columns = forest.data.summary.crp[,8:10],
        add.colnames = c('AC','LLT', 'Design'),
        #favorlabs = "Combined therapy vs Monotherapy",
        diamond.pdf = nrow(forest.data.summary.crp), diamond.col = 'red',
        column.spacing = 5, margin.right = 15, font = 'Corbel', shade.every = 1, shade.col = 4, em = 'RR',
        extra.line = "P(RR < 1) = 0.99      tau = 0.16      I² = 47%",extra.line.xshift = -2.97,extra.line.yshift = 0.25
)
dev.off()

##combined outcomes plots
##dataframes list
forest.data.summary.ddimer
forest.data.summary.crp

#combine dfs
bio_combo<-rbind(
  cbind(forest.data.summary.ddimer, Outcome = "ddimer"),
  cbind(forest.data.summary.crp, Outcome = "crp")
)

#assign groups with ID
bio_combo <-cbind(bio_combo, group_id = c(1,1,1,1,1,1,2,2,2,2,2,2))

##assign groups names
bio_combo$group_name <- ifelse(
  bio_combo$group_id == 1, "D-dimer",
  ifelse(
    bio_combo$group_id == 2, "C-reactive Protein",
    NA  # (optional: handles unexpected values)
  )
)


#plot
png(file = "forest_biomarkers.png", width = 4000, height = 3000, res = 300)
forplo2(bio_combo[,c(2:4)], 
        normal.pdf = 1:nrow(bio_combo), ci.edge=F, row.labels = bio_combo$Study,
        normal.skew = bio_combo$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "Biomarker Levels",
        add.columns = bio_combo[,8:10],
        add.colnames = c('AC','LLT','Design'),
        groups = as.numeric(factor(bio_combo$group_id)),#group
        grouplabs = unique(bio_combo$group_name),#group labels
        group.space = 2,
        diamond.pdf = c(6,12), diamond.col = 'red',
        column.spacing = 5, margin.right = 15, font = 'Corbel', shade.every = 1, shade.col = 4, em = 'log SMD',
        extra.line = c("P(SMD<0) = 0.88      tau = 0.2      I² = 27%", "P(SMD<0) = 0.98      tau = 0.71      I² = 81%")
        extra.line.y = c(11, 2),
        extra.line.xshift = -0.1
)
dev.off()

#------------------------------------------------------------------------------#
# Metaregression Model for biomarkers
#------------------------------------------------------------------------------#
bio_metareg<- read_xlsx("outcome_lltAC.xlsx", sheet="biomarker_age")

##cal log values
bio_metareg <- bio_metareg %>%
  mutate(
    # Calculate Mean
    log_mean.e = log(mean.e),
    log_mean.c = log(mean.c),
    # Calculate Standard dev
    log_sd.e = log(sd.e),
    log_sd.c = log(sd.c)
  )

#cal log-md and log-smd and SE
bio_metareg <- bio_metareg %>%
  mutate(
    # Calculate Mean Difference (MD)
    MD = log_mean.e - log_mean.c,
    # Calculate Standard Error for MD (MD_SE)
    MD_SE = sqrt((log_sd.e^2 / n.e) + (log_sd.c^2 / n.c)),
    # Calculate Standardized Mean Difference (SMD) using pooled SD
    SMD = (log_mean.e - log_mean.c) / sqrt(((n.e - 1) * log_sd.e^2 + (n.c - 1) * log_sd.c^2) / (n.e + n.c - 2)),
    # Calculate Standard Error for SMD (SMD_SE)
    SMD_SE = sqrt((1 / n.e) + (1 / n.c) + ((SMD^2) / (2 * (n.e + n.c)))),
    SMD_var = SMD_SE^2
  )

ddimer <-subset(bio_metareg, Outcome == "D-dimer")
CRP <-subset(bio_metareg, Outcome == "CRP")


# Fit the Bayesian model using brms
model_metareg <- brm(
  formula = SMD | se(SMD_SE) ~  age_diff + (1 | Study),#for meta-regression model
  data = CRP,  # Your prepared data
  family = gaussian(),  # Gaussian distribution for the continuous log-transformed IRR
  prior = c(
    prior(normal(0,1), class = "Intercept"),  # Prior for the intercept
    prior(cauchy(0, 0.5), class = "sd")# Prior for the random effect SD
  ),
  chains = 4,  # Number of chains for MCMC sampling
  iter = 2000,  # Number of iterations per chain
  warmup = 1000,  # Number of warm-up iterations (burn-in)
  control = list(adapt_delta = 0.999), # Control parameter for better convergence
  save_pars = save_pars(all = TRUE)
)

summary(model_metareg)

# Get predicted values with uncertainty
newdata <- CRP %>% 
  mutate(pred = fitted(model_metareg, newdata = ., re_formula = NULL, summary = TRUE)[, "Estimate"],
         lower = fitted(model_metareg, newdata = ., re_formula = NULL, summary = TRUE)[, "Q2.5"],
         upper = fitted(model_metareg, newdata = ., re_formula = NULL, summary = TRUE)[, "Q97.5"])

# Plot categorical variables
ggplot(newdata, aes(x = AC, y = log_IRR)) +
  geom_point(aes(size = 1 / SE_log_IRR), alpha = 0.7) +  # effect sizes
  geom_point(aes(y = pred), color = "blue", shape = 18, size = 3) +  # predicted values
  geom_text(aes(label = Study), hjust = -0.1, vjust = 0.5, size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "blue") +  # CIs
  labs(
    title = "Meta-regression of log(IRR) by Anticoagulant Type",
    y = "log(IRR)",
    x = "Anticoagulant Type (AC)"
  ) +
  theme_minimal()

# Get predictions with uncertainty
metareg_data_pred <- CRP %>%
  mutate(predicted = fitted(model_metareg, newdata = ., re_formula = NA)[, "Estimate"],
         lower = fitted(model_metareg, newdata = ., re_formula = NA)[, "Q2.5"],
         upper = fitted(model_metareg, newdata = ., re_formula = NA)[, "Q97.5"])

# Plot numerical variables
ggplot(metareg_data_pred, aes(x = age_diff, y = log_IRR)) +
  geom_point(aes(size = 1 / SE_log_IRR), alpha = 0.7) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_text_repel(aes(label = Study), size = 3) +  # study labels
  labs(
    title = "Meta-regression: log(IRR) vs Mean Age Difference",
    x = "Mean Age  Difference",
    y = "log(IRR)"
  ) +
  theme_minimal()

#------------------------------------------------------------------------------#
# Funnel plot for publication bias for biomarkers
#------------------------------------------------------------------------------#

##Funnel plot for publication bias

pooled<--0.85 ## pooled estmate of outcome of interest
# Calculate the pseudo 95% CI lines for the funnel using outcome dataset of interest
se_seq <- seq(0, max(vte_data$SE_log_IRR), length.out = 100)
upper <- pooled + 1.96 * se_seq
lower <- pooled - 1.96 * se_seq

#Create a separate data frame for the funnel lines
funnel_lines <- data.frame(
  x = c(upper, lower),
  y = rep(se_seq, 2),
  group = rep(c("upper", "lower"), each = 100)
)
png(file = "funnel_bayes_rvte.png", width = 5000, height = 3000, res = 300)

ggplot(CRP, aes(x = log_IRR, y = SE_log_IRR)) +
  geom_point() +
  geom_vline(xintercept = pooled, linetype = "dashed", color = "blue") +
  geom_line(data = funnel_lines, aes(x = x, y = y, group = group), color = "red") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = Study), vjust = -0.5, size = 3) +
  scale_y_reverse() +
  labs(x = "Effect Size", y = "Standard Error")

dev.off()

#------------------------------------------------------------------------------#
# Leave one out analysis for biomarkers
#------------------------------------------------------------------------------#

# Function to fit brms model on data subset and extract summaries
fit_and_extract <- function(data_subset, full_data_size) {
  cat("Fitting model on data with", nrow(data_subset), "rows...\n")
  
  model <- brm(
    formula = SMD | se(SMD_SE) ~ 1 + (1 | Study),
    data = data_subset,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "Intercept"),
      prior(cauchy(0, 0.5), class = "sd")
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.999),
    save_pars = save_pars(all = TRUE),
    refresh = 0  # suppress sampling progress output
  )
  
  # Extract posterior draws and rename for clarity
  post.samples <- as_draws_df(model) %>%
    rename(
      smd = b_Intercept,
      tau = sd_Study__Intercept
    )
  
  cat("Model fitting done. Extracting summaries...\n")
  
  # Posterior probability pooled effect < 0
  p_less_0 <- mean(post.samples$smd < 0)
  
  # Calculate posterior I²
  post.samples$I2 <- (post.samples$tau^2) / 
    (post.samples$tau^2 + mean(data_subset$SMD_SE^2)) * 100
  
  I2_q <- quantile(post.samples$I2, probs = c(0.025, 0.5, 0.975))
  
  # Skewness of the effect size posterior
  skew_val <- skewness(post.samples$smd)
  
  # Summaries for intercept including 95% credible intervals
  smd_mean <- mean(post.samples$smd)
  smd_lower <- quantile(post.samples$smd, 0.025)
  smd_upper <- quantile(post.samples$smd, 0.975)
  
  # Mean tau (between-study heterogeneity)
  tau_mean <- mean(post.samples$tau)
  
  # Study label: "All" if full data, else LOO name of excluded study
  study_label <- ifelse(
    nrow(data_subset) == full_data_size,
    "All",
    paste0("LOO: ", setdiff(unique(ddimer$Study), unique(data_subset$Study)))
  )
  
  cat("Summaries for", study_label, "computed.\n\n")
  
  tibble(
    Study = study_label,
    Intercept = smd_mean,
    Lower_CI = smd_lower,
    Upper_CI = smd_upper,
    Skewness = skew_val,
    P_RR_less_0 = p_less_0,
    Tau = tau_mean,
    I2_lower = I2_q[1],
    I2_median = I2_q[2],
    I2_upper = I2_q[3]
  )
}

# Get study names for LOO
studies <- unique(ddimer$Study)
full_data_size <- nrow(ddimer)

cat("Starting full model fit with all data...\n")

# Fit full model once (no study left out)
full_result <- fit_and_extract(ddimer, full_data_size)

# Initialize list to store LOO results
loo_results <- vector("list", length(studies))

# Loop over studies
for (i in seq_along(studies)) {
  cat("Running LOO iteration leaving out study:", studies[i], "\n")
  
  # Subset data excluding current study
  loo_data <- subset(ddimer, Study != studies[i])
  
  # Fit model and extract summaries
  loo_results[[i]] <- fit_and_extract(loo_data, full_data_size)
}

cat("LOO fitting completed. Combining results...\n")

# Combine all results: full data + LOO iterations
all_results <- bind_rows(full_result, bind_rows(loo_results))


cat("Final results:\n")
print(all_results)

# Optionally, to save or export the results table here
# write.csv(all_results, "loo_results_summary.csv", row.names = FALSE)


