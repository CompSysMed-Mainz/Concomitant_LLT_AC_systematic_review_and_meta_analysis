################################################################################
# Concomitant lipid-lowering and anticoagulation therapy versus                #
# anticoagulation monotherapy for venous thromboembolism prevention:           #
# a systematic review and Bayesian meta-analysis                               #
#Code for event incidence                                                                              #
#         Created: 24.07.2025                                                  #
#         Last Updated: 05.05.2026                                             #
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
vte_data <- read_xlsx("outcome_lltAC.xlsx", sheet="VTE_IRR")
pts_data <- read_xlsx("outcome_lltAC.xlsx", sheet="PTS_IRR")

# Calculate the log-transformed IRR and standard error (SE) for each study
##for VTE
vte_data <- vte_data %>%
  mutate(
    log_IRR = log(IRR),
    SE_log_IRR = (log(IRR_upper) - log(IRR_lower)) / 3.92
  )


pts_data$se_logRR <- sqrt(
  (1/pts_data$event.e) - (1/pts_data$n.e) +
    (1/pts_data$event.c) - (1/pts_data$n.c)
)
pts_data <- pts_data %>%
  mutate(
    log_RR = log(RR)
  )

#------------------------------------------------------------------------------#
# Bayesian Meta Analysis Model
#------------------------------------------------------------------------------#

# Fit the Bayesian model using brms
model <- brm(
  #by default, the se terms replace the sigma, resulting in sigma 0 and NA for Rhat and ESS
  formula = log_RR | se(se_logRR) ~  1+(1 | Study), #for PTS
  #formula = log_IRR | se(SE_log_IRR) ~  1 + (1 | Study),#for VTE
  data = pts_data,  # Your prepared data - pts_data or vte_data
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
  (post.samples$tau^2 + mean(vte_data$SE_log_IRR^2)) * 100

# For calculating I² from the posterior for PTS data
post.samples$I2 <- (post.samples$tau^2) / 
  (post.samples$tau^2 + mean(pts_data$se_logRR^2)) * 100

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
  left_join(select(pts_data, Study, AC, LLT, Design), by = "Study")
forest.data.summary<- forest.data.summary%>%
  left_join(select(skewness_table, Study, Skewness), by = "Study")

##complete pooled effect data
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "AC"] <- "ACs"
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "LLT"] <- "Statins"
forest.data.summary[forest.data.summary$Study == "Pooled Effect", "Design"] <- " "

###plotting###
loadfonts(device = "win", quiet = TRUE) 

##convert to dataframe for forplo
##convert forest.data.summary to specific outcome VTE or PTS forest plot data summary
##eg forest.data.summary.ich<-as.data.frame(forest.data.summary)
forest.data.summary.vte<-as.data.frame(forest.data.summary)
#or
forest.data.summary.pts<-as.data.frame(forest.data.summary)

#------------------------------------------------------------------------------#
# Forest plots
#------------------------------------------------------------------------------#

##forest plot single outcome
#save image
png(file = "forest_rvte.png", width = 4000, height = 2000, res = 300)
forplo2(exp(forest.data.summary.major[,c(2:4)]), 
        normal.pdf = 1:nrow(forest.data.summary.major), ci.edge=F, row.labels = forest.data.summary.major$Study,
        normal.skew = forest.data.summary.major$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "Recurrent Venous Thromboembolism",
        add.columns = forest.data.summary.major[,8:10],
        add.colnames = c('AC','LLT', 'Design'),
        #favorlabs = "Combined therapy vs Monotherapy",
        diamond.pdf = nrow(forest.data.summary.major), diamond.col = 'red',
        column.spacing = 5, margin.right = 15, font = 'Corbel', shade.every = 1, shade.col = 4, em = 'RR',
        extra.line = "P(RR < 1) = 0.99      tau = 0.16      I² = 47%",extra.line.xshift = -2.97,extra.line.yshift = 0.25
)
dev.off()

##combined outcomes plots
##dataframes list
forest.data.summary.vte
forest.data.summary.pts

#combine dfs
events_combo<-rbind(
  cbind(forest.data.summary.vte, Outcome = "Recurrent Venous Thromboembolism"),
  cbind(forest.data.summary.pts, Outcome = "Post-Thrombotic Syndrome")
)

#assign groups with ID
events_combo <-cbind(events_combo, group_id = c(1,1,1,1,1,1,1,2,2,2,2))

##assign groups names
events_combo$group_name <- ifelse(
  events_combo$group_id == 1, "Recurrent Venous Thromboembolism",
  ifelse(
    events_combo$group_id == 2, "Post-Thrombotic Syndrome",
      NA  # (optional: handles unexpected values)
    )
  )
)

#plot
png(file = "forest_rvte_pts_incidence.png", width = 4000, height = 3000, res = 300)
forplo2(exp(events_combo[,c(2:4)]), 
        normal.pdf = 1:nrow(events_combo), ci.edge=F, row.labels = events_combo$Study,
        normal.skew = events_combo$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "Recurrent Venous Thromboembolism and Post-Thrombotic Syndrome Incidence",
        add.columns = events_combo[,8:10],
        add.colnames = c('AC','LLT','Design'),
        groups = as.numeric(factor(events_combo$group_id)),#group
        grouplabs = unique(events_combo$group_name),#group labels
        group.space = 2,
        diamond.pdf = c(7,11), diamond.col = 'red',
        column.spacing = 5, margin.right = 15,margin.left = 15, font = 'Corbel', shade.every = 1, shade.col = 4, em = 'SMD',
        extra.line = c("P(RR < 1) = 0.98      tau = 0.16      I² = 7%", "P(RR < 1) = 0.73      tau = 0.44      I² = 45%"),
        extra.line.y = c(9, 2),
        extra.line.xshift = -2.5
)
dev.off()

#------------------------------------------------------------------------------#
# Metaregression Model for rVTE
#------------------------------------------------------------------------------#


# Fit the Bayesian model using brms
model_metareg <- brm(
  formula = log_IRR | se(SE_log_IRR) ~  age_diff + (1 | Study),#for meta-regression model
  data = vte_data,  # Your prepared data
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
newdata <- vte_data %>% 
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
metareg_data_pred <- vte_data %>%
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
# Funnel plot for publication bias for rVTE
#------------------------------------------------------------------------------#

##Funnel plot for publication bias

pooled<--0.3 ## pooled estmate of outcome of interest
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

ggplot(vte_data, aes(x = log_IRR, y = SE_log_IRR)) +
  geom_point() +
  geom_vline(xintercept = pooled, linetype = "dashed", color = "blue") +
  geom_line(data = funnel_lines, aes(x = x, y = y, group = group), color = "red") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = Study), vjust = -0.5, size = 3) +
  scale_y_reverse() +
  labs(x = "Effect Size", y = "Standard Error")

dev.off()

#------------------------------------------------------------------------------#
# Leave one out analysis for rVTE
#------------------------------------------------------------------------------#
# Function to fit brms model on data subset and extract summaries
##change dataset to outcome of interest
fit_and_extract <- function(data_subset, full_data_size) {
  cat("Fitting model on data with", nrow(data_subset), "rows...\n")
  
  model <- brm(
    formula = log_IRR | se(SE_log_IRR) ~ 1 + (1 | Study),
    data = data_subset,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "Intercept"),
      prior(cauchy(0, 0.5), class = "sd")
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.995),
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
    (post.samples$tau^2 + mean(data_subset$SE_log_IRR^2)) * 100
  
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
    paste0("LOO: ", setdiff(unique(vte_data$Study), unique(data_subset$Study)))##change dataset to outcome of interest
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
#use dataset of outcome of interest
studies <- unique(vte_data$Study)
full_data_size <- nrow(vte_data)

cat("Starting full model fit with all data...\n")

# Fit full model once (no study left out)
full_result <- fit_and_extract(vte_data, full_data_size)

# Initialize list to store LOO results
loo_results <- vector("list", length(studies))

# Loop over studies
for (i in seq_along(studies)) {
  cat("Running LOO iteration leaving out study:", studies[i], "\n")
  
  # Subset data excluding current study
  loo_data <- subset(vte_data, Study != studies[i])
  
  # Fit model and extract summaries
  loo_results[[i]] <- fit_and_extract(loo_data, full_data_size)
}

cat("LOO fitting completed. Combining results...\n")

# Combine all results: full data + LOO iterations
all_results <- bind_rows(full_result, bind_rows(loo_results))

# Optional: reorder columns
all_results <- all_results %>%
  select(Study, Intercept, Lower_CI, Upper_CI, Skewness, P_RR_less_0, Tau, I2_lower, I2_median, I2_upper)

cat("Final results:\n")
print(all_results)

# Optionally, to save or export the results table here
# write.csv(all_results, "loo_results_summary.csv", row.names = FALSE)
