################################################################################
# Concomitant lipid-lowering and anticoagulation therapy versus                #
# anticoagulation monotherapy for venous thromboembolism prevention:           #
# a systematic review and Bayesian meta-analysis                               #
# Pharmacokinetics and pharmacodynamics                                        #
#         Created: 25.07.2025                                                  #
#         Last Updated: 20.08.2025                                             #
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

data <- read_xlsx("outcome_lltAC.xlsx", sheet="pkpd")

#seperate outcomes
AUC <- subset(data, Outcome == "AUC")
Cmax <- subset(data, Outcome == "Cmax")
Tmax <- subset(data, Outcome == "Tmax")
thalf <- subset(data, Outcome == "t1/2")
PT<- subset(data, Outcome == "PT")

##cal log smd for AUC, Cmax and PT
AUC<- AUC %>%
  mutate(
    # Calculate Mean
    log_mean.e = log(mean.e),
    log_mean.c = log(mean.c),
    # Calculate Standard dev
    log_sd.e = log(sd.e),
    log_sd.c = log(sd.c)
  )
AUC <- AUC %>%
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

##cal SMD for Tmax and Thalf
Tmax <- Tmax %>%
  mutate(
    # Calculate Mean Difference (MD)
    MD = mean.e - mean.c,
    # Calculate Standard Error for MD (MD_SE)
    MD_SE = sqrt((sd.e^2 / n.e) + (sd.c^2 / n.c)),
    # Calculate Standardized Mean Difference (SMD) using pooled SD
    SMD = (mean.e - mean.c) / sqrt(((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (n.e + n.c - 2)),
    # Calculate Standard Error for SMD (SMD_SE)
    SMD_SE = sqrt((1 / n.e) + (1 / n.c) + ((SMD^2) / (2 * (n.e + n.c)))),
    SMD_var = SMD_SE^2
  )

#------------------------------------------------------------------------------#
# Bayesian Meta Analysis Function
#------------------------------------------------------------------------------#
##run function
bma <- function(data){
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(brms)
  library(posterior)
  
  print("Starting model fit")
  model <- brm(
    #formula = SMD ~ 1 + (1 | Study),  # Log-transformed IRR with random effects for study,removed the se() term
    formula = SMD | se(SMD_SE) ~  1+(1 | Study), #by default, the se terms replace the sigma, resulting in sigma 0 and NA for Rhat and ESS
    data = data,  # Your prepared data
    family = gaussian(),  # Gaussian distribution for the continuous log-transformed IRR
    prior = c(
      prior(normal(0,1), class = "Intercept"),  # Prior for the intercept
      prior(cauchy(0,0.5), class = "sd")# Prior for the random effect SD
    ),
    chains = 4,  # Number of chains for MCMC sampling
    iter = 2000,  # Number of iterations per chain
    warmup = 1000,  # Number of warm-up iterations (burn-in)
    control = list(adapt_delta = 0.995), # Control parameter for better convergence
    save_pars = save_pars(all = TRUE)
  )
  print("Model fit complete")
  # Check the model summary
  print(summary(model))
  print("Starting posterior processing")
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
  print(smd.ecdf(0))
  
  # For calculating I² from the posterior for VTE data
  post.samples$I2 <- (post.samples$tau^2) / 
    (post.samples$tau^2 + mean(data$SMD_SE^2)) * 100
  
  
  # Summarize I²
  print(quantile(post.samples$I2, probs = c(0.025, 0.5, 0.975)))  # 95% credible interval for I²
  
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
    left_join(select(data, Study, AC, LLT, Design), by = "Study")##choose datatset major, ICH, GIB 
  forest.data.summary<- forest.data.summary%>%
    left_join(select(skewness_table, Study, Skewness), by = "Study")
  
  ##complete pooled effect data
  forest.data.summary[forest.data.summary$Study == "Pooled Effect", "AC"] <- "DOACs"
  forest.data.summary[forest.data.summary$Study == "Pooled Effect", "LLT"] <- "Statins"
  forest.data.summary[forest.data.summary$Study == "Pooled Effect", "Design"] <- " "
  
  
  ##convert to dataframe for forplo
  ##convert forest.data.summary to specific outcome major, ICH, GIB forest plot data summary
  ##eg forest.data.summary.ich<-as.data.frame(forest.data.summary)
  outcome <- AUC$Outcome[1]
  assign(paste0("forest.data.summary.", outcome), 
         as.data.frame(forest.data.summary))
  print("Completed forest.data.summary")
  return(forest.data.summary)
}

##call bma function to run bayesian meta analysis for outcome
PT_results <- bma(PT)
##or
AUC-results <- bma(AUC) ##and so on

###combo plot###combothalf plot
###combining plots
AUC_results
Cmax_results
Tmax_results
Thalf_results

combo_data <- rbind(
  cbind(AUC_results, pharmakin = "AUC"),
  cbind(Cmax_results, pharmakin = "Cmax"),
  cbind(Tmax_results, pharmakin = "Tmax"),
  cbind(Thalf_results, pharmakin = "Thalf")
)

#combo_data$group_id <- as.numeric(factor(interaction(combo_data$biomarker, combo_data$internal_group)))
combo_data <-cbind(combo_data, group_id = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4))

combo_data$group_name <- ifelse(
  combo_data$group_id == 1, "Area Under Concentration-Time Curve",
  ifelse(
    combo_data$group_id == 2, "Maximum Concentration",
    ifelse(
      combo_data$group_id == 3, "Time to Maximum Concentration",
      ifelse(
        combo_data$group_id == 4, "Terminal Half-life",
        NA  # (optional: handles unexpected values)
      )
    )
  )
)

##forest plot PT
#save image
png(file = "forest_pt.png", width = 4000, height = 2000, res = 300)
forplo2(PT_results[,c(2:4)], 
        normal.pdf = 1:nrow(PT_results), ci.edge=F, row.labels = PT_results$Study,
        normal.skew = PT_results$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "Prothrombin Time",
        add.columns = PT_results[,8:10],
        add.colnames = c('AC','LLT', 'Design'),
        #favorlabs = "Combined therapy vs Monotherapy",
        diamond.pdf = nrow(PT_results), diamond.col = 'red',
        linreg = TRUE,
        column.spacing = 5, margin.right = 15, font = 'Corbel', shade.every = 1, shade.col = 4, em = 'log SMD',
        extra.line = "P(SMD < 0) = 0.1      tau = 1.32      I² = 89%",extra.line.xshift = -0.2,extra.line.yshift = -0.1
)
dev.off()

png(file = "forest_pharmacokinetics.png", width = 4000, height = 3000, res = 300)
forplo2(combo_data[,c(2:4)], 
        normal.pdf = 1:nrow(combo_data), ci.edge=F, row.labels = combo_data$Study,
        normal.skew = combo_data$Skewness,
        normal.col = 4, normal.alpha = 0.2, 
        title = "Anticoagulant Pharmacokinetics",
        add.columns = combo_data[,8:10],
        add.colnames = c('AC','LLT', 'Design'),
        groups = as.numeric(factor(combo_data$group_id)),
        grouplabs = unique(combo_data$group_name),
        group.space = 3,
        #favorlabs = "Combined therapy vs Monotherapy",
        diamond.pdf = c(5,10,15,20), diamond.col = 'red',
        column.spacing = 5, margin.right = 15, font = 'Corbel',margin.left = 15, shade.every = 1, shade.col = 4, em = 'log SMD',
        extra.line = list("P(SMD < 0) = 0.51       τ = 0.19      I² = 16%","P(SMD < 0) = 0.53      τ = 0.19      I² = 16%",
                          list(label="SMD             95% CI            AC                        LLT                   Design", font =2), "P(SMD < 0) = 0.07     τ = 0.24      I² = 21%","P(SMD < 1) = 0.33      τ = 0.21      I² = 19%"),
        ##list for include new column headers for Tmax and Thalf, with non-log SMD
        extra.line.y = c(29.5, 20.5,18.5, 11.5, 2.5), 
        extra.line.xshift = -0.09
)
 dev.off()


