## Load packages
library(tidyverse)
library(ggplot2)
library(dplyr)

## Read dataset
prs_data <- read_csv("PRS_March2021.csv") # PRS data
brca1 <- read_csv("brca1.csv")
brca1 <- read_csv("brca2.csv") # Comment out when working on BRCA1

## Merge data by individual id
prs_data2 <- prs_data %>%
  mutate(indID = Person_ID) %>%
  dplyr::select(-"Person_ID")
brca1_prs <- left_join(brca1, prs_data2, by = "indID")
brca2_prs <- left_join(brca2, prs_data2, by = "indID")

brca1_prs <- brca1_prs |>
  mutate(df = df1) |>
  mutate(time = timeBC,
         status = BC)

##################################################################
####### Pedigree tree - who is present who is missing (PRS) ######
##################################################################
brca1_pedtree_test <- brca1_prs[1:6,] %>%
  mutate(gender = "female", 
         PRS_exist = ifelse(is.na(PRS), 0, 1),
         indID = as.character(indID), 
         famID = as.character(famID),
         fatherID = as.character(fatherID),
         motherID = as.character(motherID))

unique_fatherID <- unique(na.omit(brca1_pedtree_test$fatherID))
rows_to_add <- data.frame(indID = unique_fatherID, gender = "male",
                          famID = brca1_pedtree_test$famID, proband = 0,
                          affect = 0, PRS_exist = 0)
first_family <- bind_rows(brca1_pedtree_test, rows_to_add)[-1]
first_family <- first_family %>%
  distinct(indID, .keep_all = TRUE) %>%
  dplyr::select(c(indID, gender, famID, proband, affect, PRS_exist, fatherID, motherID)) 

library(kinship2)
ped_fam1 <- pedigree(id = first_family$indID,
                     dadid = first_family$fatherID,
                     momid = first_family$motherID,
                     sex = first_family$gender,
                     famid = first_family$famID, 
                     affected = cbind(first_family$affect, first_family$PRS_exist, first_family$proband))

## First family pedigree tree to visualize the missing PRS information
ped_toplot <- ped_fam1['10000001']
plot.pedigree(ped_toplot)

######################################################################
check_missingness <- glm(miss ~ mgeneI, data = brca1_prs)
check_missingness <- glm(miss ~ timeBC, data = brca1_prs)
check_missingness <- glm(miss ~ BC, data = brca1_prs)
check_missingness <- glm(miss ~ timeBC*BC, data = brca1_prs)
check_missingness <- glm(miss ~ proband, data = brca1_prs)
check_missingness <- glm(miss ~ mgeneI + timeBC + BC + proband, data = brca1_prs)
summary(check_missingness)


######################################################################

brca1_pedtree_test <- brca1_prs[1:6,] %>%
  mutate(gender = "female", 
         PRS_exist = ifelse(is.na(PRS), 0, 1),
         indID = as.character(indID), 
         famID = as.character(famID),
         fatherID = as.character(fatherID),
         motherID = as.character(motherID))

unique_fatherID <- unique(na.omit(brca1_pedtree_test$fatherID))
rows_to_add <- data.frame(indID = unique_fatherID, gender = "male",
                          famID = brca1_pedtree_test$famID, proband = 0,
                          BC = 0, PRS_exist = 0, mgeneI = 0)
first_family <- bind_rows(brca1_pedtree_test, rows_to_add)[-1]
first_family <- first_family %>%
  distinct(indID, .keep_all = TRUE) %>%
  dplyr::select(c(indID, gender, famID, proband, BC, PRS_exist, mgeneI, fatherID, motherID, mgeneI)) 

library(kinship2)
ped_fam1 <- pedigree(id = first_family$indID,
                     dadid = first_family$fatherID,
                     momid = first_family$motherID,
                     sex = first_family$gender,
                     famid = first_family$famID, 
                     affected = cbind(first_family$BC, first_family$PRS_exist, first_family$proband, first_family$mgeneI))

## First family pedigree tree to visualize the missing PRS information
ped_toplot <- ped_fam1['10000001']
plot.pedigree(ped_toplot)
##################################################################
####### ---------------End of Pedigree tree--------------- #######
##################################################################

## Missing mechanism specification
brca1_prs <- brca1_prs |>
  mutate(miss_index = ifelse(is.na(PRS) == TRUE, 1, 0))
model_miss <- glm(miss_index ~ mgeneI + log(timeBC) + currentage + proband, data = brca1_prs, family = binomial)
summary(model_miss)

brca1_prs_test1 <- brca1_prs |>
  mutate(gender = PRS)

######################

pen_model <- penmodel(Surv(time, status) ~ mgene + newx + gender, cluster = "famID", gvar = "mgene", data = famx, base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, design = "pop", parms = c(1/41.41327,1,0,0, 0, 1))

# Family data simulated from population-based design using a Weibull baseline hazard 

set.seed(4321)
fam <- simfam(N.fam = 200, design = "pop+", variation = "none", base.dist = "Weibull", 
              base.parms = c(0.01, 3), vbeta = c(-1.13, 2.35), agemin = 20, allelefreq = 0.02)

# Penetrance model fit for simulated family data

fit <- penmodel(Surv(time, status) ~ gender + mgene, cluster = "famID", design = "pop+",
                parms = c(0.01, 3, -1.13, 2.35, 1), data = fam, base.dist = "Weibull", frailty.dist = "lognormal")
summary(fit)

################################
# Plots can be used for thesis #
################################

## MCAR, MAR, MNAR
set.seed(123)
x = rnorm(3000, mean = 0, sd = 1)
y = 2*x + 0.5 + rnorm(3000, mean = 0, sd = 1)
data_miss <- data.frame(
  x = x,
  y = y
)

data_MCAR <- data_miss |> mutate(y = ifelse(runif(n()) < 0.6, NA, y))
data_MAR <- data_miss |> 
  mutate(y = ifelse(runif(n()) < 0.9 & x > 0, NA, y)) |>
  mutate(y = ifelse(runif(n()) < 0.6, NA, y))
data_MNAR <- data_miss |> 
  mutate(y = ifelse(runif(n()) < 0.9 & y > 0, NA, y)) |>
  mutate(y = ifelse(runif(n()) < 0.6, NA, y))

MCAR_plot <- ggplot(data_MCAR, aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = "Fitted Line")) +
  geom_line(stat = "smooth", data = data_miss, alpha = 0.5, aes(color = "True Line")) +
  ggtitle("Missing Completely at Random") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(color = "Line Type") +
  scale_color_manual(values = c("red" = "True Line", "blue" = "Fitted Line"))
  
MAR_plot <- ggplot(data_MAR, aes(x, y)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = "Fitted Line")) +
  geom_line(stat = "smooth", data = data_miss, alpha = 0.5, aes(color = "True Line")) +
  ggtitle("Missing at Random") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(color = "Line Type") +
  scale_color_manual(values = c("red" = "True Line", "blue" = "Fitted Line"))

MNAR_plot <- ggplot(data_MNAR, aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = "Fitted Line")) +
  geom_line(stat = "smooth", data = data_miss, alpha = 0.5, aes(color = "True Line")) +
  ggtitle("Missing Not at Random") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(color = "Line Type") +
  scale_color_manual(values = c("red" = "True Line", "blue" = "Fitted Line"))

gridExtra::grid.arrange(MCAR_plot, MAR_plot, MNAR_plot, ncol = 3)

data_combined <- bind_rows(
  data_miss |> mutate(type = "No Missing"),
  data_MCAR |> mutate(type = "MCAR"),
  data_MAR |> mutate(type = "MAR"),
  data_MNAR |> mutate(type = "MNAR")
)

density_plot <- ggplot(data_combined, aes(x = y, color = type, fill = type)) +
  geom_density(alpha = 0.3) +
  ggtitle("Density Plot of y for Different Missing Data Mechanisms") +
  theme_minimal() +
  labs(y = "Density") +
  scale_fill_manual(values = c("No Missing" = "blue", "MCAR" = "red", "MAR" = "green", "MNAR" = "purple")) +
  scale_color_manual(values = c("No Missing" = "blue", "MCAR" = "red", "MAR" = "green", "MNAR" = "purple")) +
  guides(
    fill = guide_legend(title = "Type"),
    color = guide_legend(title = "Type")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(density_plot)

##################################################################
####### Pedigree tree - who is present who is missing (PRS) ######
##################################################################
brca1_pedtree_test <- check_data[1:6,] %>%
  mutate(gender = "female", 
         PRS_exist = ifelse(is.na(newx), 0, 1),
         indID = as.character(indID), 
         famID = as.character(famID),
         fatherID = as.character(fatherID),
         motherID = as.character(motherID),
         fatherID = ifelse(fatherID == 0, NA, fatherID), 
         motherID = ifelse(motherID == 0, NA, motherID))

unique_fatherID <- unique(na.omit(brca1_pedtree_test$fatherID))
rows_to_add <- data.frame(indID = unique_fatherID, gender = "male",
                          famID = unique(brca1_pedtree_test$famID), proband = 0,
                          PRS_exist = 0, mgene = 0)
first_family <- bind_rows(brca1_pedtree_test, rows_to_add)
first_family <- first_family %>%
  distinct(indID, .keep_all = TRUE) %>%
  dplyr::select(c(indID, gender, famID, proband, PRS_exist, fatherID, motherID, mgene, status)) 
first_family <- first_family |>
  mutate(status = ifelse(gender == "male", 0, status))

library(kinship2)
ped_fam1 <- pedigree(id = first_family$indID,
                     dadid = first_family$fatherID,
                     momid = first_family$motherID,
                     sex = first_family$gender,
                     famid = first_family$famID, 
                     affected = cbind(first_family$status, first_family$proband, first_family$mgene,
                                      first_family$PRS_exist))

## First family pedigree tree to visualize the missing PRS information
ped_toplot <- ped_fam1['1']
plot.pedigree(ped_toplot)
pedigree.legend(ped_toplot, location = "topleft",
                labels = list("Disease Status", "Proband", "Mutation Gene Status", "PRS Observed"),
                cex = 0.6, radius = 0.2
)
##################################################################
####### ---------------End of Pedigree tree--------------- #######
##################################################################

########################## Variants of logistic distribution function ###########################
# Load necessary library
library(ggplot2)

# Create a sequence of values for the x-axis
x_values <- seq(-5, 5, length.out = 100)

# Define the logistic functions for each type with additional spacing
logistic_function <- function(x, alpha, beta) {
  1 / (1 + exp(-alpha * (x - beta)))
}

# Define the parameters for each type
alpha <- 1  # Slope parameter
right_beta <- 0    # Shift to the right
left_beta <- 0     # Shift to the left
mid_beta <- 1     # Centered
tail_beta <- 1     # Centered

# Generate y-values for each type with some shift to avoid overlap
right_y <- logistic_function(x_values, alpha, right_beta)
left_y <- 1 - logistic_function(x_values, alpha, left_beta)  # Inverted logistic function
mid_y <- logistic_function(abs(x_values), alpha, mid_beta)   # Symmetric around zero
tail_y <- 1 - logistic_function(abs(x_values), alpha, tail_beta)  # Inverted symmetric

# Combine into a data frame for plotting
data <- data.frame(
  x = rep(x_values, 4),
  y = c(right_y, left_y, mid_y, tail_y),
  type = factor(rep(c("RIGHT", "LEFT", "MID", "TAIL"), each = length(x_values)))
)

# Plot using ggplot2 with different linetypes and avoid overlapping
ggplot(data, aes(x = x, y = y, linetype = type)) +
  geom_line(size = 0.5) +
  labs(title = "Logistic Distribution Functions Variants",
       x = "WSS",
       y = "Probability of Missing",
       linetype = "Type") +
  scale_linetype_manual(values=c("solid", "dashed", "dotdash", "dotted")) +
  theme_minimal() +
  theme(element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(legend.position = "bottom")

##################################################################
###################### Plot the penetrance ######################
##################################################################

##########
penf <- function(est, x, age, base.dist="Weibull", frailty.dist=NULL, agemin=18, cuts=NULL){
  nbase <- length(est)-length(x)
  base.est <- exp(est[1:nbase])
  if(base.dist=="lognormal") base.est[1] <- est[1] 
  if(is.null(frailty.dist) || frailty.dist=="none"){
    xbeta <- sum(est[-c(1:nbase)]*x)
    H <- cumhaz(base.dist, age-agemin, base.est, cuts=cuts)*exp(xbeta) 
    pen <- 1-exp(-H)
  }
  else{    
    k <- est[nbase+1]
    xbeta <- sum(est[-c(1:(nbase+1))]*x)
    H <- cumhaz(base.dist, age-agemin, base.est, cuts=cuts)*exp(xbeta) 
    pen <- 1-laplace(dist=frailty.dist, g=H, k=k)
  }
  #  names(pen) <- age
  return(pen)
}

frailty_dist_data <- "gamma"
frailty_dist_model <- "lognormal"

est_true <- c(-4.71, 0.804, 2.13, 1.00)
est_kappa <- log(1)

est_CCA <- est_true + c(48.7, -2, -31, -17.8)/100
est_CCA_kappa <- est_kappa + 72.7/100

est_MI_LOGT <- est_true + c(-4.6, 7.4, -4.4, 83.2)/100
est_MI_LOGT_kappa <- est_kappa + 133.8/100

est_MI_H0T <- est_true + c(-6.5, 7.1, -4.5, 80.8)/100
est_MI_H0T_kappa <- est_kappa + 134/100

est_MI_K_LOGT <- est_true + c(2.8, -2.6, -2.6, -10.4)/100
est_MI_K_LOGT_kappa <- est_kappa + 23.4/100

est_MI_K_H0T <- est_true + c(4.2, -3.1, -5.2, -9.3)/100
est_MI_K_H0T_kappa <- est_kappa + 23.1/100

est_true <- c(est_true, exp(est_kappa))
est_CCA <- c(est_CCA, exp(est_CCA_kappa))
est_MI_LOGT <- c(est_MI_LOGT, exp(est_MI_LOGT_kappa))
est_MI_H0T <- c(est_MI_H0T, exp(est_MI_H0T_kappa))
est_MI_K_LOGT <- c(est_MI_K_LOGT, exp(est_MI_K_LOGT_kappa))
est_MI_K_H0T <- c(est_MI_K_H0T, exp(est_MI_K_H0T_kappa))


age <- seq(18, 70, by=1)
data <- expand.grid(age=age, mgene=c(0, 1), PRS=c(-2, 0,2))

## est_true
data$penetrance_true <- mapply(function(mgene, PRS, age) {
  penf(est = est_true, c(mgene, PRS), age, frailty.dist = frailty_dist_data, agemin = 18)
}, data$mgene, data$PRS, data$age)
## est_CCA
data$penetrance_CCA <- mapply(function(mgene, PRS, age) {
  penf(est = est_CCA, c(mgene, PRS), age, frailty.dist = frailty_dist_model, agemin = 18)
}, data$mgene, data$PRS, data$age)
## est_MI_LOGT
data$penetrance_MI_LOGT <- mapply(function(mgene, PRS, age) {
  penf(est = est_MI_LOGT, c(mgene, PRS), age, frailty.dist = frailty_dist_model, agemin = 18)
}, data$mgene, data$PRS, data$age)
## est_MI_H0T
data$penetrance_MI_H0T <- mapply(function(mgene, PRS, age) {
  penf(est = est_MI_H0T, c(mgene, PRS), age, frailty.dist = frailty_dist_model, agemin = 18)
}, data$mgene, data$PRS, data$age)
## est_MI_K_LOGT
data$penetrance_MI_K_LOGT <- mapply(function(mgene, PRS, age) {
  penf(est = est_MI_K_LOGT, c(mgene, PRS), age, frailty.dist = frailty_dist_model, agemin = 18)
}, data$mgene, data$PRS, data$age)
## est_MI_K_H0T
data$penetrance_MI_K_H0T <- mapply(function(mgene, PRS, age) {
  penf(est = est_MI_K_H0T, c(mgene, PRS), age, frailty.dist = frailty_dist_model, agemin = 18)
}, data$mgene, data$PRS, data$age)


## Combine the data for plotting
data_long <- data |>
  pivot_longer(cols = c("penetrance_true", "penetrance_CCA", "penetrance_MI_LOGT",
                        "penetrance_MI_H0T", "penetrance_MI_K_LOGT", "penetrance_MI_K_H0T"),
               names_to = "Methods",
               values_to = "Penetrance")

data_long$Methods <- factor(data_long$Methods, levels = c("penetrance_true", "penetrance_CCA",
                                                          "penetrance_MI_LOGT", "penetrance_MI_H0T",
                                                          "penetrance_MI_K_LOGT", "penetrance_MI_K_H0T"),
                             labels = c("True Value", "CCA", "MI-LOGT",
                                        "MI-H0T", "MI-K-LOGT", "MI-K-H0T"))

custom_colors <- c("True Value" = "black",   
                   "CCA" = "#9467bd",           
                   "MI-LOGT" = "green",       
                   "MI-H0T" = "blue",        
                   "MI-K-LOGT" = "#ff7f0e",     
                   "MI-K-H0T" = "red")      
custom_alpha <- c("True Value" = 0.8,   
                  "CCA" = 0.5,           
                  "MI-LOGT" = 0.5,       
                  "MI-H0T" = 0.5,        
                  "MI-K-LOGT" = 0.5,     
                  "MI-K-H0T" = 0.5)  
custom_linetype <- c("2" = "solid", 
                     "0"="dashed",
                     "-2" = "dotted")

data_long <- data_long |>
  mutate(mgene = factor(ifelse(mgene == 1, "Mutation Carrier", "Non Mutation Carrier"),
                        levels = c("Non Mutation Carrier", "Mutation Carrier")))


brca1_gammaCCA <- ggplot(data_long, aes(x = age, y = Penetrance, 
                                        color = Methods, linetype = factor(PRS), 
                                        alpha = Methods)) +
  geom_line(size = 0.9) +
  labs(title = "Penetrance Plot for Simulations",
       subtitle = "Matched - Log Normal Frailty Data, 80% Missing, High Familial Correlation",
       x = "Age",
       y = "Penetrance",
       color = "Methods") +
  scale_color_manual(values = custom_colors) +
  scale_alpha_manual(values = custom_alpha) +
  scale_linetype_manual("PRS",values = custom_linetype) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "bold")) +
  facet_wrap(~ mgene)

print(brca1_gammaCCA)


###########################################
library(ggplot2)

Bias_CCA <- c(-2.4, -5.5, -8.9, -13, -2.8, -7.1, -12.7, -17.8)
Bias_MI_LOGT <- c(-3.6, 4.6, 22.3, 81.3, -5.5, 1.8, 20.9, 83.2)
Bias_MI_H0T <- c(-3.7, 6, 22.5, 83.7, -5.3, 2.3, 19.5, 80.8)
Bias_MI_K_LOGT <- c(1.6, 4, 7, 10.8, -0.2, -1.5, -5, -8.3)
Bias_MI_K_H0T <- c(-2.9, -2.1, -1.6, -6.5, -0.5, -2.5, -7.1, -9.3)

df <- data.frame(
  Methods = rep(c("CCA", "MI-LOGT", "MI-H0T", "MI-K-LOGT", "MI-K-H0T"), each = 8), 
  Correlation = rep(c(rep("Low", 4), rep("High", 4)), times = 5), 
  Missing_Proportion = factor(rep(c("20%", "40%", "60%", "80%"), times = 10)), 
  Bias = c(Bias_CCA, Bias_MI_LOGT, Bias_MI_H0T, Bias_MI_K_LOGT, Bias_MI_K_H0T)
)

df$Methods <- factor(df$Methods, levels = c("CCA", "MI-LOGT", "MI-H0T", "MI-K-LOGT", "MI-K-H0T"))
blue_gradient <- colorRampPalette(c("lightblue", "darkblue"))(length(levels(df$Missing_Proportion)))

ggplot(df, aes(x = Methods, y = Bias, fill = Missing_Proportion)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  facet_wrap(~ Correlation, scales = "free_x", labeller = labeller(Correlation = c("Low" = "Low Correlation", "High" = "High Correlation"))) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dashed", color = "black") + 
  scale_y_continuous(limits = c(-100, 100)) +
  scale_fill_manual(values = blue_gradient, name = "Missing Proportion") + 
  theme_minimal() +
  labs(x = "Methods", y = "Bias (%)",
       title = expression("Bias (%) of " * beta[PRS]), 
       subtitle = "Matched Log Normal Frailty Data") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) 
  










