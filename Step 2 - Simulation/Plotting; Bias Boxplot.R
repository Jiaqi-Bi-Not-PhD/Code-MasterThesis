library(ggtext)
library(ggplot2)
library(Hmisc)

show_and_save <- function(plot, filename) {
  print(plot)
  ggsave(filename, plot = plot)
}

##################################################################################################################
################################## Gamma Well-Specified - High Corr ##############################################
##################################################################################################################

## 20% missing
Well_Gamma_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/20 Missing/Gamma/MAR20 Gamma WellSpec; Simulation Results.RData")
Well_Gamma_20 <- Well_Gamma_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_20$Method <- factor(Well_Gamma_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_Gamma_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_20 <- Well_Gamma_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_20$SimulationID))
tab_sim <- Well_Gamma_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_20_HighCorr <- Well_Gamma_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 20% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_Gamma_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/20 Missing/Gamma/Well_Gamma_20_HighCorr.jpeg")


## 40% missing
Well_Gamma_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/40 Missing/Gamma/MAR40 Gamma WellSpec; Simulation Results.RData")
Well_Gamma_40 <- Well_Gamma_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_40$Method <- factor(Well_Gamma_40$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_Gamma_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_40 <- Well_Gamma_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_40$SimulationID))
tab_sim <- Well_Gamma_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_40_HighCorr <- Well_Gamma_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 40% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_Gamma_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/40 Missing/Gamma/Well_Gamma_40_HighCorr.jpeg")


## 60% missing
Well_Gamma_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/60 Missing/Gamma/MAR60 Gamma WellSpec; Simulation Results.RData")
Well_Gamma_60 <- Well_Gamma_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))

Well_CCA <- Well_Gamma_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_60 <- Well_Gamma_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_60$SimulationID))
tab_sim <- Well_Gamma_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_60_HighCorr <- Well_Gamma_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 60% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_Gamma_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/60 Missing/Gamma/Well_Gamma_60_HighCorr.jpeg")


## 80% missing
Well_Gamma_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/80 Missing/Gamma/MAR80 Gamma WellSpec; Simulation Results.RData")
Well_Gamma_80 <- Well_Gamma_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_80$Method <- factor(Well_Gamma_80$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_Gamma_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_80 <- Well_Gamma_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_80$SimulationID))
tab_sim <- Well_Gamma_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_80_HighCorr <- Well_Gamma_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 80% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_Gamma_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/80 Missing/Gamma/Well_Gamma_80_HighCorr.jpeg")

##################################################################################################################
################################### Gamma Mis-Specified - High Corr ##############################################
##################################################################################################################

## 20% missing
Mis_Gamma_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/MAR20 GL Missspecified; Simulation Results.RData")
Mis_Gamma_20 <- Mis_Gamma_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_20$Method <- factor(Mis_Gamma_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_20 <- Mis_Gamma_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_Gamma_20$SimulationID))
tab_sim <- Mis_Gamma_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_20 <- Mis_Gamma_20 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_20_HighCorr <- Mis_Gamma_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 20% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/Mis_Gamma_20_HighCorr.jpeg")


## 40% missing
Mis_Gamma_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/MAR40 GL Misspecified; Simulation Results.RData")
Mis_Gamma_40 <- Mis_Gamma_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_40$Method <- factor(Mis_Gamma_40$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_40 <- Mis_Gamma_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_Gamma_40$SimulationID))
tab_sim <- Mis_Gamma_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_40 <- Mis_Gamma_40 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_40_HighCorr <- Mis_Gamma_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 40% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/Mis_Gamma_40_HighCorr.jpeg")


## 60% missing
Mis_Gamma_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/MAR60 GL Misspecified; Simulation Results.RData")
Mis_Gamma_60 <- Mis_Gamma_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_60$Method <- factor(Mis_Gamma_60$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_60 <- Mis_Gamma_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_Gamma_60$SimulationID))
tab_sim <- Mis_Gamma_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_60 <- Mis_Gamma_60 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_60_HighCorr <- Mis_Gamma_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 60% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/Mis_Gamma_60_HighCorr.jpeg")


## 80% missing
Mis_Gamma_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/MAR80 GL Misspecified; Simulation Results.RData")
Mis_Gamma_80 <- Mis_Gamma_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_80$Method <- factor(Mis_Gamma_80$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_80 <- Mis_Gamma_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_Gamma_80$SimulationID))
tab_sim <- Mis_Gamma_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_80 <- Mis_Gamma_80 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_80_HighCorr <- Mis_Gamma_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 80% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/Mis_Gamma_80_HighCorr.jpeg")


##################################################################################################################
############################################# LogNo Well-Specified - High Corr ###################################
##################################################################################################################
## 20% missing
Well_LogNo_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/20 Missing/Log Normal/MAR20 LogNormal WellSpec; Simulation Results.RData")
Well_LogNo_20 <- Well_LogNo_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_20$Method <- factor(Well_LogNo_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_20 <- Well_LogNo_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_LogNo_20$SimulationID))
tab_sim <- Well_LogNo_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_LogNo_20_HighCorr <- Well_LogNo_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 20% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_LogNo_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/20 Missing/Log Normal/Well_LogNo_20_HighCorr.jpeg")

## 40% missing
Well_LogNo_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/40 Missing/Log Normal/MAR40 LogNormal WellSpec; Simulation Results.RData")
Well_LogNo_40 <- Well_LogNo_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_40$Method <- factor(Well_LogNo_40$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_40 <- Well_LogNo_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_LogNo_40$SimulationID))
tab_sim <- Well_LogNo_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_LogNo_40_HighCorr <- Well_LogNo_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 40% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_LogNo_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/40 Missing/Log Normal/Well_LogNo_40_HighCorr.jpeg")


## 60% missing
Well_LogNo_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/60 Missing/Log Normal/MAR60 LogNormal WellSpec; Simulation Results.RData")
Well_LogNo_60 <- Well_LogNo_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_60$Method <- factor(Well_LogNo_60$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_60 <- Well_LogNo_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_LogNo_60$SimulationID))
tab_sim <- Well_LogNo_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_LogNo_60_HighCorr <- Well_LogNo_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 60% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_LogNo_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/60 Missing/Log Normal/Well_LogNo_60_HighCorr.jpeg")

## 80% missing
Well_LogNo_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/80 Missing/Log Normal/MAR80 LogNormal WellSpec; Simulation Results.RData")
Well_LogNo_80 <- Well_LogNo_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_80$Method <- factor(Well_LogNo_80$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_80 <- Well_LogNo_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_LogNo_80$SimulationID))
tab_sim <- Well_LogNo_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_LogNo_80_HighCorr <- Well_LogNo_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 80% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Well_LogNo_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Well-Specified/80 Missing/Log Normal/Well_LogNo_80_HighCorr.jpeg")

##################################################################################################################
############################################# LogNo Mis-Specified - High Corr ###################################
##################################################################################################################
## 20% missing
Mis_LogNo_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Log Normal Using Gamma/MAR20 LG Misspecified; Simulation Results.RData")
Mis_LogNo_20 <- Mis_LogNo_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_20$Method <- factor(Mis_LogNo_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_20 <- Mis_LogNo_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_LogNo_20$SimulationID))
tab_sim <- Mis_LogNo_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_20 <- Mis_LogNo_20 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_LogNo_20_HighCorr <- Mis_LogNo_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 20% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Mis_LogNo_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Log Normal Using Gamma/Mis_LogNo_20_HighCorr.jpeg")

## 40% missing
Mis_LogNo_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Log Normal Using Gamma/MAR40 LG Misspecified; Simulation Results.RData")
Mis_LogNo_40 <- Mis_LogNo_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_40$Method <- factor(Mis_LogNo_40$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_40 <- Mis_LogNo_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_LogNo_40$SimulationID))
tab_sim <- Mis_LogNo_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_40 <- Mis_LogNo_40 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_LogNo_40_HighCorr <- Mis_LogNo_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 40% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Mis_LogNo_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Log Normal Using Gamma/Mis_LogNo_40_HighCorr.jpeg")


## 60% missing
Mis_LogNo_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Log Normal Using Gamma/MAR60 LG Misspecified; Simulation Results.RData")
Mis_LogNo_60 <- Mis_LogNo_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_60$Method <- factor(Mis_LogNo_60$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_60 <- Mis_LogNo_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_LogNo_60$SimulationID))
tab_sim <- Mis_LogNo_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_60 <- Mis_LogNo_60 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_LogNo_60_HighCorr <- Mis_LogNo_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 60% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Mis_LogNo_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Log Normal Using Gamma/Mis_LogNo_60_HighCorr.jpeg")

## 80% missing
Mis_LogNo_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Log Normal Using Gamma/MAR80 LG Misspecified; Simulation Results.RData")
Mis_LogNo_80 <- Mis_LogNo_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_80$Method <- factor(Mis_LogNo_80$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_80 <- Mis_LogNo_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Mis_LogNo_80$SimulationID))
tab_sim <- Mis_LogNo_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_80 <- Mis_LogNo_80 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_LogNo_80_HighCorr <- Mis_LogNo_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Absolute Bias*100") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 80% Missing, High Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )

show_and_save(Mis_LogNo_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Log Normal Using Gamma/Mis_LogNo_80_HighCorr.jpeg")

##################################################################################################################
############################################# Gamma Well-Specified -  Low Corr ###################################
##################################################################################################################
## 20% missing
Well_Gamma_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/20 Missing/Gamma/MAR20 Gamma Well Specified; Simulation Results.RData")
Well_Gamma_20 <- Well_Gamma_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_20$Method <- factor(Well_Gamma_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))
Well_CCA <- Well_Gamma_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_20 <- Well_Gamma_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_20$SimulationID))
tab_sim <- Well_Gamma_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
  

## Boxplot
Well_Gamma_20_LowCorr <- Well_Gamma_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 20% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_Gamma_20_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/20 Missing/Gamma/Well_Gamma_20_LowCorr.jpeg")

## 40% missing
Well_Gamma_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/40 Missing/Gamma/MAR40 Gamma Well Specified; Simulation Results.RData")
Well_Gamma_40 <- Well_Gamma_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_40$Method <- factor(Well_Gamma_40$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))
Well_CCA <- Well_Gamma_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_40 <- Well_Gamma_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_40$SimulationID))
tab_sim <- Well_Gamma_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_40_LowCorr <- Well_Gamma_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 40% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_Gamma_40_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/40 Missing/Gamma/Well_Gamma_40_LowCorr.jpeg")

## 60% missing
Well_Gamma_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/60 Missing/Gamma/MAR60 Gamma Well Specified; Simulation Results.RData")
Well_Gamma_60 <- Well_Gamma_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_60$Method <- factor(Well_Gamma_60$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_Gamma_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_60 <- Well_Gamma_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_60$SimulationID))
tab_sim <- Well_Gamma_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_60_LowCorr <- Well_Gamma_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 60% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_Gamma_60_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/60 Missing/Gamma/Well_Gamma_60_LowCorr.jpeg")

## 80% missing
Well_Gamma_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/80 Missing/Gamma/MAR80 Gamma Well Specified; Simulation Results.RData")
Well_Gamma_80 <- Well_Gamma_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_Gamma_80$Method <- factor(Well_Gamma_80$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_Gamma_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_Gamma_80 <- Well_Gamma_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = Bias)
n_sim <- length(unique(Well_Gamma_80$SimulationID))
tab_sim <- Well_Gamma_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

## Boxplot
Well_Gamma_80_LowCorr <- Well_Gamma_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty, Well-Specified, 80% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_Gamma_80_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/80 Missing/Gamma/Well_Gamma_80_LowCorr.jpeg")
##################################################################################################################
############################################# LogNo Mis-Specified -  Low Corr ####################################
##################################################################################################################
## 20% missing
Mis_LogNo_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/20 Missing/Log Normal Using Gamma/MAR20 Log Normal Mis Specified; Simulation Results.RData")
Mis_LogNo_20 <- Mis_LogNo_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_20$Method <- factor(Mis_LogNo_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_20 <- Mis_LogNo_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_LogNo_20$SimulationID))
tab_sim <- Mis_LogNo_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_20 <- Mis_LogNo_20 |>
  filter(Parameters == "mgene" | Parameters == "PRS")

## Boxplot
Mis_LogNo_20_LowCorr <- Mis_LogNo_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 20% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_LogNo_20_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/20 Missing/Log Normal Using Gamma/Mis_LogNo_20_LowCorr.jpeg")

## 40% missing
Mis_LogNo_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/40 Missing/Log Normal Using Gamma/MAR40 Log Normal Mis Specified; Simulation Results.RData")
Mis_LogNo_40 <- Mis_LogNo_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_40$Method <- factor(Mis_LogNo_40$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_40 <- Mis_LogNo_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_LogNo_40$SimulationID))
tab_sim <- Mis_LogNo_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_40 <- Mis_LogNo_40 |>
  filter(Parameters == "mgene" | Parameters == "PRS")

## Boxplot
Mis_LogNo_40_LowCorr <- Mis_LogNo_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 40% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_LogNo_40_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/40 Missing/Log Normal Using Gamma/Mis_LogNo_40_LowCorr.jpeg")

## 60% missing
Mis_LogNo_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/60 Missing/Log Normal Using Gamma/MAR60 Log Normal Mis Specified; Simulation Results.RData")
Mis_LogNo_60 <- Mis_LogNo_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_60$Method <- factor(Mis_LogNo_60$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_60 <- Mis_LogNo_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_LogNo_60$SimulationID))
tab_sim <- Mis_LogNo_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_60 <- Mis_LogNo_60 |>
  filter(Parameters == "mgene" | Parameters == "PRS")

## Boxplot
Mis_LogNo_60_LowCorr <- Mis_LogNo_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 60% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_LogNo_60_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/60 Missing/Log Normal Using Gamma/Mis_LogNo_60_LowCorr.jpeg")

## 80% missing
Mis_LogNo_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/80 Missing/Log Normal Using Gamma/MAR80 Log Normal Mis Specified; Simulation Results.RData")
Mis_LogNo_80 <- Mis_LogNo_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_LogNo_80$Method <- factor(Mis_LogNo_80$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_LogNo_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_LogNo_80 <- Mis_LogNo_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_LogNo_80$SimulationID))
tab_sim <- Mis_LogNo_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_LogNo_80 <- Mis_LogNo_80 |>
  filter(Parameters == "mgene" | Parameters == "PRS")

## Boxplot
Mis_LogNo_80_LowCorr <- Mis_LogNo_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty Data, Mis-Specified, 80% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_LogNo_80_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/80 Missing/Log Normal Using Gamma/Mis_LogNo_80_LowCorr.jpeg")


##################################################################################################################
############################################# LogNo Well-Specified -  Low Corr ###################################
##################################################################################################################
## 20% missing
Well_LogNo_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/20 Missing/Log Normal/MAR20 Log Normal Well Specified; Simulation Results.RData")
Well_LogNo_20 <- Well_LogNo_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_20$Method <- factor(Well_LogNo_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_20 <- Well_LogNo_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Well_LogNo_20$SimulationID))
print( Well_LogNo_20 |> 
         dplyr::group_by(Parameters, Method) |> 
         dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                          Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                          MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                          Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                          RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                          MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                          Coverage = mean(Coverage, na.rm = TRUE),
                          EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                          Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
         mutate(percent_model_SE = (Ave.ModSE/EmpSE - 1)*100,
                Bias = Bias*100) |>
         dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) , n = 30)

## Boxplot
Well_LogNo_20_LowCorr <- Well_LogNo_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 20% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_LogNo_20_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/20 Missing/Log Normal/Well_LogNo_20_LowCorr.jpeg")

## 40% missing
Well_LogNo_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/40 Missing/Log Normal/MAR40 Log Normal Well Specified; Simulation Results.RData")
Well_LogNo_40 <- Well_LogNo_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_40$Method <- factor(Well_LogNo_40$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))
Well_CCA <- Well_LogNo_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_40 <- Well_LogNo_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Well_LogNo_40$SimulationID))
print( Well_LogNo_40 |> 
         dplyr::group_by(Parameters, Method) |> 
         dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                          Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                          MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                          Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                          RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                          MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                          Coverage = mean(Coverage, na.rm = TRUE),
                          EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                          Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
         mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
                Bias = Bias*100,
                Coverage = Coverage*100) |>
         dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) , n = 30)

## Boxplot
Well_LogNo_40_LowCorr <- Well_LogNo_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 40% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_LogNo_40_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/40 Missing/Log Normal/Well_LogNo_40_LowCorr.jpeg")

## 60% missing
Well_LogNo_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/60 Missing/Log Normal/MAR60 Log Normal Well Specified; Simulation Results.RData")
Well_LogNo_60 <- Well_LogNo_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_60$Method <- factor(Well_LogNo_60$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_60 <- Well_LogNo_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Well_LogNo_60$SimulationID))
print( Well_LogNo_60 |> 
         dplyr::group_by(Parameters, Method) |> 
         dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                          Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                          MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                          Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                          RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                          MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                          Coverage = mean(Coverage, na.rm = TRUE),
                          EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                          Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
         mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
                Bias = Bias*100,
                Coverage = Coverage*100) |>
         dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) , n = 30)

## Boxplot
Well_LogNo_60_LowCorr <- Well_LogNo_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 60% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_LogNo_60_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/60 Missing/Log Normal/Well_LogNo_60_LowCorr.jpeg")

## 80% missing
Well_LogNo_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/80 Missing/Log Normal/MAR80 Log Normal Well Specified; Simulation Results.RData")
Well_LogNo_80 <- Well_LogNo_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Well_LogNo_80$Method <- factor(Well_LogNo_80$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Well_CCA <- Well_LogNo_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Well_LogNo_80 <- Well_LogNo_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Well_LogNo_80$SimulationID))
tab_sim <- Well_LogNo_80 |> 
         dplyr::group_by(Parameters, Method) |> 
         dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                          Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                          MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                          Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                          RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                          MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                          Coverage = mean(Coverage, na.rm = TRUE),
                          EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                          Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
         mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
                Bias = Bias*100,
                Coverage = Coverage*100) |>
         mutate(Coverage = round(Coverage, 1)) |>
         dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
         pivot_wider(
           id_cols = Parameters,
           names_from = Method,
           values_from = c(Bias, percent_model_SE, Coverage),
           names_sep = "_"
         ) |>
         ungroup() |>
         dplyr::select(
           Parameters,
           Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
           `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
           `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
           `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
           `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
         ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
as.vector(unname(unlist(tab_sim[1, -1])))
as.vector(unname(unlist(tab_sim[2, -1])))
as.vector(unname(unlist(tab_sim[3, -1])))
as.vector(unname(unlist(tab_sim[4, -1])))
as.vector(unname(unlist(tab_sim[5, -1])))

## Boxplot
Well_LogNo_80_LowCorr <- Well_LogNo_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Log Normal Frailty, Well-Specified, 80% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Well_LogNo_80_LowCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Well-Specified/80 Missing/Log Normal/Well_LogNo_80_LowCorr.jpeg") ## Change number here

##################################################################################################################
############################################# Gamma Mis-Specified  - High Corr ###################################
##################################################################################################################
## 20% missing
Mis_Gamma_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/MAR20 GL Missspecified; Simulation Results.RData")
Mis_Gamma_20 <- Mis_Gamma_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_20$Method <- factor(Mis_Gamma_20$Method, 
                               levels = c("CCA", "MI-logT", 
                                          "MI-H0T", "MI-K-logT", 
                                          "MI-K-H0T"),
                               labels = c("CCA", "MI-logT", 
                                          "MI-H0T", "**MI-K-logT**", 
                                          "**MI-K-H0T**"))

Mis_Gamma_20 <- Mis_Gamma_20 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_20_HighCorr <- Mis_Gamma_20 |>
  ggplot(aes(x = Parameters, y = Bias, fill = Method)) +
  geom_boxplot(aes(group = interaction(Parameters, Method)), 
               position = position_dodge(width = 0.8), size = 0.2,
               width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  labs(fill = "Method") +
  theme_minimal() +
  theme(legend.text = element_markdown()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  ggtitle("Biases - Gamma Using LogNormal, Mis-Specified, 20% Missing, High Corr") + 
  theme(plot.title = element_text(size = 17, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_markdown(size = 14))
show_and_save(Mis_Gamma_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/Mis_Gamma_20_HighCorr.jpeg")

## 40% missing
Mis_Gamma_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/MAR40 GL Misspecified; Simulation Results.RData")
Mis_Gamma_40 <- Mis_Gamma_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_40$Method <- factor(Mis_Gamma_40$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Mis_Gamma_40 <- Mis_Gamma_40 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_40_HighCorr <- Mis_Gamma_40 |>
  ggplot(aes(x = Parameters, y = Bias, fill = Method)) +
  geom_boxplot(aes(group = interaction(Parameters, Method)), 
               position = position_dodge(width = 0.8), size = 0.2,
               width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  labs(fill = "Method") +
  theme_minimal() +
  theme(legend.text = element_markdown()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  ggtitle("Biases - Gamma Using LogNormal, Mis-Specified, 40% Missing, High Corr") + 
  theme(plot.title = element_text(size = 17, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_markdown(size = 14))
show_and_save(Mis_Gamma_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/Mis_Gamma_40_HighCorr.jpeg")

## 60% missing
Mis_Gamma_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/MAR60 GL Misspecified; Simulation Results.RData")
Mis_Gamma_60 <- Mis_Gamma_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_60$Method <- factor(Mis_Gamma_60$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Mis_Gamma_60 <- Mis_Gamma_60 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_60_HighCorr <- Mis_Gamma_60 |>
  ggplot(aes(x = Parameters, y = Bias, fill = Method)) +
  geom_boxplot(aes(group = interaction(Parameters, Method)), 
               position = position_dodge(width = 0.8), size = 0.2,
               width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  labs(fill = "Method") +
  theme_minimal() +
  theme(legend.text = element_markdown()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  ggtitle("Biases - Gamma Using LogNormal, Mis-Specified, 60% Missing, High Corr") + 
  theme(plot.title = element_text(size = 17, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_markdown(size = 14))
show_and_save(Mis_Gamma_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/Mis_Gamma_60_HighCorr.jpeg")

## 80% missing
Mis_Gamma_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/MAR80 GL Misspecified; Simulation Results.RData")
Mis_Gamma_80 <- Mis_Gamma_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_80$Method <- factor(Mis_Gamma_80$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Mis_Gamma_80 <- Mis_Gamma_80 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_80_HighCorr <- Mis_Gamma_80 |>
  ggplot(aes(x = Parameters, y = Bias, fill = Method)) +
  geom_boxplot(aes(group = interaction(Parameters, Method)), 
               position = position_dodge(width = 0.8), size = 0.2,
               width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  labs(fill = "Method") +
  theme_minimal() +
  theme(legend.text = element_markdown()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  ggtitle("Biases - Gamma Using LogNormal, Mis-Specified, 80% Missing, High Corr") + 
  theme(plot.title = element_text(size = 17, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_markdown(size = 14))
show_and_save(Mis_Gamma_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/High Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/Mis_Gamma_80_HighCorr.jpeg")


##################################################################################################################
############################################# Gamma Mis-Specified  -  Low Corr ###################################
##################################################################################################################
## 20% missing
Mis_Gamma_20 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/MAR20 Gamma Mis Specified; Simulation Results.RData")
Mis_Gamma_20 <- Mis_Gamma_20 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_20$Method <- factor(Mis_Gamma_20$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_20 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_20 <- Mis_Gamma_20 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_Gamma_20$SimulationID))
tab_sim <- Mis_Gamma_20 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_20 <- Mis_Gamma_20 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_20_HighCorr <- Mis_Gamma_20 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 20% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_20_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/20 Missing/Gamma Using LogNormal/Mis_Gamma_20_LowCorr.jpeg")

## 40% missing
Mis_Gamma_40 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/MAR40 Gamma Mis Specified; Simulation Results.RData")
Mis_Gamma_40 <- Mis_Gamma_40 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_40$Method <- factor(Mis_Gamma_40$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_40 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_40 <- Mis_Gamma_40 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_Gamma_40$SimulationID))
tab_sim <- Mis_Gamma_40 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_40 <- Mis_Gamma_40 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_40_HighCorr <- Mis_Gamma_40 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 40% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_40_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/40 Missing/Gamma Using LogNormal/Mis_Gamma_40_LowCorr.jpeg")

## 60% missing
Mis_Gamma_60 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/MAR60 Gamma Mis Specified; Simulation Results.RData")
Mis_Gamma_60 <- Mis_Gamma_60 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_60$Method <- factor(Mis_Gamma_60$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_60 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_60 <- Mis_Gamma_60 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_Gamma_60$SimulationID))
tab_sim <- Mis_Gamma_60 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_60 <- Mis_Gamma_60 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_60_HighCorr <- Mis_Gamma_60 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-Specified, 60% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_60_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/60 Missing/Gamma Using LogNormal/Mis_Gamma_60_LowCorr.jpeg")

## 80% missing
Mis_Gamma_80 <- readRDS("~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/MAR80 Gamma Mis Specified; Simulation Results.RData")
Mis_Gamma_80 <- Mis_Gamma_80 |>
  mutate(Parameters = recode(Parameters, "newx_I" = "PRS", "newx" = "PRS"))
Mis_Gamma_80$Method <- factor(Mis_Gamma_80$Method, 
                              levels = c("CCA", "MI-logT", 
                                         "MI-H0T", "MI-K-logT", 
                                         "MI-K-H0T"),
                              labels = c("CCA", "MI-logT", 
                                         "MI-H0T", "**MI-K-logT**", 
                                         "**MI-K-H0T**"))

Well_CCA <- Mis_Gamma_80 |>
  filter(Method == "CCA") |>
  dplyr::select(Parameters, true_value) |>
  arrange(Parameters) |>
  distinct(Parameters, .keep_all = TRUE) 

Mis_Gamma_80 <- Mis_Gamma_80 |>
  left_join(Well_CCA, by = "Parameters", suffix = c("", "_CCA")) |>
  mutate(true_value = ifelse(Method != "CCA", true_value_CCA, true_value)) |>
  dplyr::select(-true_value_CCA) |>
  mutate(Bias2 = (Estimates - true_value)/true_value)
n_sim <- length(unique(Mis_Gamma_80$SimulationID))
tab_sim <- Mis_Gamma_80 |> 
  dplyr::group_by(Parameters, Method) |> 
  dplyr::summarise(True_Value = mean(true_value, na.rm = TRUE),
                   Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias2, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 ) ) |>
  mutate(percent_model_SE = round((Ave.ModSE/EmpSE - 1)*100),
         Bias = Bias*100,
         Coverage = Coverage*100) |>
  mutate(Coverage = round(Coverage, 1)) |>
  dplyr::select(Parameters, Method, Bias, Coverage, percent_model_SE) |>
  pivot_wider(
    id_cols = Parameters,
    names_from = Method,
    values_from = c(Bias, percent_model_SE, Coverage),
    names_sep = "_"
  ) |>
  ungroup() |>
  dplyr::select(
    Parameters,
    Bias_CCA, percent_model_SE_CCA, Coverage_CCA,
    `Bias_MI-logT`, `percent_model_SE_MI-logT`, `Coverage_MI-logT`,
    `Bias_MI-H0T`, `percent_model_SE_MI-H0T`, `Coverage_MI-H0T`,
    `Bias_**MI-K-logT**`, `percent_model_SE_**MI-K-logT**`, `Coverage_**MI-K-logT**`,
    `Bias_**MI-K-H0T**`, `percent_model_SE_**MI-K-H0T**`, `Coverage_**MI-K-H0T**`
  ) |>
  dplyr::mutate(Parameters = as.character(Parameters)) |>
  dplyr::arrange(match(Parameters, c("log.lambda", "log.rho", "mgene", "PRS", "log.kappa")))
cat(format(as.vector(unname(unlist(tab_sim[1, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[2, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[3, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[4, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")
cat(format(as.vector(unname(unlist(tab_sim[5, -1]))), trim = TRUE, drop0trailing = TRUE), sep = " & ")

Mis_Gamma_80 <- Mis_Gamma_80 |>
  filter(Parameters == "mgene" | Parameters == "PRS")
## Boxplot
Mis_Gamma_80_HighCorr <- Mis_Gamma_80 |>
  ggplot(aes(x = Parameters, y = Bias2*100, fill = Method)) +
  geom_boxplot(
    aes(group = interaction(Parameters, Method)), 
    position = position_dodge(width = 0.8), 
    size = 0.2,
    width = 0.5
  ) +
  stat_summary(
    aes(group = Method),
    fun = mean, 
    geom = "point", 
    position = position_dodge(width = 0.8), 
    shape = 23,           
    size = 1,           
    color = "black",     
    fill = "white"       
  ) +
  ylab("Bias (%)") +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "black", 
    size = 0.3
  ) +
  labs(fill = "Method") +
  ggtitle("Biases - Gamma Frailty Data, Mis-specified, 80% Missing, Low Correlations") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_markdown(size = 14)
  ) +
  scale_y_continuous(
    limits = c(-140, 140), 
    breaks = seq(-140, 140, by = 20)
  )
show_and_save(Mis_Gamma_80_HighCorr, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Bias Bar Plots Folder/Low Corr/Mis-Specified/80 Missing/Gamma Using LogNormal/Mis_Gamma_80_LowCorr.jpeg")
