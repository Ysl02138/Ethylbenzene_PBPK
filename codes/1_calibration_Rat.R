# 
#   This code will produce Figure 4B for rat calibration dataset
#

# Set working directory (need to be customized) -------------------------------------
if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("plots")) dir.create("plots")
if(!dir.exists("plots/suppl")) dir.create("plots/suppl")

# set file path ----------------------------------------------
source("MCSim/function.R")
set_PATH()
##############################################################

# package
library(bayesplot)
library(bayestestR)
library(coda)
library(cowplot)
library(data.table)
library(doParallel)
library(dplyr)
library(foreach)
library(GGally)
library(ggplot2)
library(ggpubr)
library(kableExtra)
library(LaplacesDemon)
library(magrittr)
library(PKNCA)
library(pksensi)
library(PowerTOST)
library(purrr)
library(rstan)
library(scales)
library(sensitivity)
library(tidyverse)
library(ggpattern)
source("MCSim/function.R")

# Posterior check 
# The file with suffixes “3365”, “4880”, “5916”, and “6734” referring to the seed numbers used in the MCMC analysis
# The readers can refer to the input files used for the MCMC analysis located in  the folder of "MCMC replication analysis"

Rat_out <- c("outputs/EBRatMCMC_3365.out",
               "outputs/EBRatMCMC_4880.out",
               "outputs/EBRatMCMC_5916.out",
               "outputs/EBRatMCMC_6734.out")
Rat_data <- Rat_out |> map(fread) |> map(as.data.frame)
n_chains <- length(Rat_data)
sample_number <- dim(Rat_data[[1]])[1]
dim <- c(sample_number, n_chains, dim(Rat_data[[1]])[2])
n_iter <- dim(Rat_data[[1]])[1]
n_param <- dim(Rat_data[[1]])[2]
Rat_x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  Rat_x[, i, ] <- as.matrix(Rat_data[[i]][1:n_iter, ])
}
dimnames(Rat_x)[[3]] <- names(Rat_data[[1]])
dim(Rat_x)

# Save to RData
Rat_mcmc <- Rat_x[seq(100001, 200000, 10), , ]
save(Rat_mcmc, file = "outputs/EBRat_mcmc.RData")

# tidy
rm(list = ls())

# data manipulate (random sample 125 iterations from 4 chains)
load("outputs/EBRat_mcmc.RData")
no_sample <- 125
set.seed(12345)
sample_iters <- sample(seq_len(dim(Rat_mcmc)[1]), no_sample)
sample_Rat_mcmc <- Rat_mcmc[sample_iters, , ]
nd2 <- dim(sample_Rat_mcmc)[3]
dim(sample_Rat_mcmc) <- c(4 * no_sample, nd2)
dim(sample_Rat_mcmc)

# posterior predictive simulation
model <- "EBPBPK.model"
if (!file.exists("mcsim.EBPBPK.model.exe")) {
  RMCSim::makemcsim(model, dir = "MCSim")
}
for (iter in seq(dim(sample_Rat_mcmc)[1])){
  head(sample_Rat_mcmc, iter) |> tail(1) |>
  write.table(file = "MCMCRat.check.dat", row.names = FALSE, sep = "\t")
  vld <- "./mcsim.EBPBPK.model.exe MCSim/EBRat.MCMC.check.in" 
  # This input file ".in" is the same as MCMC input files to use MCMC outputs for posterior predictive sampling
  system(vld)
  Rat_out <- read.delim("MCMCRat.check.out")
  Rat_out$iter <- iter
  if (iter == 1) Rat_xx <- Rat_out
  else Rat_xx <- rbind(Rat_xx, Rat_out)
}
Rat_xx$Output_Var |> unique()

# output manipulate
Rat_xx <- Rat_xx |>
  mutate(organs = ifelse(Output_Var == "CalvPPM", "Alveloar",
                         ifelse(Output_Var == "CV", "Blood (venous)",
                                ifelse(Output_Var == "Cart", "Arterial blood",
                                       ifelse(Output_Var == "CFat", "Fat",
                                              ifelse(Output_Var == "CLung", "Lung",
                                                     ifelse(Output_Var == "CLiv", "Liver",
                                                            ifelse(Output_Var == "AUrineMAmg", "Urine MA", "Total metabolites"))))))))
Rat_xx <- Rat_xx |>
  mutate(label = ifelse(Simulation == 1, "50 ppm, 4 hrs, male (Haddad, 99)",
                        ifelse(Simulation == 2, "100 ppm, 4 hrs, male (Haddad, 99)",
                               ifelse(Simulation == 3, "200 ppm, 4 hrs, male (Haddad, 99)",
                                      ifelse(Simulation == 4, "500 ppm, 4 hrs, male (Haddad, 99)",
                                             ifelse(Simulation == 5, "25 ppm, 6 hrs, female (Take, 2020)",
                                                    ifelse(Simulation == 6, "50 ppm, 6 hrs, male (Take, 2020)",
                                                           ifelse(Simulation == 7, "100 ppm, 6 hrs, male (Take, 2020)",
                                                                  ifelse(Simulation == 8, "200 ppm, 6 hrs, male (Take, 2020)",
                                                                         ifelse(Simulation == 9, "75 ppm, 6 hrs, male (Fuciarelli, 2000)",
                                                                                "750 ppm, 6 hrs, female (Fuciarelli, 2000)"))))))))))

Rat_xx$Data[Rat_xx$Data == -1] <- NA
adj_level <- Rat_xx$label |> unique()
Rat_xx$label <- factor(Rat_xx$label, level = adj_level)
Rat_xx |> tail()

# define plotting element
set_theme <- theme(
  axis.text.y      = element_text(color = "black"),
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  axis.title       = element_blank(),
  panel.background = element_blank()
)
options(warn=-1)

p1median <- Rat_xx |>
  filter(Simulation %in% c(1:4) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025)) 
p1 <- Rat_xx |>
  filter(Simulation %in% c(1:4) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-2, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p1median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p2median <- Rat_xx |>
  filter(Simulation %in% c(5:8) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- Rat_xx |>
  filter(Simulation %in% c(5:8) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-2, 10^1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p3median <- Rat_xx |>
  filter(Simulation %in% c(9:10) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- Rat_xx |>
  filter(Simulation %in% c(9:10) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-3.5, 10^2.5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol = 2) +
  theme_bw() +
  set_theme


# add the title and axis label
title <- ggdraw() +
  draw_label(
    # "(4B)"
    "Rat",
    fontface = "bold",
    x = 0,
    size = 18,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 1)
  )
xlab <- ggdraw() +
  draw_label(
    "Time (hr)",
    fontface = "bold", size = 15, hjust = 0,
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )
ylab <- ggdraw() +
  draw_label(
    "EB concentration in tissues (mg/L) / MA excreted, urine (mg)",
    fontface = "bold", size = 15, vjust = 0, angle = 90
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )

# plot
plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(
      plot_grid(p1, p2, nrow = 2, labels = c("I", "II"),
                rel_heights = c(2 / 4, 2/ 4)),
      plot_grid(
        p3, nrow = 1, labels = c("III")
      ),
      nrow = 1, rel_widths = c(2, 2)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/Figure_4B_calibration_Rat.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()
