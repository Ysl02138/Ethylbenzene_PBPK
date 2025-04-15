# 
#   This code will produce Figure 4C for human calibration dataset
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

Human_out <- c("outputs/EBHumanMCMC_3365.out",
               "outputs/EBHumanMCMC_4880.out",
               "outputs/EBHumanMCMC_5916.out",
               "outputs/EBHumanMCMC_6734.out")
Human_data <- Human_out |> map(fread) |> map(as.data.frame)
n_chains <- length(Human_data)
sample_number <- dim(Human_data[[1]])[1]
dim <- c(sample_number, n_chains, dim(Human_data[[1]])[2])
n_iter <- dim(Human_data[[1]])[1]
n_param <- dim(Human_data[[1]])[2]
Human_x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  Human_x[, i, ] <- as.matrix(Human_data[[i]][1:n_iter, ])
}
dimnames(Human_x)[[3]] <- names(Human_data[[1]])
dim(Human_x)

# Save to RData
Human_mcmc <- Human_x[seq(100001, 200000, 10), , ]
save(Human_mcmc, file = "outputs/EBHuman_mcmc.RData")

# tidy
rm(list = ls())

# data manipulate (a total random sample of 500 iters from 4 chains combined)
load("outputs/EBHuman_mcmc.RData")
no_sample <- 125
set.seed(12345)
sample_iters <- sample(seq_len(dim(Human_mcmc)[1]), no_sample)
sample_Human_mcmc <- Human_mcmc[sample_iters, , ]
nd2 <- dim(sample_Human_mcmc)[3]
dim(sample_Human_mcmc) <- c(4 * no_sample, nd2)
dim(sample_Human_mcmc)

# posterior predictive simulation


model <- "EBPBPK.model"
if (!file.exists("mcsim.EBPBPK.model.exe")) {
  RMCSim::makemcsim(model, dir = "MCSim")
}
for (iter in seq(dim(sample_Human_mcmc)[1])){
  head(sample_Human_mcmc, iter) |> tail(1) |>
  write.table(file = "MCMCHuman.check.dat", row.names = FALSE, sep = "\t")
  vld <- "./mcsim.EBPBPK.model.exe MCSim/EBHuman.MCMC.check.in"
  # This input file ".in" is the same as MCMC input files to use MCMC outputs for posterior predictive sampling
  system(vld)
  Human_out <- read.delim("MCMCHuman.check.out")
  Human_out$iter <- iter
  if (iter == 1) Human_xx <- Human_out
  else Human_xx <- rbind(Human_xx, Human_out)
}
Human_xx$Output_Var |> unique()

# output manipulate
Human_xx <- Human_xx |>
  mutate(organs = ifelse(Output_Var == "CalvPPM", "Alveloar",
                       ifelse(Output_Var == "CV", "Blood (venous)",
                              ifelse(Output_Var == "Cart", "Arterial blood",
                                     ifelse(Output_Var == "CFat", "Fat",
                                            ifelse(Output_Var == "CLung", "Lung",
                                                   ifelse(Output_Var == "CLiv", "Liver",
                                                         ifelse(Output_Var == "AUrineMAmg", "Urine MA", "Total metabolites"))))))))
Human_xx <- Human_xx |>
  mutate(label = ifelse(Simulation == 1, "33 ppm, 7 hrs (Tardiff, 1997)",
                        ifelse(Simulation == 2, "12.5 ppm, 6 hrs (Marchand, 2015)",
                               ifelse(Simulation == 3, "25 ppm, 6 hrs (Marchand, 2015)",
                                                    "100 ppm, 8 hrs (Knecht, 2000)"))))

Human_xx$Data[Human_xx$Data == -1] <- NA
adj_level <- Human_xx$label |> unique()
Human_xx$label <- factor(Human_xx$label, level = adj_level)
Human_xx |> tail()

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

plmedian <- Human_xx |>
  filter(Simulation == 1 & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p1 <- Human_xx |>
  filter(Simulation == 1 & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-1, 10^1),
                breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = plmedian, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p2median <- Human_xx |>
  filter(Simulation %in% c(2,3) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- Human_xx |>
  filter(Simulation %in% c(2,3) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-2, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol = 2) +
  theme_bw() +
  set_theme


p3median <- Human_xx |>
  filter(Simulation %in% c(4) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- Human_xx |>
  filter(Simulation %in% c(4) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-1, 10^3),
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
    #"(4C)",
    "Human",
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
    "EB concentration in blood or alveolar air (mg/L) / MA excreted, urine (mg)",
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
      plot_grid(p1, p3, nrow = 2, labels = c("I", "II"),
                rel_heights = c(2 / 4, 2/ 4)),
      plot_grid(
        p2, nrow = 1, labels = c("III")
      ),
      nrow = 1, rel_widths = c(2, 2)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/Figure_4C_calibration_Human.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()

