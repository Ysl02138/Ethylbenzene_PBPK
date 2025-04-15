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
source("MCSim/function.R")

# posterior check
Mouse_out <- c("outputs/EBMouseMCMC_3365.out",
             "outputs/EBMouseMCMC_4880.out",
             "outputs/EBMouseMCMC_5916.out",
             "outputs/EBMouseMCMC_6734.out")
Mouse_data <- Mouse_out |> map(fread) |> map(as.data.frame)
n_chains <- length(Mouse_data)
sample_number <- dim(Mouse_data[[1]])[1]
dim <- c(sample_number, n_chains, dim(Mouse_data[[1]])[2])
n_iter <- dim(Mouse_data[[1]])[1]
n_param <- dim(Mouse_data[[1]])[2]
Mouse_x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  Mouse_x[, i, ] <- as.matrix(Mouse_data[[i]][1:n_iter, ])
}
dimnames(Mouse_x)[[3]] <- names(Mouse_data[[1]])
dim(Mouse_x)

# Save to RData
Mouse_mcmc <- Mouse_x[seq(100001, 200000, 10), , ]
save(Mouse_mcmc, file = "outputs/EBMouse_mcmc.RData")

# tidy
rm(list = ls())

# data manipulate (a total random sample of 500 iters from 4 chains combined)
load("outputs/EBMouse_mcmc.RData")
no_sample <- 125
set.seed(12345)
sample_iters <- sample(seq_len(dim(Mouse_mcmc)[1]), no_sample)
sample_Mouse_mcmc <- Mouse_mcmc[sample_iters, , ]
nd2 <- dim(sample_Mouse_mcmc)[3]
dim(sample_Mouse_mcmc) <- c(4 * no_sample, nd2)
dim(sample_Mouse_mcmc)

# posterior predictive simulation
model <- "EBPBPK.model"
if (!file.exists("mcsim.EBPBPK.model.exe")) {
  RMCSim::makemcsim(model, dir = "MCSim")
}
for (iter in seq(dim(sample_Mouse_mcmc)[1])){
  head(sample_Mouse_mcmc, iter) |> tail(1) |>
    write.table(file = "MCMCMouse.check.dat", row.names = FALSE, sep = "\t")
  vld <- "./mcsim.EBPBPK.model.exe MCSim/EBMouse.MCMC.check.in"
  system(vld)
  
  Mouse_out <- read.delim("MCMCMouse.check.out")
  Mouse_out$iter <- iter
  if (iter == 1) Mouse_xx <- Mouse_out
  else Mouse_xx <- rbind(Mouse_xx, Mouse_out)
}
Mouse_xx$Output_Var |> unique()


# output manipulate
Mouse_xx <- Mouse_xx |>
  mutate(organs = ifelse(Output_Var == "CalvPPM", "Alveloar",
                       ifelse(Output_Var == "CV", "Blood (venous)",
                              ifelse(Output_Var == "Cart", "Arterial blood",
                                     ifelse(Output_Var == "CFat", "Fat",
                                            ifelse(Output_Var == "CLung", "Lung",
                                                   ifelse(Output_Var == "CLiv", "Liver",
                                                         ifelse(Output_Var == "AUrineMAmg", "Urine MA", "Total metabolites"))))))))
Mouse_xx <- Mouse_xx |>
  mutate(label = ifelse(Simulation == 1, "75 ppm, 4 hrs, female (Tardif, 06)",
                        ifelse(Simulation == 2, "200 ppm, 4 hrs, female (Tardif, 06)",
                               ifelse(Simulation == 3, "500 ppm, 4 hrs, female (Tardif, 06)",
                                      ifelse(Simulation == 4, "1000 ppm, 4 hrs, female (Tardif, 06)",
                                             ifelse(Simulation == 5, "75 ppm, 6 hrs, female (Tardif, 06)",
                                                    ifelse(Simulation == 6, "750 ppm, 6 hrs, female (Tardif, 06)",
                                                           ifelse(Simulation == 7, "75 ppm, 6 hrs, male (Tardif, 06)",
                                                                  ifelse(Simulation == 8, "750 ppm, 6 hrs, male (Tardif, 06)",
                                                                         ifelse(Simulation == 9, "75 ppm, 6 hrs, female (Fuciarelli, 00)",
                                                                                "750 ppm, 6 hrs, male (Fuciarelli, 00)"))))))))))
Mouse_xx$Data[Mouse_xx$Data == -1] <- NA
adj_level <- Mouse_xx$label |> unique()
Mouse_xx$label <- factor(Mouse_xx$label, level = adj_level)
Mouse_xx |> tail()

# define plotting element
set_theme <- theme(
  axis.text.y      = element_text(color = "black", size = 10),
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black", size = 10),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  axis.title       = element_blank(),
  panel.background = element_blank()
)
options(warn=-1)

plmedian <- Mouse_xx |>
  filter(Simulation %in% c(1:4) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p1 <- Mouse_xx |>
  filter(Simulation %in% c(1:4) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-2, 10^2.5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = plmedian, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p2median <- Mouse_xx |>
  filter(Simulation %in% c(5:8) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- Mouse_xx |>
  filter(Simulation %in% c(5:8) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-2, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p3median <- Mouse_xx |>
  filter(Simulation %in% c(9:10) & Time > 0) %>% group_by(Output_Var, Time, organs, label)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- Mouse_xx |>
  filter(Simulation %in% c(9:10) & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-4, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


# add the title and axis label
title <- ggdraw() +
  draw_label(
    "(4A) Mouse",
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
    "EB concentration in tissues (ug/L)",
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
                rel_heights = c(1 / 2, 1/ 2)),
      plot_grid(
        p3, nrow = 1, labels = c("III")
      ),
      nrow = 1, rel_widths = c(2, 2)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/Figure_4A_calibration_Mouse.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()

