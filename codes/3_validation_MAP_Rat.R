# 
#   This code will produce Figure 5B for rat validation dataset to evaluate model performance after MCMC calibration
#

# Load packages ----------------------------------------------------------------

library(rstan)
library(tidyverse)
library(bayesplot)
library(GGally)
library(ggpubr)
library(scales)
library(cowplot)
library(magrittr)
library(LaplacesDemon)
library(bayestestR)
library(PowerTOST)
library(foreach)
library(doParallel) 
library(PKNCA)
library(kableExtra)
library(data.table)
library(formattable)
source("MCSim/function.R")


# Validation -------------------------------------------------------------------
Ratsim1.1 <- fread("outputs/EBRatMCMC_3365.out") 
Ratsim2.1 <- fread("outputs/EBRatMCMC_6734.out") 
Ratsim3.1 <- fread("outputs/EBRatMCMC_4880.out") 
Ratsim4.1 <- fread("outputs/EBRatMCMC_5916.out") 

Rat_x <-mcmc_array(list(Ratsim1.1, Ratsim2.1, Ratsim3.1, Ratsim4.1))
pars_name <- dimnames(Rat_x)[[3]]
str <- which(pars_name == "Ve_CV(1)")
end <- which(pars_name == "lnKpgaC(1.3)") 
parms <- pars_name[str:end]

j <- seq(100001, 200000, 10) 
sum_chains <- length(j)*4 #

n = 500 # 500 virtual study
d <- Rat_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Rat_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Rat_x <- d[i, parms] 
tmp.Rat_x %>% write.table(file="EB_validation_MAP_Rat.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/EB_validation_MAP_Rat.in"
system(vld)

df <- fread("EB_validation_MAP_Rat.out") %>% as.data.frame()
file.remove(c("EB_validation_MAP_Rat.dat", "EB_validation_MAP_Rat.out"))

# Import Rat EB validation dataset --- same datasets as "EB_validation_MAP_Rat.in"
Time <- c(6.384 , 6.402 , 6.42 , 6.584 , 6.602 , 6.636 , 6.784 , 6.819 , 6.886 , 7.151 , 7.169 , 7.186 , 8.034 , 8.085 , 8.103 , 9.034 , 9.069 , 9.086 , 12.184 , 12.185 , 12.236 , 15.201 , 15.269 , 15.336 , 18.118 , 18.202 , 18.27 , 24.084 , 24.119 , 24.136 , 29.818 , 29.851 , 29.852 , 
          
          6.383 , 6.4 , 6.417 , 6.583 , 6.6 , 6.633 , 6.783 , 6.817 , 6.883 , 7.15 , 7.167 , 7.183 , 8.033 , 
          
          6.384 , 6.402 , 6.42 , 6.584 , 6.602 , 6.636 , 6.784 , 6.819 , 6.886 , 7.151 , 7.169 , 7.186 , 8.034 , 8.085 , 8.103 , 9.034 , 9.069 , 9.086 , 12.184 , 12.185 , 
          12.236 , 15.201 , 15.269 , 15.336 , 18.118 , 18.202 , 18.27 , 24.119 , 24.136 , 29.818 , 
          
          6.384 , 6.402 , 6.42 , 6.584 , 6.602 , 6.636 , 6.784 , 6.819 , 6.886 , 7.151 , 7.169 , 7.186 , 8.034 , 8.085 , 8.103 , 9.034 , 9.069 , 9.086 , 12.184 , 12.236 , 15.201 , 15.269 , 15.336 , 18.118 , 18.202 , 18.27 , 24.084 , 
          
          12.001 , 12.002 , 12.003 , 12.004 , 12.005 , 24.001 , 24.002 , 24.003 , 24.004 , 24.005 , 48.001 , 48.002 , 48.003 , 48.004 , 48.005 , 
          
          6.268 , 6.319 , 6.353 , 6.518 , 6.535 , 6.553 , 6.968 , 7.002 , 7.003 , 7.484 , 7.485 , 7.486 , 8.768 , 8.869 , 9.003 , 11.984 , 12.002 , 12.003 , 14.918 , 14.969 , 15.036 , 17.918 , 17.969 , 18.036 , 23.951 , 23.952 , 23.953 , 29.934 , 29.952 , 29.97 , 41.951 , 42.119 , 42.236 , 
          
          6.268 , 6.319 , 6.353 , 6.518 , 6.535 , 6.553 , 6.968 , 7.002 , 7.003 , 7.484 , 7.485 , 7.486 , 8.768 , 8.869 , 9.003 , 11.984 , 12.002 , 12.003 , 14.918 , 14.969 , 15.036 , 17.918 , 17.969 , 
          
          6.268 , 6.319 , 6.353 , 6.518 , 6.535 , 6.553 , 6.968 , 7.002 , 7.003 , 7.484 , 7.485 , 7.486 , 8.768 , 8.869 , 9.003 , 11.984 , 12.002 , 12.003 , 14.918 , 14.969 , 15.036 , 17.918 , 17.969 , 18.036 , 23.951 , 23.952 , 23.953 , 29.934 , 29.952 , 29.97 , 41.951 , 42.119 , 42.236 , 
          
          6.268 , 6.319 , 6.353 , 6.518 , 6.535 , 6.553 , 6.968 , 7.002 , 7.003 , 7.484 , 7.485 , 7.486 , 8.768 , 8.869 , 9.003 , 11.984 , 12.002 , 12.003 , 14.918 , 14.969 , 15.036 , 17.918 , 17.969 , 18.036 , 23.951 , 23.952 , 23.953 , 29.934 , 29.952 , 29.97 , 42.119 , 42.236 , 
          
          12.001 , 12.002 , 12.003 , 12.005 , 24.001 , 24.002 , 24.003 , 24.004 , 24.005 , 48.001 , 48.002 , 48.003 , 48.004 , 48.005 , 
          
          6 , 24,  48 , 6 , 24,  48 ,
          8.05 , 8.2) # Observation time 

Data <- c(0.145 , 0.174 , 0.133 , 0.0892 , 0.141 , 0.118 , 0.114 , 0.0694 , 0.0911 , 0.0735 , 0.0768 , 0.107 , 0.0446 , 0.0693 , 0.053 , 0.0462 , 0.0426 , 0.061 , 0.0188 , 0.0165 , 0.0162 , 0.0112 , 0.014 , 0.0194 , 0.00883 , 0.00949 , 0.00763 , 0.00309 , 0.00435 , 0.00408 , 0.00246 , 0.0036 , 0.00277, 
          
          0.0749 , 0.274 , 0.0616 , 0.0192 , 0.0476 , 0.0638 , 0.0509 , 0.0176 , 0.0296 , 0.04 , 0.0237 , 0.0507 , 0.0611 , 
          
          9.37 , 3.44 , 3.33 , 14.1 , 4.58 , 6.37 , 22.3 , 4.28 , 5.66 , 18.6 , 1.91 , 5.91 , 9.27 , 2.3 , 2.71 , 1.98 , 0.882 , 3.69 , 0.886 , 0.636 , 0.369 , 0.937 , 0.936 , 2.16 , 1.73 , 0.772 , 0.175 , 0.256 , 0.239 , 1.55 , 
          
          1.1 , 0.387 , 0.25 , 0.279 , 0.323 , 0.141 , 0.545 , 0.108 , 0.286 , 0.343 , 0.257 , 0.114 , 0.05 , 0.0261 , 0.0535 , 0.0488 , 0.0285 , 0.0963 , 0.0213 , 0.015 , 0.0207 , 0.0308 , 0.0364 , 0.0179 , 0.015 , 0.0164 , 0.0126 , 
          
          1.29 , 1.11 , 1.05 , 1.04 , 0.925 , 1.395 , 1.216 , 1.18 , 1.199 , 0.9898 , 1.474 , 1.397 , 1.2446 , 1.36 , 1.0724 , 
          
          
          6.98 , 5.05 , 4.98 , 4.83 , 4.85 , 4.66 , 3.68 , 3.43 , 3.26 , 2.93 , 2.82 , 2.6 , 1.22 , 0.952 , 1.28 , 0.371 , 0.482 , 0.428 , 0.214 , 0.185 , 0.185 , 0.141 , 0.19 , 0.0855 , 0.0807 , 0.0688 , 0.0553 , 0.0543 , 0.0504 , 0.0262 , 0.0128 , 0.0106 , 0.00906 , 
          
          8.61 , 8.44 , 6.7 , 6.67 , 7.64 , 8.21 , 6.73 , 6.56 , 5.1 , 5.12 , 4.38 , 3.11 , 1.76 , 1.59 , 1.44 , 0.296 , 0.421 , 0.276 , 0.079 , 0.0566 , 0.0344 , 0.0638 , 0.0948 , 
          
          89.4 , 178 , 285 , 315 , 112 , 293 , 199 , 214 , 164 , 173 , 121 , 113 , 97.9 , 86.4 , 214 , 30 , 44.8 , 18.8 , 13.9 , 30.7 , 6.82 , 38 , 8.75 , 11.6 , 9.8 , 2.13 , 1.61 , 0.866 , 2.7 , 0.89 , 0.0974 , 0.183 , 0.144 , 
          
          7.06 , 4.13 , 5.47 , 3.69 , 3.8 , 3.84 , 13 , 2.75 , 2.46 , 1.96 , 1.76 , 1.66 , 0.689 , 0.984 , 0.824 , 0.417 , 0.474 , 0.321 , 0.278 , 0.172 , 0.361 , 0.247 , 0.0963 , 0.128 , 0.0881 , 0.0739 , 0.0547 , 0.0322 , 0.0923 , 0.084 , 0.0216 , 0.0259 , 
          
          5.29 , 3.44 , 2.19 , 1.52 , 7.06 , 7.79 , 5.88 , 6.06 , 4.31 , 7.856 , 8.83 , 6.712 , 9.1 , 4.84, 
          
          0.288, 7.093, 9.266, 0.255612, 9.649353, 16.550877 ,  
          23.23 , 5.7) # Observation concentration

formattable(Time , digits = 3, format = "f")
val_data <- data.frame(Data, Time)
study <- c("Fuciarelli, 2000", "Engstrom, 1984", "Cappaert, 2002")
dose <- c("75 ppm, 6h, female", "750 ppm, 6h, male", "300 ppm, 6h, male","600 ppm, 6h, male", "500 ppm, 8h, male")

organs <- c(rep("Blood (venous)", 33), rep("Liver", 13), rep("Fat", 30), rep("Lung", 27), rep("Urine MA", 15),
            rep("Blood (venous)", 33), rep("Liver", 23), rep("Fat", 33), rep("Lung", 32), rep("Urine MA", 14),
            rep("Urine MA", 6), rep("Blood (venous)", 2))

Output_Var <- c(rep("CV", 33), rep("CLiv", 13), rep("CFat", 30), rep("CLung", 27), rep("AUrineMAmg", 15),
                rep("CV", 33), rep("CLiv", 23), rep("CFat", 33), rep("CLung", 32), rep("AUrineMAmg", 14),
                rep("AUrineMAmg", ,6), rep("CV", 2))

study <- c(rep("Fuciarelli, 2000", 253),
           rep("Engstrom, 1984", 6), 
           rep("Cappaert, 2002", 2))

dose <- c(rep("75 ppm, 6h, female", 118),
          rep("750 ppm, 6h, male", 135),
          rep("300 ppm, 6h, male", 3),
          rep("600 ppm, 6h, male", 3),
          rep("500 ppm, 8h, male", 2))

Simulation <- c(
  rep("1", 118),
  rep("2", 135),
  rep("3", 3),
  rep("4", 3),
  rep("5", 2))


val_data$organs <- organs
val_data$study  <- study
val_data$dose <- dose
val_data$Simulation <- Simulation
val_data  %<>% unite(label, dose, study, sep = " ", remove = FALSE)
label <- c(val_data$label)
val_data <- select(val_data, -dose, -study)

str_1.1 <- which(names(df) == "CV_1.1")
end_5.2 <- which(names(df) == "CV_5.2")

qt.line <- df[,c(str_1.1:end_5.2)] %>% gather() %>% 
  `colnames<-`(c("var","Prediction")) %>%
  add_column(Time =  rep(rep(Time, each = n),1)) %>%
  add_column(Data =  rep(rep(Data, each = n),1)) %>%
  add_column(label = rep(rep(label, each = n),1)) %>%
  add_column(Simulation =  rep(rep(Simulation, each = n),1)) %>%
  add_column(Output_Var =  rep(rep(Output_Var, each = n),1)) %>%
  add_column(organs = rep(rep(organs, each = n),1))

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
  panel.background = element_blank())




p1median <- qt.line  |>
  filter(Simulation == 1 & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p1 <- qt.line  |>
  filter(Simulation == 1 & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-15, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p1median, aes(Time, median), color = "black") +
  geom_line(data = p1median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p1median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol=1) +
  theme_bw() +
  set_theme


p2median <- qt.line  |>
  filter(Simulation == 2 & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- qt.line  |>
  filter(Simulation == 2 & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-15, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_line(data = p2median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p2median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol=1) +
  theme_bw() +
  set_theme


p3median <- qt.line  |>
  filter(Simulation %in% c(3:4) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- qt.line  |>
  filter(Simulation %in% c(3:4) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-1, 10^2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_line(data = p3median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p3median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol=2) +
  theme_bw() +
  set_theme

p4median <- qt.line  |>
  filter(Simulation %in% c(5) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p4 <- qt.line  |>
  filter(Simulation %in% c(5) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^0, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p4median, aes(Time, median), color = "black") +
  geom_line(data = p4median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p4median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data, group = 1)) +
  facet_wrap(organs ~ label, scales = "free", ncol=2) +
  theme_bw() +
  set_theme


# add the title and axis label
title <- ggdraw() +
  draw_label(
    #(5B) 
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
      plot_grid(p3, p4, nrow = 2, labels = c("I", "II"),
                rel_heights = c(2 / 4, 2/ 4)),
      plot_grid(
        p1,p2, nrow = 1, labels = c("III","IV")
      ),
      nrow = 1, rel_widths = c(2, 3)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/Figure_5B_MAP_Validation_Rat.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()