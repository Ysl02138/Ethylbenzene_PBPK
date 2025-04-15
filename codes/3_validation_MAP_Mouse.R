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
Mousesim1.1 <- fread("outputs/EBMouseMCMC_3365.out") 
Mousesim2.1 <- fread("outputs/EBMouseMCMC_6734.out") 
Mousesim3.1 <- fread("outputs/EBMouseMCMC_4880.out") 
Mousesim4.1 <- fread("outputs/EBMouseMCMC_5916.out") 

Mouse_x <-mcmc_array(list(Mousesim1.1, Mousesim2.1, Mousesim3.1, Mousesim4.1))
pars_name <- dimnames(Mouse_x)[[3]]
str <- which(pars_name == "Ve_CV(1)")
end <- which(pars_name == "lnKMRptC(1.2)")
parms <- pars_name[str:end]

j <- seq(100001, 200000, 10) 
sum_chains <- length(j)*4

n = 500 # 500 virtual study
d <- Mouse_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Mouse_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Mouse_x <- d[i, parms] 
tmp.Mouse_x %>% write.table(file="EB_validation_MAP_Mouse.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/EB_validation_MAP_Mouse.in"
system(vld)

df <- fread("EB_validation_MAP_Mouse.out") %>% as.data.frame()
file.remove(c("EB_validation_MAP_Mouse.dat", "EB_validation_MAP_Mouse.out"))

# Import Mouse EB validation datasets
Time <- c(2.022 , 3.027 , 4, 4.018 , 4.043 , 4.086 , 4.129 , 4.168 , 4.204 , 
          
          2.022 , 3.032 , 4, 4.026 , 4.176 , 4.348 , 4.663 , 4.999 , 
          
          2.021 , 3.021 , 4, 4.035 , 4.33 , 4.66 , 5.015 , 5.506 , 
          
          2.02 , 3.02 , 4, 4.03 , 4.5 , 5 , 5.5 , 6 , 
          
          3 , 5 , 6, 6.04 , 6.08 , 6.13 , 6.2 , 
          
          3 , 5 , 6 , 6.667 , 7.333 , 8 , 
          
          3 , 5 , 6, 6.04 , 6.08 , 6.13 , 6.2 , 
          
          3 , 5 , 6 , 6.67 , 7.34 , 8 , 

    
          6.268 , 6.269 , 6.27 , 6.368 , 6.369 , 6.37 , 6.484 , 6.485 , 6.586 , 6.601 , 6.602 , 6.734 , 6.736 , 6.752 , 7.051 , 7.053 , 7.069 , 7.518 , 7.519 , 7.52 , 8.018 , 8.019 , 8.02 , 9.036 , 9.068 , 9.085 , 10.02 , 10.034 , 10.069 , 11.034 , 11.036 , 11.052 , 12.053 , 12.068 , 12.069 , 14.051 , 14.052 , 14.053 , 15.084 , 15.086,
        
          6.268 , 6.269 , 6.27 , 6.368 , 6.369 , 6.37 , 6.484 , 6.485 , 6.52 , 6.584 , 6.602 , 6.603 , 6.734 , 6.735 , 6.753 , 7.051 , 7.052 , 7.07 , 7.518 , 7.519 , 7.52 , 8.018 , 8.019 , 8.02 , 9.034 , 9.069 , 9.086 , 10.018 , 10.035 , 10.07 , 11.034 , 11.035 , 12.051 , 12.069 , 12.07 , 14.051 , 15.084 , 15.085,
          
          6.268 , 6.269 , 6.27 , 6.368 , 6.369 , 6.37 , 6.484 , 6.486 , 6.519 , 6.586 , 6.601 , 6.602 , 6.734 , 6.736 , 6.752 , 7.051 , 7.053 , 7.069 , 7.518 , 7.519 , 8.02 , 
          
          6.268 , 6.269 , 6.27 , 6.368 , 6.369 , 6.37 , 6.484 , 6.485 , 6.52 , 6.584 , 6.602 , 6.603 , 6.734 , 6.735 , 6.753 , 7.051 , 7.052 , 7.07 , 7.518 , 7.519 , 8.018,
          
          6.251 , 6.252 , 6.253 , 6.334 , 6.335 , 6.37 , 6.418 , 6.419 , 6.436 , 6.501 , 6.502 , 6.52 , 6.684 , 6.685 , 6.686 , 6.834 , 6.835 , 6.853 , 7.018 , 7.019 , 7.02 , 7.334 , 7.335 , 7.336 , 7.668 , 7.685 , 7.686 , 7.934 , 7.935,
            
          6.251 , 6.252 , 6.253 , 6.334 , 6.37 , 6.418 , 6.419 , 6.501 , 6.502 , 6.52 , 6.684 , 6.685 , 6.834 , 6.835 , 6.853 , 7.018 , 7.019 , 7.02 , 7.334 , 7.935) # Observation time 

Data <- c( 0.606 , 0.623 , -1, 0.711 , 0.566 , 0.283 , 0.156 , 0.128 , 0.088 , 
           
           2.191 , 2.281 , -1, 2.148 , 0.672 , 0.283 , 0.186 , 0.128 , 
           
           15.993 , 19.767 , -1, 10.431 , 2.913 , 0.587 , 1.014 , 0.149 , 
           
           61.689 , 75.813 , -1, 66.705 , 41.15 , 28.356 , 22.328 , 13.742 , 
           
           0.553 , 0.419 , -1, 0.364 , 0.137 , 0.071 , 0.036 , 
           
           37.13 , 39.375 , 32.927 , 13.125 , 10.407 , 2.329 , 
           
           0.57 , 0.72 , -1,  0.54 , 0.26 , 0.12 , 0.06 , 
           
           46 , 34.73 , 27.98 , 5.33 , 0.57 , 0.02 , 
           
           11.1 , 9.37 , 8.05 , 7.58 , 6.41 , 5.64 , 6.65 , 5.1 , 4.67 , 3.43 , 3.37 , 2.6 , 2.67 , 4.34 , 1.55 , 3.25 , 2.59 , 0.433 , 0.429 , 0.671 , 0.152 , 0.0412 , 0.294 , 0.0349 , 0.102 , 0.027 , 0.0145 , 0.0252 , 0.0102 , 0.00302 , 0.00453 , 0.00535 , 0.00179 , 0.0072 , 0.00233 , 0.00162 , 0.00148 , 0.00218 , 0.00188 , 0.00197,
           
           608 , 442 , 233 , 441 , 288 , 238 , 235 , 151 , 502 , 182 , 238 , 168 , 130 , 114 , 111 , 84.6 , 74.2 , 170 , 50.4 , 43.1 , 4.02 , 35.5 , 16.1 , 1.09 , 0.381 , 12.1 , 0.131 , 0.437 , 0.705 , 0.456 , 0.554 , 0.435 , 1.02 , 0.783 , 0.622 , 0.541 , 0.543 , 0.315,
           
           10.8 , 8.9 , 8.54 , 6.76 , 6.16 , 4.45 , 4.94 , 4.48 , 6.81 , 2.84 , 2.38 , 1.53 , 1.32 , 0.871 , 3.79 , 0.529 , 0.431 , 1.53 , 0.0574 , 0.0466 , 0.0559 , 

           7.79 , 1.72 , 0.825 , 12.1 , 2.14 , 0.458 , 5.75 , 1.27 , 1.66 , 0.852 , 0.499 , 0.186 , 0.343 , 0.247 , 0.328 , 0.276 , 0.236 , 0.941 , 0.18 , 0.0214 , 0.265, 
           
           0.0624 , 0.0458 , 0.0261 , 0.0395 , 0.0301 , 0.0218 , 0.0277 , 0.0174 , 0.0206 , 0.0173 , 0.0132 , 0.0119 , 0.0163 , 0.0122 , 0.00872 , 0.0159 , 0.00519 , 0.00701 , 0.00502 , 0.00369 , 0.00346 , 0.00344 , 0.00228 , 0.00212 , 0.002 , 0.00261 , 0.00144 , 0.00335 , 0.00164,
           
           15 , 6.98 , 3.05 , 13.9 , 9.62 , 9.74 , 9.73 , 5.79 , 3.15 , 6.32 , 5.72 , 2.76 , 6.47 , 0.485 , 3.38 , 3.65 , 0.989 , 0.401 , 0.289 , 0.295)



 
formattable(Time , digits = 3, format = "f")
val_data <- data.frame(Data, Time)
study <- c("Tardiff, 2006", "Fuciarelli, 2000")
dose <- c("75 ppm, 4h, male", "200 ppm, 4h, male", "500 ppm, 4h, male","1000 ppm, 4h, male","75 ppm, 6h, male",
"750 ppm, 6h, male", "75 ppm, 6h, female","750 ppm, 6h, female")


organs <- c(rep("Blood (venous)", 99), rep("Fat", 38), rep("Liver", 21),rep("Lung", 21),
            rep("Blood (venous)", 29), rep("Fat", 20))
Output_Var <- c(rep("CV", 99), rep("CFat", 38), rep("CLiv", 21),rep("CLung", 21),
            rep("CV", 29), rep("CFat", 20))
study <- c(rep("Tardiff, 2006", 59), 
           rep("Fuciarelli, 2000", 169))
dose <- c(rep("75 ppm, 4h, male", 9),
          rep("200 ppm, 4h, male", 8),
          rep("500 ppm, 4h, male", 8),
          rep("1000 ppm, 4h, male", 8),
          rep("75 ppm, 6h, male", 7),
          rep("750 ppm, 6h, male", 6),
          rep("75 ppm, 6h, female", 7),
          rep("750 ppm, 6h, female", 6),
          rep("750 ppm, 6h, female", 120),
          rep("75 ppm, 6h, male", 49))
Simulation <- c(
          rep("1", 9),
          rep("2", 8),
          rep("3", 8),
          rep("4", 8),
          rep("5", 7),
          rep("6", 6),
          rep("7", 7),
          rep("8", 6),
          rep("9", 120),
          rep("10", 49))

val_data$organs <- organs
val_data$study  <- study
val_data$dose <- dose
val_data$Simulation <- Simulation
val_data  %<>% unite(label, dose, study, sep = " ", remove = FALSE)
label <- c(val_data$label)
val_data <- select(val_data, -dose, -study)

str_1.1 <- which(names(df) == "CV_1.1")
end_10.20 <- which(names(df) == "CFat_10.20")

qt.line <- df[,c(str_1.1:end_10.20)] %>% gather() %>% 
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
options(warn=-1)

p1median <- qt.line  |>
  filter(Simulation %in% c(1:4) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p1 <- qt.line  |>
  filter(Simulation %in% c(1:4) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-4.5, 10^3.5),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p1median, aes(Time, median), color = "black") +
  geom_line(data = p1median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p1median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p2median <- qt.line  |>
  filter(Simulation %in% c(5:8) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- qt.line  |>
  filter(Simulation %in% c(5:8) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-4.5, 10^3.5),
    breaks = trans_breaks("log10", function(x) 10^x, n = 4),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p2median, aes(Time, median), color = "black") +
  geom_line(data = p2median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p2median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme


p3median <- qt.line  |>
  filter(Simulation %in% c(9:10) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- qt.line  |>
  filter(Simulation %in% c(9:10) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-15, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_line(data = p3median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p3median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free", ncol=2) +
  theme_bw() +
  set_theme

# add the title and axis label
title <- ggdraw() +
  draw_label(
    #"(5A)"
    "Mouse",
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
    "EB concentration in tissues (mg/L)",
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
ggsave(file = "plots/Figure_5A_MAP_Validation_Mouse.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()

