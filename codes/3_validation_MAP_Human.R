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
Humansim1.1 <- fread("outputs/EBHumanMCMC_3365.out") 
Humansim2.1 <- fread("outputs/EBHumanMCMC_6734.out") 
Humansim3.1 <- fread("outputs/EBHumanMCMC_4880.out") 
Humansim4.1 <- fread("outputs/EBHumanMCMC_5916.out") 

Human_x <-mcmc_array(list(Humansim1.1, Humansim2.1, Humansim3.1, Humansim4.1))
pars_name <- dimnames(Human_x)[[3]]
str <- which(pars_name == "Ve_CalvPPM(1)")
end <- which(pars_name == "lnKpgaC(1.3)") 
parms <- pars_name[str:end]

j <- seq(100001, 200000, 10) 
sum_chains <- length(j)*4

n = 500 # 500 virtual study
d <- Human_x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(Human_x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.Human_x <- d[i, parms] 
tmp.Human_x %>% write.table(file="EB_validation_MAP_Human.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/EB_validation_MAP_Human.in"
system(vld)

df <- fread("EB_validation_MAP_Human.out") %>% as.data.frame()
file.remove(c("EB_validation_MAP_Human.dat", "EB_validation_MAP_Human.out"))

# Import Human EB validation datasets
Time <- c(1 , 4 , 8.75 , 1.776649746 , 2.131979695 , 3.819796954 , 3.829796954 , 3.839796954 , 4.61928934 , 4.974619289 , 5.152284264 , 6.395939086 , 7.017766497 , 7.817258883 , 8.350253807 , 8.439086294 , 8.972081218 , 8.973081218 , 9.682741117 , 10.4822335 , 10.5822335 , 11.37055838 , 11.63705584 , 11.9035533 , 12.6142132 , 13.0887176 , 13.76903553 , 14.39086294 , 14.47969543 , 15.19035533 , 15.27918782 , 15.45685279 , 17.67766497 , 19.18781726 , 21.23096447 , 21.40862944 , 21.94162437 , 22.47461929 , 23.54060914 ,4	, 8, 	24, 48) # Observation concentration
Data <- c(0.105, 0.184, 0.186, 10.7257732 , 2.383505155 , 21.84879725 , 24.62955326 , 28.60206186 , 32.97182131 , 28.60206186 , 69.51890034 , 42.10859107 , 48.06735395 , 58.39587629 , 65.94364261 , 63.56013746 , 70.31340206 , 60.38213058 , 84.21718213 , 80.24467354 , 78.25841924 , 80.24467354 , 77.46391753 , 86.99793814 , 92.16219931 , 91.36769759 , 112.819244 , 88.18969072 , 101.2989691 , 94.54570447 , 96.92920962 , 99.31271478 , 112.819244 , 114.444 , 106.8604811 , 109.2439863 , 110.038488 , 103.2852234 , 112.0247423, 132.1460483,	421.8508466,	635.31754, 0.0000417) # Observation time 







formattable(Time , digits = 3, format = "f")
val_data <- data.frame(Data, Time)
study <- c("Knecht, 2000", "Engstrom, 1984","Lin, 2017" )
dose <- c("25 ppm, 8h, male", "150 ppm, 4h, male", "2 ppb, population")

organs <- c(rep("Blood (venous)", 3),rep("Urine MA", 39),rep("Blood (venous)", 1))
Output_Var <- c(rep("CV", 3), rep("AUrineMAmg", 39), rep("CV", 1))
study <- c(rep("Knecht, 2000", 3),
           rep("Knecht, 2000", 36), 
           rep("Engstrom, 1984", 3),
           rep("Lin, 2017", 1))
dose <- c(rep("25 ppm, 8h, male", 3),
          rep("25 ppm, 8h, male", 36),
          rep("150 ppm, 4h, male", 3),
          rep("2 ppb, population", 1))

Simulation <- c(
  rep("1", 3),
  rep("2", 36),
  rep("3", 3),
  rep("4", 1))

val_data$organs <- organs
val_data$study  <- study
val_data$dose <- dose
val_data$Simulation <- Simulation
val_data  %<>% unite(label, dose, study, sep = " ", remove = FALSE)
label <- c(val_data$label)
val_data <- select(val_data, -dose, -study)

str_1.1 <- which(names(df) == "CV_1.1")
end_3.1 <- which(names(df) == "CV_3.1")
qt.line <- df[,c(str_1.1:end_3.1)] %>% gather() %>% 
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
  filter(Simulation %in% c(1) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p1 <- qt.line  |>
  filter(Simulation %in% c(1) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-2, 10^0), 
    breaks = trans_breaks("log10", function(x) 10^x, n = 2),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p1median, aes(Time, median), color = "black") +
  geom_line(data = p1median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p1median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme

p2median <- qt.line  |>
  filter(Simulation %in% c(2) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p2 <- qt.line  |>
  filter(Simulation %in% c(2) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^-1.5, 10^2.5), 
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
  filter(Simulation %in% c(3) & Time > 0) %>% group_by(label, organs, Time, Simulation)  %>% 
  summarise(Data = median(Data), median = median(Prediction),
            UCL = quantile(Prediction, 0.975), LCL = quantile(Prediction, 0.025))  
p3 <- qt.line  |>
  filter(Simulation %in% c(3) & Time > 0) |>
  ggplot() +  
  scale_y_log10(lim = c(10^1, 10^3),
    breaks = trans_breaks("log10", function(x) 10^x, n = 3),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(data = p3median, aes(Time, median), color = "black") +
  geom_line(data = p3median, aes(Time, UCL), color = "blue",lty=2) + 
  geom_line(data = p3median, aes(Time, LCL), color = "blue",lty=2) +
  geom_point(aes(x = Time, y = Data)) +
  facet_wrap(organs ~ label, scales = "free") +
  theme_bw() +
  set_theme

## Histogram plot for NHANES
NHANESPred <- qt.line  |> filter(Simulation %in% c(4) & Time == 48) %>% group_by(label, organs, Time, Simulation)
Keeps <- c("Prediction", "organs", "label")
NHANESPred2 <-NHANESPred[ , (names(NHANESPred) %in% Keeps)]# remove average EB and replace with individual data
names(NHANESPred2)[1] <- "Blood_EB"
NHANESPred2$Type <- "Prediction"

NHANESData<-read.csv("outputs/HumanNHANESData.csv", header = TRUE)
names(NHANESData)[1] <- "Blood_EB"
NHANESData$Type <- "Data"
NHANESData$organs <- "Blood (venous)"
NHANESData$label <- "2 ppb, population Lin, 2017"
NHANES <- merge(NHANESPred2, NHANESData, all = TRUE) 
NHANES$Type <- as.factor(NHANES$Type)

p4 <- ggplot(NHANES, aes(Type, Blood_EB)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + 
  scale_y_log10(lim = c(10^-6, 10^-2),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                label = trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(organs ~ label, scales = "free") +
  geom_point(aes(color = Type), position = position_jitter(seed = 1, width = 0.075)) +
  scale_color_manual(values = c("Data" = "black", "Prediction" = "blue")) +
  theme_bw() +
  set_theme

# add the title and axis label
title <- ggdraw() +
  draw_label(
    #(5C) 
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
      plot_grid(p1,p2, nrow = 2, labels = c("I"),
                rel_heights = c(2 / 4, 2/ 4)),
      plot_grid(
        p4, p3, nrow = 2, labels = c("II", "III")
      ),
      nrow = 1, rel_widths = c(2/4, 2/ 4)
    ),
    xlab, nrow = 3, rel_heights = c(0.025, 1, 0.025)),
  nrow = 1, rel_widths = c(0.02, 1)
)
ggsave(file = "plots/Figure_5C_MAP_Validation_Human.eps", device = cairo_ps, height = 12, width = 20, dpi = 600)
dev.off()

