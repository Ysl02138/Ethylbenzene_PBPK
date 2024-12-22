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

###########################################################################
# Mouse -------------------------------------------------------------------
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
tmp.Mouse_x %>% write.table(file="Mouse.HEC.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/Mouse_HEC.in"
system(vld)

Mouse_out <- read.delim("Mouse.HEC.out")

# Tidy data with median and 95% confidence interval
vars <- names(Mouse_out)
index <- which(vars == "AUCArt_1.1" | vars == "DailyAMetTotalBW34_1.1")
Mouse_Summary <- apply(Mouse_out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(Mouse_Summary) <- c("median", "LCL", "UCL")
df <- as.data.frame(Mouse_Summary)


###########################################################################
# Human -------------------------------------------------------------------
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
tmp.Human_x %>% write.table(file="Human.HEC.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/Human_HEC.in"
system(vld)

Human_out <- read.delim("Human.HEC.out")

# Tidy data with median and 95% confidence interval
vars <- names(Human_out)
index <- which(vars == "AUCArt_1.1" | vars == "DailyAMetTotalBW34_1.1")
Human_Summary <- apply(Human_out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(Human_Summary) <- c("median", "LCL", "UCL")
df <- as.data.frame(Human_Summary)


#####################################################################
## Calculated the HEC 

## Mouse 
AUCArt75       = ((Mouse_out$AUCArt_1.1)/(Human_out$AUCArt_1.1))*75             # AUCArt,  Mouse / Human, NOAEL 75 ug/kg/d
AUCR75          = ((Mouse_out$AUCR_1.1)/(Human_out$AUCR_1.1))*75                   # AUCR, Mouse / Human, NOAEL 75 ug/kg/d
DailyAMetLiverVol         = ((Mouse_out$DailyAMetLiverVol_1.1)/(Human_out$DailyAMetLiverVol_1.1))*75  # DailyAMetLiverVol,  Mouse / Human, NOAEL 75 ug/kg/d
DailyAMetTotalBW3475 = ((Mouse_out$DailyAMetTotalBW34_1.1)/(Human_out$DailyAMetTotalBW34_1.1))*75 # DailyAMetTotalBW34,  Mouse / Human, NOAEL 75 ug/kg/d

# Summary output 
df_Mouse <-data.frame(Mouse_HEC_AUCArt=c(AUCArt75), Mouse_HEC_AUCR= c(AUCR75), Mouse_HEC_DailyAMetLiverVol =c(DailyAMetLiverVol), Mouse_HEC_DailyAMetTotalBW34 = c(DailyAMetTotalBW3475)) 
sub_df_Mouse <-df_Mouse[,c('Mouse_HEC_AUCArt','Mouse_HEC_AUCR','Mouse_HEC_DailyAMetLiverVol','Mouse_HEC_DailyAMetTotalBW34')] 
Mouse_HEC <- apply(sub_df_Mouse, 2, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))) %>% as.data.frame
write.csv(Mouse_HEC, "outputs/Mouse_HEC.csv", row.names=TRUE)

# tidy
# rm(list = ls())


print(Mouse_AUCArt_HEC   <- quantile(AUCArt75, probs = c(0.025, 0.5, 0.975)))
print(Mouse_AUCR_HEC      <- quantile(AUCR75 , probs = c(0.025, 0.5, 0.975)))

print(Mouse_DailyAMetLiverVol_HEC     <- quantile(DailyAMetLiverVol, probs = c(0.025, 0.5, 0.975)))
print(Mouse_DailyAMetTotalBW34_HEC   <- quantile(DailyAMetTotalBW3475 , probs = c(0.025, 0.5, 0.975)))






