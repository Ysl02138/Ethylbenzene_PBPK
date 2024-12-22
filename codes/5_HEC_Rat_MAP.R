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
# Rat -------------------------------------------------------------------

Ratsim1.1 <- fread("outputs/EBRatMCMC_3365.out") 
Ratsim2.1 <- fread("outputs/EBRatMCMC_6734.out") 
Ratsim3.1 <- fread("outputs/EBRatMCMC_4880.out") 
Ratsim4.1 <- fread("outputs/EBRatMCMC_5916.out") 

# Ratsim1.1 <- fread("outputs/EBRatMCMC_2896.out")
# Ratsim2.1 <- fread("outputs/EBRatMCMC_5656.out")
# Ratsim3.1 <- fread("outputs/EBRatMCMC_7099.out")
# Ratsim4.1 <- fread("outputs/EBRatMCMC_1375.out")

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
tmp.Rat_x %>% write.table(file="Rat.HEC.dat", row.names=T, sep="\t")

vld <- "./mcsim.EBPBPK.model.exe MCSim/Rat_HEC.in"
system(vld)

Rat_out <- read.delim("Rat.HEC.out")

# Tidy data with median and 95% confidence interval
vars <- names(Rat_out)
index <- which(vars == "AUCArt_1.1" | vars == "DailyAMetTotalBW34_1.1")
Rat_Summary <- apply(Rat_out[index[1]:index[2]], 2, quantile,  c(0.025, 0.5, 0.975)) %>% t()
colnames(Rat_Summary) <- c("median", "LCL", "UCL")
df <- as.data.frame(Rat_Summary)


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
Human_out = as.data.frame(Human_out)
# Tidy data with median and 95% confidence interval
vars <- names(Human_out)
index <- which(vars == "AUCArt_1.1" | vars == "DailyAMetTotalBW34_1.1")
Human_Summary <- apply(Human_out[index[1]:index[2]], 2, quantile,  c(0.025, 0.5, 0.975)) %>% t()
colnames(Human_Summary) <- c("median", "LCL", "UCL")
df <- as.data.frame(Human_Summary)


#####################################################################
## Calculated the HEC 

## Rat 
AUCArt75       = ((Rat_out$AUCArt_1.1)/(Human_out$AUCArt_1.1))*75             # AUCArt,  Rat / Human, NOAEL 75 ug/kg/d
AUCR75          = ((Rat_out$AUCR_1.1)/(Human_out$AUCR_1.1))*75                   # AUCR, Rat / Human, NOAEL 75 ug/kg/d
DailyAMetLiverVol75   = ((Rat_out$DailyAMetLiverVol_1.1)/(Human_out$DailyAMetLiverVol_1.1))*75     # DailyAMetLiverVol,  Rat / Human, NOAEL 75 ug/kg/d
DailyAMetTotalBW3475 = ((Rat_out$DailyAMetTotalBW34_1.1)/(Human_out$DailyAMetTotalBW34_1.1))*75 # DailyAMetTotalBW34,  Rat / Human, NOAEL 75 ug/kg/d

# Summary output 
df_Rat <-data.frame(Rat_HEC_AUCArt=c(AUCArt75), Rat_HEC_AUCR= c(AUCR75), Rat_HEC_DailyAMetLiverVol =c(DailyAMetLiverVol75), Rat_HEC_DailyAMetTotalBW34 = c(DailyAMetTotalBW3475)) 
sub_df_Rat <-df_Rat[,c('Rat_HEC_AUCArt','Rat_HEC_AUCR','Rat_HEC_DailyAMetLiverVol','Rat_HEC_DailyAMetTotalBW34')] 
Rat_HEC <- apply(sub_df_Rat, 2, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))) %>% as.data.frame
write.csv(Rat_HEC, "outputs/Rat_HEC.csv", row.names=TRUE)

# tidy
# rm(list = ls())


print(Rat_AUCArt_HEC   <- quantile(AUCArt75, probs = c(0.025, 0.5, 0.975)))
print(Rat_AUCR_HEC  <- quantile(AUCR75 , probs = c(0.025, 0.5, 0.975)))

print(Rat_DailyAMetLiverVol_HEC   <- quantile(DailyAMetLiverVol75, probs = c(0.025, 0.5, 0.975)))
print(Rat_DailyAMetTotalBW34_HEC   <- quantile(DailyAMetTotalBW3475 , probs = c(0.025, 0.5, 0.975)))






