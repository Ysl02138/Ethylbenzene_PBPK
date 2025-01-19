# Tidy
rm(list=ls())

# Load packages ----------------------------------------------------------------
pacman::p_load(pksensi, ggplot2, ggpubr,cowplot, reshape2, dplyr, viridis)

# Create mod.exe ---------------------------------------------------------------
source("MCSim/function.R")
set_PATH()
# makemod()

# Compile model code -----------------------------------------------------------
mName <- "EBPBPK.model" 
# makemcsim(mName)

# Parameter range
load("outputs/EBMouse_mcmc.RData")
out_vars <- dimnames(Mouse_mcmc)[[3]] 
nd1 <- dim(Mouse_mcmc)[1] * 4
nd2 <- dim(Mouse_mcmc)[3]
dim(Mouse_mcmc) <- c(nd1,nd2)
colnames(Mouse_mcmc) <- out_vars
str <- which(colnames(Mouse_mcmc) == "M_lnPBC(1)")
end <- which(colnames(Mouse_mcmc) == "M_lnKMRptC(1)")
pop_min <- Mouse_mcmc[,str:end] |> apply(2, min)
pop_max <- Mouse_mcmc[,str:end] |> apply(2, max)
pop_max - pop_min

# Parameter space generating
params <- c("lnPBC", "lnPFatC", "lnPLivC", "lnPRptC", 
  "lnPSptC", "lnPLungC", "lnVmaxC", "lnKMC", "lnVmax2C", 
  "lnKM2C", "lnVmaxLungC", "lnKMLungC", "lnVmaxRptC", "lnKMRptC")
q <- rep("qunif", 14)
q.arg <- list(list(pop_min[1], pop_max[1]), 
  list(pop_min[2], pop_max[2]), 
  list(pop_min[3], pop_max[3]), 
  list(pop_min[4], pop_max[4]), 
  list(pop_min[5], pop_max[5]), 
  list(pop_min[6], pop_max[6]), 
  list(pop_min[7], pop_max[7]), 
  list(pop_min[8], pop_max[8]), 
  list(pop_min[9], pop_max[9]), 
  list(pop_min[10], pop_max[10]), 
  list(pop_min[11], pop_max[11]), 
  list(pop_min[12], pop_max[12]), 
  list(pop_min[13], pop_max[13]), 
  list(pop_min[14], pop_max[14])
)
set.seed(15)
x <- rfast99(params = params, n = 10000, q = q, q.arg = q.arg, replicate = 1)

# Single dose simulation
conditions <- c("BWmeas = 0.03;",
                "Species = 3;", 
                "Sex = 2;", 
                "ExpoInduc = 1;",
                "expoday = PerDose (75, 24, 0, 24);",
                "expowk  = PerDose (1.0, 168,  0, 168);",
                "expodur = PerDose (1.0, 336, 0, 336);"
                )

vars <- c("AUCArt", "AUCR",  "DailyAMetLiverVol","DailyAMetTotalBW34")
times <- c(336)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
mSI <- out$mSI
iSI <- out$iSI
dim(mSI) <- dim(iSI) <- c(14,4)
colnames(mSI) <- colnames(iSI) <- vars
rownames(mSI) <- rownames(iSI) <- params
x1 <- melt(mSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "main effect", specice = "Mouse", exposure = "75 ppm")
x2 <- melt(iSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "interaction", specice = "Mouse", exposure = "75 ppm")

conditions[5] <- "expoday = PerDose (750, 24, 0, 24);"
conditions[4] <- "ExpoInduc = 3;"
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
mSI <- out$mSI
iSI <- out$iSI
dim(mSI) <- dim(iSI) <- c(14,4)
colnames(mSI) <- colnames(iSI) <- vars
rownames(mSI) <- rownames(iSI) <- params
x3 <- melt(mSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "main effect", specice = "Mouse", exposure = "750 ppm")
x4 <- melt(iSI) |> `colnames<-`(c("parameter", "variable", "index")) |>
  mutate(order = "interaction", specice = "Mouse", exposure = "750 ppm")

x_total <- do.call(rbind, list(x1, x2, x3 ,x4))
x_total$order <- factor(x_total$order, level=c("interaction", "main effect"))

set_theme <- theme(
  legend.position  = "top",
  legend.title     = element_blank(),
  axis.title       = element_blank(),
)


ggplot(x_total, aes(x = parameter, y = index, fill=order)) +
  geom_bar(position = "stack", stat = "identity") + 
  coord_flip() +
  facet_grid(exposure ~ variable) + 
  theme_bw() +
  ggtitle("Mouse") +
  scale_fill_viridis(discrete = TRUE, direction = -1, end = 0.8) + 
  set_theme
ggsave(file = "plots/suppl/Supplemental_Figure_S10_Mouse EB GSA.jpg", height = 12, width = 20, dpi = 600)

dev.off()  
