library(foreign)
library(dplyr)
# Data Preparation

# Load Data
setwd("/Users/qg251/Google Drive/Quant3_TA/TA/anes_timeseries_2012_dta")
nes2012 <- read.dta("anes_timeseries_2012_Stata12.dta")
nes2012 <- data.frame(nes2012, stringsAsFactors = )

# Keep data with RD vote
nes2012_RD <- filter(nes2012, postvote_presvtwho == '2. [preload: rep_pcname'|postvote_presvtwho == '1. [preload: dem_pcname]')
# 0 dem, 1 rep
nes2012_RD$vote <- as.numeric(nes2012_RD$postvote_presvtwho) - 6


for (i in 1:2250){
  nes2012_RD[,i][nes2012_RD[,i] == ""] <- NA
}

nes2012_RD$vote
write.dta(nes2012_RD, "nes2012_RDvote.dta")

