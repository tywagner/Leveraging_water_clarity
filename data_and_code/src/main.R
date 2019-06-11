
### main.R
###
### Meridith Bartley, Ephraim Hanks, Erin Schliep, Tyler Wagner, Nathan Wikle
###
### The main script to recreate the analysis and results from the manuscript,
###   "Increasing accuracy of lake nutrient predictions in thousands of 
###     lakes by leveraging water clarity data," 2019.


### 1) Source in necessary functions:
source("./src/functions.R")

### 2) Perform the analysis on the seven data sets. MCMC samples are saved
###     to a run#.rds file in ./data, while full response post. pred. samples
###     are in ./data/pred#_180627 and conditional preds are in 
###     ./data/pred#_cond.rds.

###   Note: This takes a long time!

# holdout set 1: random25 
source("./src/run1.R")
# holdout set 2: random75
source("./src/run2.R")

rm(list = ls())

### 3) Compute root mean squared error (RMSE) for TN, TP, and CHL.

# load full prediction results
pred1 <- readRDS("pred1_180627.Rdata")
pred2 <- readRDS("pred2_180627.Rdata")

# load conditional prediction results
pred1.cond <- readRDS("pred1_cond_180627.Rdata")
pred2.cond <- readRDS("pred2_cond_180627.Rdata")

rmse.values.full <- rbind(rmse(pred1),
                          rmse(pred2))

rmse.values.cond <- rbind(rmse(pred1.cond),
                          rmse(pred2.cond))

mpe.values.full <- rbind(percent.error(pred1),
                         percent.error(pred2))

mpe.values.cond <- rbind(percent.error(pred1.cond),
                         percent.error(pred2.cond))

saveRDS(rmse.values.full, "./output/rmse_full.rds")
saveRDS(rmse.values.cond, "./output/rmse_cond.rds")
saveRDS(mpe.values.full, "./output/mpe_full.rds")
saveRDS(mpe.values.cond, "./output/mpe_cond.rds")
