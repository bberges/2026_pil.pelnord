## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

# setwd('C:/git/2026_pil.pelnord')

# .libPaths( c( .libPaths()[2], "C:/git/2025_her.27.3a47d_forecast/renv/library/windows/R-4.4/x86_64-w64-mingw32") )  

rm(list=ls())

library(icesTAF)

library(FLCore)
library(stockassessment)
library(FLSAM)
library(FLBRP)
# library(FLSRTMB)
# library(FLRef)

# taf.session()

library(methods)


mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'pil_pelnord2026_sf'

source('utilities_model_config.R')


pil.sams <- new('FLSAMs')

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_pil.pelnord.all.Rdata'))

pg <- range(pil.all)[2]

pil.all@catch.n[,ac(2020)] <- NA
pil.all@landings.n[,ac(2020)] <- NA

# stk <- pil.all
# tun <- pil.all.tun
# pil.all.tun <- tun
pil.all.tun$MA_idx <- pil.all.tun$MA_idx[-c(1,2),]

#pil.all@catch[] <- computeCatch(pil.all)

pil.all.tun <- pil.all.tun[c("MA_idx","nansen","rec_russ",'caton')]#,'caton',"rec_russ",,"MA_ssb","nansen",'caton',"rec_russ"

pil.ctrl <- config_sf_pil.all_cecaf2025(pil.all,
                                           pil.all.tun,
                                           pg)

pil.all@catch.n[ac(0),ac(c(1990:1992,2010))] <- NA
pil.all@catch.n[ac(5),ac(2013)] <- NA
pil.all@catch.n[ac(2),ac(2007:2008)] <- NA

# temp[['caton']] <- pil.sam
# temp[['base']] <- pil.sam
# 
# plot(temp)

### ============================================================================
### Run single fleet model then derive ref points
### ============================================================================

pil.sams[['base']]     <- FLSAM(pil.all,
                     pil.all.tun,
                     pil.ctrl)# return.fit=T



### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_pil.pelnord.all.Rdata'))

pg <- range(pil.all)[2]

pil.all@catch.n[,ac(2020)] <- NA
pil.all@landings.n[,ac(2020)] <- NA

# stk <- pil.all
# tun <- pil.all.tun
# pil.all.tun <- tun
# pil.all.tun$MA_idx <- pil.all.tun$MA_idx[-c(1,2),]

#pil.all@catch[] <- computeCatch(pil.all)

pil.all.tun <- pil.all.tun[c("MA_idx","nansen")]#,'caton',"rec_russ",,"MA_ssb","nansen",'caton',"rec_russ"

pil.ctrl <- config_sf_pil.all_cecaf2025_v4(pil.all,
                                        pil.all.tun,
                                        pg)

pil.all@catch.n[ac(0),ac(c(1990:1992,2010))] <- NA
pil.all@catch.n[ac(5),ac(2013)] <- NA
pil.all@catch.n[ac(2),ac(2007:2008)] <- NA

### ============================================================================
### Run single fleet model then derive ref points
### ============================================================================

pil.sams[['surv_alt']]     <- FLSAM(pil.all,
                                pil.all.tun,
                                pil.ctrl)# return.fit=T


