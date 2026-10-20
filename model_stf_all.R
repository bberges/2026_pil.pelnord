## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

rm(list=ls())

library(icesTAF)

library(FLCore)
library(stockassessment)
library(FLSAM)
library(FLBRP)
library(FLSRTMB)
library(FLRef)

# taf.session()

library(methods)


mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'pil_pelnord2026_sf'

load(file.path(modelPath,paste0(runName,'.RData')))

# forecast parameters
niters <- 200
nFutureyrs          <- 5
histMinYr           <- dims(pil.all)$minyear
histMaxYr           <- dims(pil.all)$maxyear
yearCurrent         <- histMinYr:histMaxYr # vector the years
futureMaxYr         <- histMaxYr + nFutureyrs
histPeriod          <- ac(histMinYr:histMaxYr)
projPeriod          <- ac((histMaxYr+1):futureMaxYr)

out.SRR <- attr(pil.sr.rick,'SV')

rec_residuals <- FLQuant(NA, dimnames = list(year=projPeriod, iter=1:niters))
rec_residuals <- rec_residuals %=% out.SRR$sigmaR
rec_residuals <-   rlnorm(niters,
                          meanlog=0,
                          sdlog=rec_residuals[,projPeriod])

pil.stf <- stf(pil.all, nyears = nFutureyrs)
pil.stf <- propagate(pil.stf, niters)

run.fwd <- new('FLStocks')

#-------------------------------------------------------------------------------
# stf no fish
#-------------------------------------------------------------------------------
df.ctrl <- data.frame(year = an(projPeriod),
                      quant = "catch",
                      value = 0)
fwd.ctrl <- fwdControl(df.ctrl)


run.fwd[['nf']] <- fwd(pil.stf,
                        control = fwd.ctrl,
                        sr = pil.sr.rick,
                        residuals = rec_residuals)

#-------------------------------------------------------------------------------
# stf with F status quo
#-------------------------------------------------------------------------------
f_status_quo <- mean(fbar(pil.all)[,as.character(2022:2024)])

df.ctrl <- data.frame(year = an(projPeriod),
                          quant = "f",
                          value = f_status_quo)
fwd.ctrl <- fwdControl(df.ctrl)
fwd.ctrl <- propagate(fwd.ctrl,niters)
fwd.ctrl@iters[,"value",] <- fwd.ctrl@iters[,"value",]*rlnorm(nFutureyrs*niters, meanlog = 0, sdlog=0.15)

run.fwd[['fsq']] <- fwd(pil.stf,
                        control = fwd.ctrl,
                        sr = pil.sr.rick,
                        residuals = rec_residuals)

#-------------------------------------------------------------------------------
# stf with catch status quo
#-------------------------------------------------------------------------------
C_status_quo <- mean(catch(pil.all)[,as.character(2024)])

df.ctrl <- data.frame(year = an(projPeriod),
                      quant = "catch",
                      value = C_status_quo)
fwd.ctrl <- fwdControl(df.ctrl)
fwd.ctrl <- propagate(fwd.ctrl,niters)
fwd.ctrl@iters[,"value",] <- fwd.ctrl@iters[,"value",]*rlnorm(nFutureyrs*niters, meanlog = 0, sdlog=0.15)

run.fwd[['Csq']] <- fwd(pil.stf,
                        control = fwd.ctrl,
                        sr = pil.sr.rick,
                        residuals = rec_residuals)

#-------------------------------------------------------------------------------
# stf with recovery to SBMSY
#-------------------------------------------------------------------------------

fwd.ctrl  <- fwdControl(data.frame(year = an(projPeriod),
                                   quant = "ssb_spawn",
                                   value = c(refpts['SBMSY'])))

run.fwd[['SBMSY']] <- fwd(pil.stf,
                        control = fwd.ctrl,
                        sr = pil.sr.rick,
                        residuals = rec_residuals)

#plot(window(run.fwd[['SBMSY']], start=2010))

#-------------------------------------------------------------------------------
# stf with recovery to SBMSY with 10% down cap  on catches
#-------------------------------------------------------------------------------

rel_catch_bound <- 0.10

df.ctrl <- data.frame(  year = rep(an(projPeriod), 2),
                        relYear =c(rep(NA,nFutureyrs), an(projPeriod)-1),
                        quant = c(rep("ssb_spawn",nFutureyrs), rep("catch",nFutureyrs)),
                        value = c(rep(refpts['SBMSY'],nFutureyrs),rep(NA,nFutureyrs)),
                        max = c(rep(NA,nFutureyrs), rep(NA, nFutureyrs)),
                        min = c(rep(NA,nFutureyrs), rep(1-rel_catch_bound, nFutureyrs)))
fwd.ctrl <- fwdControl(df.ctrl)

run.fwd[['SBMSY_C.cap10']] <- fwd(pil.stf,
                                control = fwd.ctrl,
                                sr = pil.sr.rick,
                                residuals = rec_residuals)


#-------------------------------------------------------------------------------
# stf with recovery to SBMSY with 10% down cap  on catches
#-------------------------------------------------------------------------------

rel_catch_bound <- 0.20

df.ctrl <- data.frame(  year = rep(an(projPeriod), 2),
                        relYear =c(rep(NA,nFutureyrs), an(projPeriod)-1),
                        quant = c(rep("ssb_spawn",nFutureyrs), rep("catch",nFutureyrs)),
                        value = c(rep(refpts['SBMSY'],nFutureyrs),rep(NA,nFutureyrs)),
                        max = c(rep(NA,nFutureyrs), rep(NA, nFutureyrs)),
                        min = c(rep(NA,nFutureyrs), rep(1-rel_catch_bound, nFutureyrs)))
fwd.ctrl <- fwdControl(df.ctrl)

run.fwd[['SBMSY_C.cap20']] <- fwd(pil.stf,
                                  control = fwd.ctrl,
                                  sr = pil.sr.rick,
                                  residuals = rec_residuals)

#-------------------------------------------------------------------------------
# stf with recovery to SBMSY with 10% down cap  on catches
#-------------------------------------------------------------------------------

rel_catch_bound <- 0.30

df.ctrl <- data.frame(  year = rep(an(projPeriod), 2),
                        relYear =c(rep(NA,nFutureyrs), an(projPeriod)-1),
                        quant = c(rep("ssb_spawn",nFutureyrs), rep("catch",nFutureyrs)),
                        value = c(rep(refpts['SBMSY'],nFutureyrs),rep(NA,nFutureyrs)),
                        max = c(rep(NA,nFutureyrs), rep(NA, nFutureyrs)),
                        min = c(rep(NA,nFutureyrs), rep(1-rel_catch_bound, nFutureyrs)))
fwd.ctrl <- fwdControl(df.ctrl)

run.fwd[['SBMSY_C.cap30']] <- fwd(pil.stf,
                                  control = fwd.ctrl,
                                  sr = pil.sr.rick,
                                  residuals = rec_residuals)


#-------------------------------------------------------------------------------
# saving results
#-------------------------------------------------------------------------------

save(run.fwd,
     file=file.path(modelPath,paste0(runName,'_stf.RData')))

#-------------------------------------------------------------------------------
# stf with recovery to SBMSY with 10% down cap  on catches with stochasticity
#-------------------------------------------------------------------------------

# rel_catch_bound <- 0.10
# 
# df.ctrl <- data.frame(
#   year = rep(an(projPeriod), 2),
#   relYear =c(rep(NA,nFutureyrs), an(projPeriod)-1),
#   quant = c(rep("ssb_spawn",nFutureyrs), rep("catch",nFutureyrs)),
#   value = c(rep(refpts['SBMSY'],nFutureyrs),rep(NA,nFutureyrs)),
#   max = c(rep(NA,nFutureyrs), rep(NA, nFutureyrs)),
#   min = c(rep(NA,nFutureyrs), rep(1-rel_catch_bound, nFutureyrs)))
# fwd.ctrl <- fwdControl(df.ctrl)
# fwd.ctrl <- propagate(fwd.ctrl,niters)
# 
# idxCol <- iter(fwd.ctrl,1)@target$quant == 'catch'
# fwd.ctrl@iters[idxCol,"min",] <- fwd.ctrl@iters[idxCol,"min",]*rlnorm(nFutureyrs*niters, meanlog = 0, sdlog=0.15)
# 
# run.fwd[['SBMSY_C.cap.rand']] <- fwd(pil.stf,
#                                 control = fwd.ctrl,
#                                 sr = pil.sr.rick,
#                                 residuals = rec_residuals)

# plot(window(run.fwd[['SBMSY_C.cap.rand']], start=2010))
# 
# plot(window(run.fwd, start=2010))