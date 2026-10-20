## Extract results of interest, write TAF output tables

## Before: results_sf.RData (model)
## After:  fatage.csv, natage.csv, results_sf.RData, summary.csv (output)

rm(list=ls())

library(icesTAF)
library(tidyverse)

library(FLCore)
library(stockassessment)
library(FLSAM)

runName <- 'pil_pelnord2025_sf'

mkdir("output")

data.save       <- file.path(".", "data")
model.save      <- file.path(".",'model')
output.save     <- file.path(".", "output")

terminal.year <- 2024

#cp("model/results_sf.RData", "output")
#cp("model/results_mf.RData", "output")

load(file.path(model.save,paste0(runName,'.RData')))
load(file.path(model.save,paste0(runName,'_retro.RData')))

load(file.path(model.save,paste0('pil_pelnord2025_AB_sf','.RData')))
load(file.path(model.save,paste0('pil_pelnord2025_C_sf','.RData')))


#------------------------------------------
### get survey from age 0
pil.age.all <- read.csv(file.path(data.save,'MA_idx.csv'))

ageRange <- ac(0:5)
yearRange <- an(unique(pil.age.all$year))
FLR_temp <- FLQuant(  array(0,dim = c(length(unique(pil.age.all$age)),
                                      length(min(yearRange):max(yearRange)),
                                      1,1,2,1)),
                      dimnames=list(age=ac(unique(pil.age.all$age)),
                                    year=ac(min(yearRange):max(yearRange)),
                                    unit=c("1e6"),
                                    season='all',
                                    area=c('AB','C')))

FLR_temp[,ac(unique(pil.age.all$year)),,,'AB']  <- 
  as.matrix(subset(pil.age.all,region == 'AB') %>% 
              select(-c(region)) %>%
              pivot_wider(names_from = year, values_from = idx) %>%
              select(-c(age)))

FLR_temp[,ac(unique(pil.age.all$year)),,,'C']   <- 
  as.matrix(subset(pil.age.all,region == 'C') %>% 
              select(-c(region)) %>%
              pivot_wider(names_from = year, values_from = idx) %>%
              select(-c(age)))

MA_idx <- FLIndex(name = 'MA_idx',
                  desc = 'Marocan survey index',
                  index=setPlusGroup(FLR_temp,5,na.rm=T)) # set plus group to 4 to avoid NA in the C area
range(MA_idx)[c('startf', 'endf')] <- c(0.85,0.9) # Nov-Dec
#MA_idx@index@.Data[which(is.na(MA_idx@index))] <- -1
MA_idx <- MA_idx[ageRange,] # start from age 1 as age 0 is flimsy
MA_idx@index[,ac(2020)] <- NA

# all regions combined
FLR_temp <- FLQuant(  array(0,dim = c(length(dimnames(MA_idx)$age),
                                      length(dimnames(MA_idx)$year),
                                      1,1,1,1)),
                      dimnames=list(age=dimnames(MA_idx)$age,
                                    year=dimnames(MA_idx)$year,
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))
MA_idx.all <- FLIndex(name = 'MA_idx',
                      desc = 'Marocan survey index',
                      index=FLR_temp)

MA_idx.all@index[] <- areaSums(MA_idx@index[,,,,c('AB','C')])
range(MA_idx.all)[c('startf', 'endf')] <- c(0.5,0.6)
#MA_idx.all@index@.Data[which(is.na(MA_idx.all@index))] <- -1
MA_idx.all <- MA_idx.all[ageRange,] # start from age 1 as age 0 is flimsy
MA_idx.all@index[,ac(2020)] <- NA

### ============================================================================
### mohn rho
### ============================================================================

df.mr <- data.frame(matrix(ncol = 3, nrow = 0))

df.mr <- rbind(df.mr,cbind(seq((terminal.year-5),terminal.year),
                           t(t(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="ssb")$rho)),
                           rep('ssb',10,1)))
df.mr <- rbind(df.mr,cbind(seq((terminal.year-5),terminal.year),
                           t(t(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="fbar")$rho)),rep('fbar',5,1)))
df.mr <- rbind(df.mr,cbind(seq((terminal.year-5),terminal.year),
                           t(t(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="rec")$rho)),rep('rec',5,1)))

df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="ssb")$rho),'ssb'))
df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="fbar")$rho),'fbar'))
df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(pil.retro,span=5,ref.year=terminal.year,type="rec")$rho),'rec'))

colnames(df.mr) <- c('year','mohn_rho','var')

df.mr <- df.mr %>% pivot_wider(names_from = var, values_from = mohn_rho)

write.taf(df.mr, "output/mohn_rho_sf.csv")

### ============================================================================
### caton
### ============================================================================

df.caton <- as.data.frame(pil.all@catch) %>% select(-c(area,age,unit,season,iter))
df.caton$region <- 'all'

df.temp <- as.data.frame(pil.AB@catch) %>% select(-c(area,age,unit,season,iter))
df.temp$region <- 'AB'

df.caton <- rbind(df.caton,df.temp)

df.temp <- as.data.frame(pil.C@catch) %>% select(-c(area,age,unit,season,iter))
df.temp$region <- 'C'

df.caton <- rbind(df.caton,df.temp)
df.caton <- df.caton %>% pivot_wider(names_from = region,values_from = data)

write.taf(df.caton, "output/caton.csv")

### ============================================================================
### input data
### ============================================================================

my.df <- data.frame(flr2taf(yearMeans(pil.all@mat))) %>% rename(age=Year)
write.taf(my.df, "output/maturity.csv")
my.df <- data.frame(flr2taf(yearMeans(pil.all@m))) %>% rename(age=Year)
write.taf(my.df, "output/natmort.csv")
write.taf(flr2taf(yearMeans(pil.all@harvest.spwn)), "output/propf.csv")
write.taf(flr2taf(yearMeans(pil.all@m.spwn)), "output/propm.csv")

west.df <- rbind(data.frame(plus(flr2taf(pil.all@catch.wt)),region='all',check.names = F),
                  data.frame(plus(flr2taf(pil.AB@catch.wt)),region='AB',check.names = F),
                  data.frame(plus(flr2taf(pil.C@catch.wt)),region='C',check.names = F))

write.taf(west.df, "output/wcatch.csv")
write.taf(west.df, "output/wstock.csv")

# canum
canum.df <- rbind(data.frame(plus(flr2taf(pil.all@catch.n)),region='all',check.names = F),
                  data.frame(plus(flr2taf(pil.AB@catch.n)),region='AB',check.names = F),
                  data.frame(plus(flr2taf(pil.C@catch.n)),region='C',check.names = F))
write.taf(canum.df, "output/catage.csv")

# caton
caton.df <- rbind(data.frame(flr2taf(pil.all@catch),region='all',check.names = F),
                  data.frame(flr2taf(pil.AB@catch),region='AB',check.names = F),
                  data.frame(flr2taf(pil.C@catch),region='C',check.names = F)) %>%
  pivot_wider(names_from = region,values_from = Value)
write.taf(caton.df, "output/caton.csv")

# survey indices
idx.df <- rbind(data.frame(plus(flr2taf(areaSums(MA_idx@index[,,,,c('AB','C')]))),region='all',check.names = F),
                 data.frame(plus(flr2taf(MA_idx@index[,,,,c('AB')])),region='AB',check.names = F),
                 data.frame(plus(flr2taf(MA_idx@index[,,,,c('C')])),region='C',check.names = F))
write.taf(idx.df, "output/survey_MA_idx.csv")

# nansen SSB index
idx.df <- data.frame(flr2taf(pil.all.tun$nansen@index))
write.taf(idx.df, "output/nansen_ssb.csv")

# MA SSB index
idx.df <- rbind(data.frame(flr2taf(pil.all.tun$MA_ssb@index),region='all',check.names = F),
                 data.frame(flr2taf(pil.AB.tun$MA_ssb@index),region='AB',check.names = F),
                 data.frame(flr2taf(pil.C.tun$MA_ssb@index),region='C',check.names = F)) %>%
  pivot_wider(names_from = region,values_from = Value)
write.taf(idx.df, "output/MA_ssb.csv")

# rec index
idx.rec <- rbind(data.frame(flr2taf(pil.all.tun$rec_russ@index),region='all',check.names = F),
                 data.frame(flr2taf(pil.AB.tun$rec_russ@index),region='AB',check.names = F),
                 data.frame(flr2taf(pil.C.tun$rec_russ@index),region='C',check.names = F)) %>%
            pivot_wider(names_from = region,values_from = Value)
write.taf(idx.rec, "output/survey_rec.csv")

### ============================================================================
### assessment outputs
### ============================================================================
# stock numbers
stock.n.df <- rbind(data.frame(plus(flr2taf(pil.all@stock.n)),region='all'),
                    data.frame(plus(flr2taf(pil.AB@stock.n)),region='AB'),
                    data.frame(plus(flr2taf(pil.C@stock.n)),region='C'))
write.taf(stock.n.df, "output/stock.n.csv")

# ssb numbers
stock.df <- rbind(data.frame(flr2taf(computeStock(pil.all)),region='all',check.names = F),
                    data.frame(flr2taf(computeStock(pil.AB)),region='AB',check.names = F),
                    data.frame(flr2taf(computeStock(pil.C)),region='C',check.names = F)) %>% 
            pivot_wider(names_from = region,values_from = Value)
write.taf(stock.df, "output/stock.csv")

# F
harvest.df <- rbind(data.frame(plus(flr2taf(pil.all@harvest)),region='all'),
                    data.frame(plus(flr2taf(pil.AB@harvest)),region='AB'),
                    data.frame(plus(flr2taf(pil.C@harvest)),region='C'))
write.taf(harvest.df, "output/harvest.csv")

# rec numbers
rec.df <- rbind(data.frame(flr2taf(rec(pil.all)),region='all',check.names = F),
                  data.frame(flr2taf(rec(pil.AB)),region='AB',check.names = F),
                  data.frame(flr2taf(rec(pil.C)),region='C',check.names = F)) %>% 
  pivot_wider(names_from = region,values_from = Value)
write.taf(rec.df, "output/rec.csv")

### ============================================================================
### Reference points
### ============================================================================
as.data.frame(refpts)

df.refpts <- rbind(data.frame(flr2taf(refpts),region='all',check.names = F),
              data.frame(flr2taf(refpts.AB),region='AB',check.names = F),
              data.frame(flr2taf(refpts.C),region='C',check.names = F)) %>% 
              pivot_wider(names_from = region,values_from = Value)

write.taf(df.refpts, "output/refpts.csv")

### ============================================================================
### stock status
### ============================================================================

df.ssb <- as.data.frame(ssb(pil.C)/refpts.C['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'C'
df.fbar <- as.data.frame(fbar(pil.C)/refpts.C['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'C'

df.all <- rbind(df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.AB)/refpts.AB['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'AB'
df.fbar <- as.data.frame(fbar(pil.AB)/refpts.AB['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'AB'

df.all <- rbind(df.all,df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.all)/refpts['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'all'
df.fbar <- as.data.frame(fbar(pil.all)/refpts.C['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'all'

df.all <- rbind(df.all,df.ssb,df.fbar)
df.all <- df.all %>% select(-c(age,unit,season,area,iter)) %>% pivot_wider(names_from = quant,values_from = data)

write.taf(df.all, "output/stock_status.csv")