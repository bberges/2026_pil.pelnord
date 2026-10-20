### ============================================================================
### library loading and setting up paths
### ============================================================================

rm(list=ls())

library(icesTAF)
library(FLCore)
library(tidyverse)
library(FSA)

mkdir("data")

dataPath     <- file.path(".", "data")
dataSource   <- file.path(".","bootstrap", "data")

### ============================================================================
### read data
### ============================================================================

matprop             <- read.csv(file.path(dataSource,'matprop.csv'),sep = ';')
ALK_catch           <- read.csv(file.path(dataSource,'ALK_catch.csv'),sep = ';',check.names = F)
canum_AB            <- read.csv(file.path(dataSource,'canum_AB.csv'),sep = ';',check.names = F)
canum_C             <- read.csv(file.path(dataSource,'canum_C.csv'),sep = ';',check.names = F)
catch_LF_AB         <- read.csv(file.path(dataSource,'catch_LF_AB.csv'),sep = ';',check.names = F)
catch_LF_C          <- read.csv(file.path(dataSource,'catch_LF_C.csv'),sep = ';',check.names = F)
catch_LF_canary     <- read.csv(file.path(dataSource,'catch_LF_canary.csv'),sep = ';',check.names = F)
catch_LF_north      <- read.csv(file.path(dataSource,'catch_LF_north.csv'),sep = ';',check.names = F)
caton_AB            <- read.csv(file.path(dataSource,'caton_AB.csv'),sep = ';',check.names = F)
caton_C             <- read.csv(file.path(dataSource,'caton_C.csv'),sep = ';',check.names = F)
caton_north         <- read.csv(file.path(dataSource,'caton_north.csv'),sep = ';',check.names = F)
survey_nansen_SSB   <- read.csv(file.path(dataSource,'survey_nansen_biomass.csv'),sep = ';',check.names = F)
survey_MA_SSB       <- read.csv(file.path(dataSource,'survey_MA_biomass.csv'),sep = ';',check.names = F)
survey_MA_LF_AB     <- read.csv(file.path(dataSource,'survey_MA_LF_AB.csv'),sep = ';',check.names = F)
survey_MA_LF_C      <- read.csv(file.path(dataSource,'survey_MA_LF_C.csv'),sep = ';',check.names = F)
survey_rec_northCB  <- read.csv(file.path(dataSource,'survey_rec_northCB.csv'),sep = ';',check.names = F)
survey_rec_southCB  <- read.csv(file.path(dataSource,'survey_rec_southCB.csv'),sep = ';',check.names = F)
west_AB             <- read.csv(file.path(dataSource,'west_AB.csv'),sep = ';',check.names = F)
west_C              <- read.csv(file.path(dataSource,'west_C.csv'),sep = ';',check.names = F)
natmort             <- read.csv(file.path(dataSource,'natmort.csv'),sep = ';',check.names = F)

# tidy up recruitment survey
survey_rec_southCB[is.na(survey_rec_southCB)] <- 0
survey_rec_northCB[is.na(survey_rec_northCB)] <- 0
survey_rec <- survey_rec_northCB+survey_rec_southCB
survey_rec[survey_rec == 0] <- NA
survey_rec <- survey_rec[,2:dim(survey_rec)[2]]
survey_rec_northCB <- survey_rec_northCB[,2:dim(survey_rec_northCB)[2]]
survey_rec_northCB[survey_rec_northCB == 0] <- NA
survey_rec_southCB <- survey_rec_southCB[,2:dim(survey_rec_southCB)[2]]
survey_rec_southCB[survey_rec_southCB == 0] <- NA

# combine survey at length for MA survey
survey_MA_LF_AB <- survey_MA_LF_AB %>% 
  pivot_longer(!L,names_to = "year", values_to = "idx") %>% 
  rename(len=L)
survey_MA_LF_AB$region <- 'AB'

survey_MA_LF_C <- survey_MA_LF_C %>% 
  pivot_longer(!L,names_to = "year", values_to = "idx") %>% 
  rename(len=L)
survey_MA_LF_C$region <- 'C'

survey_MA <- rbind(survey_MA_LF_C,survey_MA_LF_AB)

### ============================================================================
### compute indices at age from ALKs
### ============================================================================

flagFirst <- T
for(idxRegion in unique(survey_MA$region)){
  for(idxYear in unique(survey_MA$year)){
    idx.sub <- subset(survey_MA,year == idxYear & region == idxRegion) %>% 
                select(-c(year))
    idx.sub$LCat <- idx.sub$len
    idx.sub$idx[is.na(idx.sub$idx)] <- 0
    
    flagFirst.sub <- T
    for(idxlen in unique(idx.sub$len)){
      idx.current <- which(idx.sub$len == idxlen)
      temp <- idx.sub[rep(which(idx.sub$len == idxlen), idx.sub$idx[idx.current]),]
      
      if(length(temp) != 0){
        if(flagFirst.sub){
          ind.sub <- temp
          flagFirst.sub <- F
        }else{
          ind.sub <- rbind(ind.sub,temp)
        }
      }
    }
    
    if(dim(ind.sub)[1] != 0){
      ind.sub <- ind.sub %>% select(-c(idx,region))
      ind.sub$age <- an(NA)
      ind.sub$ID <- 1:dim(ind.sub)[1]
      
      ALK.sub <- subset(ALK_catch,region == idxRegion)
      
      if(min(ALK.sub$year) > an(idxYear)){
        ALK.sub <- subset(ALK.sub,year == min(ALK.sub$year))
      }else{
        ALK.sub <- subset(ALK.sub,year == idxYear)
      }
      
      # fill lower part
      n.fill <- length(seq(min(idx.sub$len),min(ALK.sub$L)-0.5,by=0.5))
      ALK.fill.down <- ALK.sub[rep(which(ALK.sub$L == min(ALK.sub$L)), n.fill), ]
      ALK.fill.down$L <- seq(min(idx.sub$len),min(ALK.sub$L)-0.5,by=0.5)
      
      # fill higher part
      n.fill <- length(seq(max(ALK.sub$L)+0.5,max(idx.sub$len),by=0.5))
      ALK.fill.up <- ALK.sub[rep(which(ALK.sub$L == max(ALK.sub$L)), n.fill), ]
      ALK.fill.up$L <- seq(max(ALK.sub$L)+0.5,max(idx.sub$len),by=0.5)
      
      ALK.sub <- rbind(ALK.sub,ALK.fill.down,ALK.fill.up) %>% arrange(L)
      
      ALK.sub <- ALK.sub %>%
        select(-c(region, year)) %>%
        pivot_longer(!L,names_to = "age", values_to = "Freq") %>%
        rename(LCat=L)
      ALK.sub$Freq[is.na(ALK.sub$Freq)] <- 0
      ALK.sub$Freq <- ALK.sub$Freq/100
      pil.key <- xtabs(Freq~LCat+age,data=subset(ALK.sub,!is.na(age)))
      
      brks <- unique(idx.sub$len)
      
      pil.len <- alkIndivAge(pil.key,age~len,data=ind.sub,breaks=brks)
      
      pil.age <- pil.len %>% group_by(age) %>% summarise(idx=n())
      pil.age$year    <- idxYear
      pil.age$region  <- idxRegion
      
      if(flagFirst){
        pil.age.all <- pil.age
        flagFirst <- F
      }else{
        pil.age.all <- rbind(pil.age.all,pil.age)
      }
    }
  }
}

write.csv(pil.age.all,file.path(dataPath,'MA_idx.csv'),row.names = F)

### ============================================================================
### prepare the FLStock object
### ============================================================================
minYear <- 1990
maxYear <- 2025
yearVec <- minYear:maxYear
nYears  <- length(yearVec)
minAge  <- 0
maxAge  <- 7
ageVec  <- minAge:maxAge
nAges   <- length(ageVec)
nIter   <- 1
regions <- c('north','AB','C')
nRegions <- length(regions)
pg <- 5
pg_survey <- 5

#### stock object
FLR_temp <- FLQuant( array(0,dim = c(nAges,nYears,1,1,nRegions,nIter)),
                     dimnames=list( age=as.character(ageVec),
                                    year=as.character(yearVec),
                                    unit=c("1e6"),
                                    season='all',
                                    area=regions),
                     units='1e6')

pil.pelnord <- FLStock(name = 'pil', 
                       desc = 'COPACE Pel-Nord sardine',
                       FLR_temp)

# caton
pil.pelnord@catch[,,,,'north']  <- an(caton_north[1,])
pil.pelnord@catch[,,,,'AB']     <- an(caton_AB[1,])
pil.pelnord@catch[,,,,'C']      <- an(caton_C[1,])

pil.pelnord@landings <- pil.pelnord@catch

pil.pelnord@discards[] <- 0

# weca
pil.pelnord@catch.wt[,,,,'north'] <- as.matrix(west_AB[,2:dim(west_AB)[2]])
pil.pelnord@catch.wt[,,,,'AB'] <- as.matrix(west_AB[,2:dim(west_AB)[2]])
pil.pelnord@catch.wt[,,,,'C'] <- as.matrix(west_C[,2:dim(west_C)[2]])

pil.pelnord@landings.wt <- pil.pelnord@catch.wt
pil.pelnord@discards.wt <- pil.pelnord@catch.wt

# west
pil.pelnord@stock.wt[,,,,'north'] <- as.matrix(west_AB[,2:dim(west_AB)[2]])
pil.pelnord@stock.wt[,,,,'AB'] <- as.matrix(west_AB[,2:dim(west_AB)[2]])
pil.pelnord@stock.wt[,,,,'C'] <- as.matrix(west_C[,2:dim(west_C)[2]])

# canum
pil.pelnord@catch.n[,,,,'AB'] <- as.matrix(canum_AB[,2:dim(canum_AB)[2]])
pil.pelnord@catch.n[,,,,'C'] <- as.matrix(canum_C[,2:dim(canum_C)[2]])

pil.pelnord@landings.n <- pil.pelnord@catch.n

pil.pelnord@discards.n[] <- 0

# mat
pil.pelnord@mat[] <- as.matrix(matprop[,2])

# spawning time
pil.pelnord@m.spwn[]        <- 0.5
pil.pelnord@harvest.spwn[]  <- 0.5

pil.pelnord@m[] <- natmort$m

pil.pelnord@stock.wt@units <- 'kg'
pil.pelnord@harvest@units <- 'f'

# FLStocks combining region AB and region C
# stock object
FLR_temp <- FLQuant( array(0,dim = c(nAges,nYears,1,1,1,1)),
                     dimnames=list( age=as.character(ageVec),
                                    year=as.character(yearVec),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'),
                     units='1e6')

pil.all <- FLStock(name = 'pil', 
                       desc = 'COPACE Pel-Nord sardine all regions',
                       FLR_temp)

pil.all@catch <- areaSums(pil.pelnord@catch[,,,,c('AB','C')])
pil.all@landings <- pil.all@catch
pil.all@discards[] <- 0

pil.all@catch.wt <- areaMeans(pil.pelnord@catch.wt[,,,,c('AB','C')])
pil.all@landings.wt <- pil.all@catch.wt
pil.all@discards.wt <- pil.all@catch.wt

pil.all@stock.wt <- areaMeans(pil.pelnord@stock.wt[,,,,c('AB','C')])

pil.all@catch.n <- areaSums(pil.pelnord@catch.n[,,,,c('AB','C')])
pil.all@landings.n <- pil.all@catch.n
pil.all@discards.n[] <- 0

pil.all@mat[] <- as.matrix(matprop[,2])

pil.all@m.spwn[]        <- 0.5
pil.all@harvest.spwn[]  <- 0.5

pil.all@m[] <- natmort$m # taken from Sardina pilchardus in divisions 8.c and 9.a

pil.all@stock.wt@units <- 'kg'
pil.all@harvest@units <- 'f'

# set plus group
pil.pelnord <- setPlusGroup(pil.pelnord,pg,na.rm=T)
pil.all <- setPlusGroup(pil.all,pg,na.rm=T)

range(pil.all)[c('minfbar', 'maxfbar')] <- c(2,4)
range(pil.pelnord)[c('minfbar', 'maxfbar')] <- c(2,4)

# start ages to 1
# pil.pelnord <- pil.pelnord[ac(1:pg),]
# pil.all     <- pil.all[ac(1:pg),]

# create FLStocks per region
pil.AB <- pil.pelnord[,,,,'AB']
pil.C <- pil.pelnord[,,,,'C']

### ============================================================================
### recruitment survey
### ============================================================================

### rec survey
ageRange <- ac(0)
FLR_temp <- FLQuant(  array(0,dim = c(2,length(colnames(survey_rec)),1,1,1,1)),
                              dimnames=list(age=as.character(c(0:1)),
                                            year=colnames(survey_rec),
                                            unit=c("1e6"),
                                            season='all',
                                            area='all'))

# combined rec index
rec_russ <- FLIndex(name = 'russian_rec', 
                    desc = 'Russian recruitment survey S+N',
                    index=FLR_temp)
rec_russ <- rec_russ[ac(0)]
rec_russ@index[] <- apply(survey_rec, c(2), sum,na.rm=T)#as.matrix(survey_rec)
range(rec_russ)[c('startf', 'endf')] <- c(0.8,0.9) # Oct-Nov
rec_russ@index@.Data[which(rec_russ@index == 0)] <- NA

#rec_russ@index@.Data[which(is.na(rec_russ@index))] <- -1
# rec_russ <- rec_russ[ageRange,] # reduce to age 0
# rec_russ@index[ac(0)] <- quantMeans(rec_russ@index)
# rec_russ <- rec_russ[ac(0)]

# df.plot <- as.data.frame(rec_russ@index)
# ggplot(df.plot,aes(x=year,y=data,col=as.factor(age)))+
#   geom_point()

# south rec index
rec_russ.S <- FLIndex(name = 'russian_rec', 
                    desc = 'Russian recruitment survey South',
                    index=FLR_temp)
rec_russ.S <- rec_russ.S[ac(0)]
rec_russ.S@index[] <- apply(survey_rec_southCB, c(2), sum,na.rm=T)#as.matrix(survey_rec)
range(rec_russ.S)[c('startf', 'endf')] <- c(0.8,0.9) # Oct-Nov
#rec_russ@index@.Data[which(is.na(rec_russ@index))] <- -1
rec_russ.S@index@.Data[which(rec_russ.S@index == 0)] <- NA

# north rec index
rec_russ.N <- FLIndex(name = 'russian_rec', 
                      desc = 'Russian recruitment survey North',
                      index=FLR_temp)
rec_russ.N <- rec_russ.N[ac(0)]
rec_russ.N@index[] <- apply(survey_rec_northCB, c(2), sum,na.rm=T)#as.matrix(survey_rec)
range(rec_russ.N)[c('startf', 'endf')] <- c(0.8,0.9) # Oct-Nov
#rec_russ@index@.Data[which(is.na(rec_russ@index))] <- -1
rec_russ.N@index@.Data[which(rec_russ.N@index == 0)] <- NA

### ============================================================================
### nansen survey (SSB)
### ============================================================================
FLR_temp <- FLQuant(  array(0,dim = c(1,dim(survey_nansen_SSB)[1],1,1,1,1)),
                      dimnames=list(age='all',
                                    year=ac(survey_nansen_SSB$year),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))

nansen <- FLIndex(name = 'nansen_ssb',
                  desc = 'nansen SSB index',
                  index=FLR_temp)

nansen@index[] <- survey_nansen_SSB$idx
range(nansen)[c('startf', 'endf')] <- c(0.75,0.9) # Sept-Nov
#nansen@index@.Data[which(is.na(nansen@index))] <- -1

### ============================================================================
### Marocan survey (SSB)
### ============================================================================

FLR_temp <- FLQuant(  array(0,dim = c(1,dim(survey_MA_SSB)[1],1,1,3,1)),
                      dimnames=list(age='all',
                                    year=ac(survey_MA_SSB$year),
                                    unit=c("1e6"),
                                    season='all',
                                    area=c('north','AB','C')))

MA_ssb <- FLIndex(name = 'MA_ssb',
                  desc = 'Marocan SSB index',
                  index=FLR_temp)

MA_ssb@index[,,,,'AB']      <- survey_MA_SSB$AB
MA_ssb@index[,,,,'C']       <- survey_MA_SSB$C
MA_ssb@index[,,,,'north']   <- survey_MA_SSB$north
range(MA_ssb)[c('startf', 'endf')] <- c(0.85,0.9) # Nov-Dec
#MA_ssb@index@.Data[which(is.na(MA_ssb@index))] <- -1

# all regions combined
FLR_temp <- FLQuant(  array(0,dim = c(1,dim(survey_MA_SSB)[1],1,1,1,1)),
                      dimnames=list(age='all',
                                    year=ac(survey_MA_SSB$year),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))

MA_ssb.all <- FLIndex(name = 'MA_ssb',
                      desc = 'Marocan SSB index',
                      index=FLR_temp)

MA_ssb.all@index[] <- areaSums(MA_ssb@index[,,,,c('AB','C')])
range(MA_ssb.all)[c('startf', 'endf')] <- c(0.5,0.6)
#MA_ssb.all@index@.Data[which(is.na(MA_ssb.all@index))] <- -1

### ============================================================================
### Marocan survey index at age
### ============================================================================

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

FLR_temp[,unique(pil.age.all$year),,,'AB']  <- 
                                      as.matrix(subset(pil.age.all,region == 'AB') %>% 
                                      select(-c(region)) %>%
                                      pivot_wider(names_from = year, values_from = idx) %>%
                                      select(-c(age)))

FLR_temp[,unique(pil.age.all$year),,,'C']   <- 
            as.matrix(subset(pil.age.all,region == 'C') %>% 
            select(-c(region)) %>%
            pivot_wider(names_from = year, values_from = idx) %>%
            select(-c(age)))

MA_idx <- FLIndex(name = 'MA_idx',
                  desc = 'Marocan survey index',
                  index=setPlusGroup(FLR_temp,5,na.rm=T)) # set plus group to 4 to avoid NA in the C area
range(MA_idx)[c('startf', 'endf')] <- c(0.85,0.9) # Nov-Dec
#MA_idx@index@.Data[which(is.na(MA_idx@index))] <- -1
# MA_idx <- MA_idx[ageRange,] # start from age 1 as age 0 is flimsy
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
                      index=setPlusGroup(FLR_temp,5,na.rm=T))

MA_idx.all@index[] <- areaSums(MA_idx@index[,,,,c('AB','C')])
range(MA_idx.all)[c('startf', 'endf')] <- c(0.5,0.6)
#MA_idx.all@index@.Data[which(is.na(MA_idx.all@index))] <- -1
# MA_idx.all <- MA_idx.all[ageRange,] # start from age 1 as age 0 is flimsy
MA_idx.all@index[,ac(2020)] <- NA

### ============================================================================
### caton as a separate index
### ============================================================================
FLR_temp <- FLQuant(  array(0,dim = c(1,dim(pil.all@catch)[2],1,1,1,1)),
                      dimnames=list(age='all',
                                    year=dimnames(pil.all@catch)$year,
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))

caton_idx <- FLIndex(name = 'caton',
                  desc = 'caton',
                  index=FLR_temp)
caton_idx@index[] <- pil.all@catch
range(caton_idx)[c('startf', 'endf')] <- c(0.2,0.7)

# region AB
FLR_temp <- FLQuant(  array(0,dim = c(1,dim(pil.all@catch)[2],1,1,1,1)),
                      dimnames=list(age='all',
                                    year=dimnames(pil.all@catch)$year,
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))

caton.AB_idx <- FLIndex(name = 'caton',
                     desc = 'caton',
                     index=FLR_temp)
caton.AB_idx@index[] <- pil.AB@catch
range(caton.AB_idx)[c('startf', 'endf')] <- c(0.2,0.7)

# region C
FLR_temp <- FLQuant(  array(0,dim = c(1,dim(pil.all@catch)[2],1,1,1,1)),
                      dimnames=list(age='all',
                                    year=dimnames(pil.all@catch)$year,
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))

caton.C_idx <- FLIndex(name = 'caton',
                        desc = 'caton',
                        index=FLR_temp)
caton.C_idx@index[] <- pil.C@catch
range(caton.C_idx)[c('startf', 'endf')] <- c(0.2,0.7)


### ============================================================================
### prepare the FLIndices
### ============================================================================

### indices region AB
pil.AB.tun <- FLIndices(MA_idx=MA_idx[,,,,'AB'],
                        MA_ssb=MA_ssb[,,,,'AB'],
                        rec_russ=rec_russ.N,
                        caton=caton.AB_idx)
pil.AB.tun   <- lapply(pil.AB.tun,function(x) {x@type <- "number"; return(x)})

# indices region C
pil.C.tun <- FLIndices(MA_idx=MA_idx[,,,,'C'],
                       MA_ssb=MA_ssb[,,,,'C'],
                       rec_russ=rec_russ.S,
                       caton=caton.C_idx)
pil.C.tun   <- lapply(pil.C.tun,function(x) {x@type <- "number"; return(x)})


pil.all.tun <- FLIndices(MA_idx=MA_idx.all,
                         MA_ssb=MA_ssb.all,
                         rec_russ=rec_russ,
                         nansen=nansen,
                         caton=caton_idx)
pil.all.tun   <- lapply(pil.all.tun,function(x) {x@type <- "number"; return(x)})

### ============================================================================
### save objects
### ============================================================================

save(pil.AB, pil.AB.tun, 
     file=file.path(dataPath,'data_pil.pelnord.AB.Rdata'))

save(pil.C, pil.C.tun, 
     file=file.path(dataPath,'data_pil.pelnord.C.Rdata'))

save(pil.all, pil.all.tun, 
     file=file.path(dataPath,'data_pil.pelnord.all.Rdata'))

