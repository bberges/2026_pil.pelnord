config_sf_pil.all_cecaf2025 <- function(  stk,
                                          tun,
                                          pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  nansenRow     <- grep("nansen",rownames(sam.ctrl@obs.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))

  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)

  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(nansenRow,MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat

  fleets <- sam.ctrl@fleets[]
  fleets[c(nansenRow,MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets

  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,4)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  sam.ctrl@obs.vars[idxRow,ac(2:5)]   <- c(101,103,104,104)
  # sam.ctrl@obs.vars[idxRow,ac(0:5)]   <- c(101,103,103,104,105,105)
  
  # # # #Catchabilities of the surveys
  sam.ctrl@catchabilities[idxRow,ac(2:5)]   <- c(102,103,104,104)
  # sam.ctrl@catchabilities[idxRow,ac(0:5)]   <- c(101,102,103,104,105,105)
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  # noScaledYearsFit <- 10
  # sam.ctrl@scaleNoYears[] <- as.integer(noScaledYearsFit)
  # sam.ctrl@scaleYears <-   as.numeric((range(stk)['maxyear'] - noScaledYearsFit + 1):range(stk)['maxyear'])
  # 
  # keyParScaledYA =  matrix(data = rep(c(0), 
  #                                     each = sam.ctrl@scaleNoYears * dims(stk)$age),
  #                          nrow = sam.ctrl@scaleNoYears) # one param all years
  # 
  # scalePars        <- matrix(as.integer(0),nrow=sam.ctrl@scaleNoYears,ncol=dims(stk)$age,
  #                            dimnames=list(years=sam.ctrl@scaleYears,age=dimnames(stk@catch.n)$age))
  #scalePars[1,] <- as.integer(0)
  # scalePars[2,] <- as.integer(1)
  # scalePars[3,] <- as.integer(1)
  # scalePars[4,] <- as.integer(1)
  
  #sam.ctrl@scalePars <- scalePars
  # sam.ctrl@scalePars
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}

config_sf_pil.all_cecaf2025_v2 <- function(  # confirguation with recruitment at age 1
    stk,tun,pg){
  sam.ctrl                                    <- FLSAM.control(stk,
                                                               tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  nansenRow     <- grep("nansen",rownames(sam.ctrl@obs.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  
  # set recruitment function to BH
  sam.ctrl@srr[] <- as.integer(3)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(nansenRow,MA.ssbRow)] <- 5
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(nansenRow,MA.ssbRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(2:pg)))
  # 
  # #All fishing mortality states are free except
  # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,2,3,3,3)#seq(dim(stk)[1])
  # 
  # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,4,4,5)
  sam.ctrl@obs.vars[idxRow,ac(1:5)]   <- c(101,102,104,104,105)
  # sam.ctrl@obs.vars[recRow,ac(1)]   <- c(200)
  # sam.ctrl@obs.vars[nansenRow,ac(1)]       <- c(300)
  # sam.ctrl@obs.vars[MA.ssbRow,ac(1)]       <- c(400)

  # #Catchabilities of the surveys
  sam.ctrl@catchabilities[idxRow,ac(1:5)]   <- c(101,102,103,104,104)
  # sam.ctrl@catchabilities[recRow,ac(1)]     <- c(200)
  # sam.ctrl@catchabilities[nansenRow,ac(1)]       <- c(300)
  # sam.ctrl@catchabilities[MA.ssbRow,ac(1)]       <- c(400)
  # 
  # #Add correlation correction for Q3 survey, not for HERAS
  idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  sam.ctrl@cor.obs[idx,1:4]                   <- c(rep(101,4))
  sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}

config_sf_pil.all_cecaf2025_v3 <- function(  stk,
                                          tun,
                                          pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  nansenRow     <- grep("nansen",rownames(sam.ctrl@obs.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))
  
  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(nansenRow,MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(nansenRow,MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,5)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  # sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  # sam.ctrl@obs.vars[idxRow,ac(2:5)]   <- c(101,103,104,104)
  
  # # # #Catchabilities of the surveys
  # sam.ctrl@catchabilities[idxRow,ac(2:5)]   <- c(102,103,104,104)
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  # noScaledYearsFit <- 10
  # sam.ctrl@scaleNoYears[] <- as.integer(noScaledYearsFit)
  # sam.ctrl@scaleYears <-   as.numeric((range(stk)['maxyear'] - noScaledYearsFit + 1):range(stk)['maxyear'])
  # 
  # keyParScaledYA =  matrix(data = rep(c(0), 
  #                                     each = sam.ctrl@scaleNoYears * dims(stk)$age),
  #                          nrow = sam.ctrl@scaleNoYears) # one param all years
  # 
  # scalePars        <- matrix(as.integer(0),nrow=sam.ctrl@scaleNoYears,ncol=dims(stk)$age,
  #                            dimnames=list(years=sam.ctrl@scaleYears,age=dimnames(stk@catch.n)$age))
  #scalePars[1,] <- as.integer(0)
  # scalePars[2,] <- as.integer(1)
  # scalePars[3,] <- as.integer(1)
  # scalePars[4,] <- as.integer(1)
  
  #sam.ctrl@scalePars <- scalePars
  # sam.ctrl@scalePars
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}

config_sf_pil.all_cecaf2025_v4 <- function(  stk,
                                          tun,
                                          pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  nansenRow     <- grep("nansen",rownames(sam.ctrl@obs.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))
  
  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(nansenRow,MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(nansenRow,MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,4)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  sam.ctrl@obs.vars[idxRow,ac(0:5)]   <- c(101,102,104,104,105,105)
  
  # # # #Catchabilities of the surveys
  sam.ctrl@catchabilities[idxRow,ac(0:5)]   <- c(101,102,104,104,105,105)
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  # noScaledYearsFit <- 10
  # sam.ctrl@scaleNoYears[] <- as.integer(noScaledYearsFit)
  # sam.ctrl@scaleYears <-   as.numeric((range(stk)['maxyear'] - noScaledYearsFit + 1):range(stk)['maxyear'])
  # 
  # keyParScaledYA =  matrix(data = rep(c(0), 
  #                                     each = sam.ctrl@scaleNoYears * dims(stk)$age),
  #                          nrow = sam.ctrl@scaleNoYears) # one param all years
  # 
  # scalePars        <- matrix(as.integer(0),nrow=sam.ctrl@scaleNoYears,ncol=dims(stk)$age,
  #                            dimnames=list(years=sam.ctrl@scaleYears,age=dimnames(stk@catch.n)$age))
  #scalePars[1,] <- as.integer(0)
  # scalePars[2,] <- as.integer(1)
  # scalePars[3,] <- as.integer(1)
  # scalePars[4,] <- as.integer(1)
  
  #sam.ctrl@scalePars <- scalePars
  # sam.ctrl@scalePars
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}


config_sf_pil.AB_cecaf2025 <- function(  
    stk,tun,pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))
  
  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,4)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  sam.ctrl@obs.vars[idxRow,ac(2:5)]   <- c(101,103,104,104)
  
  # # # #Catchabilities of the surveys
  sam.ctrl@catchabilities[idxRow,ac(2:5)]   <- c(102,103,104,104)
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}


config_sf_pil.AB_cecaf2025_v2 <- function(  
    stk,tun,pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))
  
  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,4)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  sam.ctrl@obs.vars[idxRow,ac(2:5)]   <- c(101,103,104,104)
  
  # # # #Catchabilities of the surveys
  sam.ctrl@catchabilities[idxRow,ac(2:5)]   <- c(102,103,104,104)
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}


config_sf_pil.C_cecaf2025 <- function(  
    stk,tun,pg){
  sam.ctrl      <- FLSAM.control(stk,
                                 tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  MA.ssbRow     <- grep("MA_ssb",rownames(sam.ctrl@obs.vars))
  recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  idxRow        <- grep("MA_idx",rownames(sam.ctrl@obs.vars))
  catonRow         <- grep("caton",rownames(sam.ctrl@f.vars))
  
  # set recruitment function
  sam.ctrl@srr[] <- as.integer(0)
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(MA.ssbRow)] <- 5
  biomassTreat[c(catonRow)] <- 1
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(MA.ssbRow)] <- 3
  fleets[c(catonRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # # #All fishing mortality states are free except
  # # #oldest ages to ensure stablity
  sam.ctrl@states[catchRow,]                  <- c(0,1,2,3,4,4)#seq(dim(stk)[1])
  # 
  # # # Group observation variances of catches to ensure stability
  sam.ctrl@obs.vars[catchRow,]        <- c(1,2,2,4,4,6)
  # sam.ctrl@obs.vars[idxRow,ac(2:5)]   <- c(101,103,104,104)
  sam.ctrl@obs.vars[idxRow,ac(0:5)]   <- c(101,102,104,104,105,105)
  
  # # # #Catchabilities of the surveys
  # sam.ctrl@catchabilities[idxRow,ac(2:5)]   <- c(102,103,104,104)
  sam.ctrl@catchabilities[idxRow,ac(0:5)]   <- c(101,102,104,104,105,105)
  
  
  # # # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,3:5]                   <- c(rep(101,3))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  
  #Variance of F random walks
  #sam.ctrl@f.vars[catchRow,]    <- c(0,1,1,1,1)
  #correlation structure among ages in F random walks
  sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "sardine cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}
