###
# create initial condition for Malrates
MM_GenInitOde <- function(maldata){
  alldata <- loadWorkbook(maldata)
  pvxy <- readWorksheet(alldata, sheet="Sheet1")
  pvxy <- pvxy[1:nrow(pvxy),]
  
  initcondfit<-c(0.4*pvxy[1,2],0,0.2*pvxy[1,2],0,0,0,0,0,0.4*pvxy[1,2],0,0.4*pvxy[1,2],0,0.2*pvxy[1,2],0,0,0,0,0,0.4*pvxy[1,2],0,0,0,0)
  for (n in 2:nrow(pvxy)){
    initcondfit<-c(initcondfit,c(0.4*pvxy[n,2],0,0.2*pvxy[n,2],0,0,0,0,0,0.4*pvxy[n,2],0,0.4*pvxy[n,2],0,0.2*pvxy[n,2],0,0,0,0,0,0.4*pvxy[n,2],0,0,0,0))
  }
  # all initial conditions must be integers
  initcondfit<-round(initcondfit)
  initodefit<-initcondfit
  
  #statefit <- c(initodefit,0)
  
  return(initodefit)
}

# ************************************************************************************* #
# LIKELIHOOD FUNCTION
# ************************************************************************************* #
# ************************************************************************************* #
# Adjust starting conditions to best availalbe data

MM_CalcLL<-function(parfit, odemethod = "lsoda", maldata, climatedata){
  
  # ************************************************************************************* #
  # define initial conditions
  # initcondfit<-c(0.4*pvxy[1,2],0,0.2*pvxy[1,2],0,0,0,0,0,0.4*pvxy[1,2],0,0.4*pvxy[1,2],0,0.2*pvxy[1,2],0,0,0,0,0,0.4*pvxy[1,2],0,0,0,0)
  # for (n in 2:N){
  #   initcondfit<-c(initcondfit,c(0.4*pvxy[n,2],0,0.2*pvxy[n,2],0,0,0,0,0,0.4*pvxy[n,2],0,0.4*pvxy[n,2],0,0.2*pvxy[n,2],0,0,0,0,0,0.4*pvxy[n,2],0,0,0,0))
  # }
  # # all initial conditions must be integers
  # initcondfit<-round(initcondfit)
  # initodefit<-initcondfit
  #
  # statefit <- c(initodefit,0)
  
  
  
  
  initodefit <- MM_GenInitOde(maldata = maldata)
  statefit <- c(initodefit,0)
  ti<-1
  # maldata = malaria prevalence data
  inp<-MM_Inputs(parfit, maldata = maldata, climatedata = climatedata)
  transitfit <- MM_Malrates(initodefit,inp,parfit,0,ti)
  
  # # SOLVE THE ODEs and get output
  timesfit <- seq(0, tyears, by = dtout) # Model run time
  #Solve ODE
  outodefit <- ode(y = statefit, times = timesfit, func = MM_EpiModel, parms = parfit, method  = odemethod, input=inp)
  # Compute transitions at each time step
  tranodefit<-matrix(0,nrow=length(outodefit[,1]),ncol=length(transitions))
  for (ti in 1:(tsteps+1)){
    tranodefit[ti,]<-t(MM_Malrates(outodefit[ti,2:(1+V)],inp, parfit,0,ti))
  }
  #Compute outputs
  #MM_Postproc <- function(parpro,out,tran)
  ppoutfit<-MM_Postproc(parfit,outodefit,tranodefit)
  modeltimes<-outodefit[,1]+startyear
  
  vmw_predf1_fit<-ppoutfit[,(2*N+1):(3*N)]
  vmw_predv1_fit<-ppoutfit[,(3*N+1):(4*N)]
  vmw_predmix_fit<-ppoutfit[,(4*N+1):(5*N)]
  his_predf1_fit<-ppoutfit[,(7*N+1):(8*N)]
  his_predv1_fit<-ppoutfit[,(8*N+1):(9*N)]
  his_predmix_fit<-ppoutfit[,(9*N+1):(10*N)]
  fatal_pred_fit<-ppoutfit[,(10*N+1):(11*N)]
  severe_predf1_fit<-ppoutfit[,(13*N+1):(14*N)]
  severe_predv1_fit<-ppoutfit[,(14*N+1):(15*N)]
  severe_predmix_fit<-ppoutfit[,(15*N+1):(16*N)]
  
  ##############################################
  # Compute Likelihood (only fits to available data values - ensure that missing data is NA in data files)
  
  # Poisson likelihood - Pf (vmw+his), Pv (vmw+his), Mixed inf, Severe, Fatalities
  L_hisf<- L_hisv<-L_hismix<-matrix(0,nrow=sum(his_inc_t%in%modeltimes),ncol=N)
  
  for (i in 1:N) {
    L_hisf[,i]<-dpois(his_incf[his_inc_t%in%modeltimes,i],his_predf1_fit[modeltimes[1:(length(modeltimes)-1)]%in%his_inc_t,i],log=TRUE)
    L_hisv[,i]<-dpois(his_incv[his_inc_t%in%modeltimes,i],his_predv1_fit[modeltimes[1:(length(modeltimes)-1)]%in%his_inc_t,i],log=TRUE)
    L_hismix[,i]<-dpois(his_incmix[his_inc_t%in%modeltimes,i],his_predmix_fit[modeltimes[1:(length(modeltimes)-1)]%in%his_inc_t,i],log=TRUE)
  }
  LL_his<-sum(L_hisf[which(L_hisf!=-Inf & L_hisf!="NA")])+sum(L_hisv[which(L_hisv!=-Inf & L_hisv!="NA")])+sum(L_hismix[which(L_hismix!=-Inf & L_hismix!="NA")])
  
  L_vmwf<- L_vmwv<-L_vmwmix<-matrix(0,nrow=sum(vmw_inc_t%in%modeltimes),ncol=N)
  for (i in 1:N) {
    L_vmwf[,i]<-dpois(vmw_incf[vmw_inc_t%in%modeltimes,i],vmw_predf1_fit[modeltimes[1:(length(modeltimes)-1)]%in%vmw_inc_t,i],log=TRUE)
    L_vmwv[,i]<-dpois(vmw_incv[vmw_inc_t%in%modeltimes,i],vmw_predv1_fit[modeltimes[1:(length(modeltimes)-1)]%in%vmw_inc_t,i],log=TRUE)
    L_vmwmix[,i]<-dpois(vmw_incmix[vmw_inc_t%in%modeltimes,i],vmw_predmix_fit[modeltimes[1:(length(modeltimes)-1)]%in%vmw_inc_t,i],log=TRUE)
  }
  LL_vmw<-sum(L_vmwf[which(L_vmwf!=-Inf & L_vmwf!="NA")])+sum(L_vmwv[which(L_vmwv!=-Inf & L_vmwv!="NA")])+sum(L_vmwmix[which(L_vmwmix!=-Inf & L_vmwmix!="NA")])
  
  #Fatalities
  L_fatal=matrix(0,nrow=sum(fatal_inc_t%in%modeltimes),ncol=N)
  for (i in 1:N) {
    L_fatal[,i]<-dpois(fatal_inc[fatal_inc_t%in%modeltimes,i],parfit$rep_fat*fatal_pred_fit[modeltimes[1:(length(modeltimes)-1)]%in%fatal_inc_t,i],log=TRUE)}
  LL_fatal<-sum(L_fatal[which(L_fatal!=-Inf& L_fatal!="NA")])
  
  #Severe Infections
  L_severe=matrix(0,nrow=sum(severe_inc_t%in%modeltimes),ncol=N)
  severe_fitdata<-severe_predf1_fit+severe_predv1_fit+severe_predmix_fit
  for (i in 1:N) {
    L_severe[,i]<-dpois(severe_inc[severe_inc_t%in%modeltimes,i],severe_fitdata[modeltimes[1:(length(modeltimes)-1)]%in%severe_inc_t,i],log=TRUE)}
  LL_severe<-sum(L_severe[which(L_severe!=-Inf& L_severe!="NA")])
  
  
  LL<-LL_his+LL_vmw+LL_fatal+LL_severe
  NLL<--1*LL
  return(list(NLL,his_predf1_fit,his_predv1_fit,his_predmix_fit,vmw_predf1_fit,vmw_predv1_fit,vmw_predmix_fit,fatal_pred_fit,severe_fitdata ,outodefit))
}


###
### for running the model
## taken from MM_CalcLL but the likielihood was removed
MM_RunMod <- function(parfit, odemethod = "lsoda", maldata, climatedata, parfile=NULL, parallel=FALSE){
  
  
  initodefit <- MM_GenInitOde(maldata = maldata)
  statefit <- c(initodefit,0)
  
  
  if(!is.null(parfile)){
    parfit <- MM_readpars(parfile)
  }
  
  if(parallel){
    ncores<-parallel::detectCores()
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
  }
  
  # maldata = malaria prevalence data
  inp<-MM_Inputs(parfit, maldata = maldata, climatedata = climatedata)
  ti<-1
  transitfit <- MM_Malrates(initodefit,inp,parfit,0,ti)

  tyears<-inp$tyears
  dtout<-inp$dtout
  transitionslength<-inp$transitionslength
  
  # # SOLVE THE ODEs and get output
  timesfit <- seq(0, tyears, by = dtout) # Model run time
  
  #Solve ODE
  outodefit <- ode(y = statefit, times = timesfit, func = MM_EpiModel, parms = parfit, method  = odemethod, input=inp)

  # Compute transitions at each time step
  tranodefit<-matrix(0,nrow=length(outodefit[,1]),ncol=transitionslength)
  
  tsteps<-inp$tsteps
  V<-inp$V
  if(parallel){
    tranodefit<-foreach(ti = 1:(tsteps+1), .combine = rbind )%dopar%{
      t(MM_Malrates(outodefit[ti,2:(1+V)],inp, parfit,0,ti))
      
    }
  }else{
    for(ti in 1:(tsteps+1)){
      tranodefit[ti,]<-t(MM_Malrates(outodefit[ti,2:(1+V)],inp, parfit,0,ti))
    }
  }
  
  
  N<-inp$N
  
  #Compute outputs
  #MM_Postproc <- function(parpro,out,tran)
  ppoutfit<-MM_Postproc(parfit,outodefit,tranodefit,npatch=N)
  startyear<-inp$startyear
  modeltimes<-outodefit[,1]+startyear

  vmw_predf1_fit<-ppoutfit[,(2*N+1):(3*N)]
  vmw_predv1_fit<-ppoutfit[,(3*N+1):(4*N)]
  vmw_predmix_fit<-ppoutfit[,(4*N+1):(5*N)]
  his_predf1_fit<-ppoutfit[,(7*N+1):(8*N)]
  his_predv1_fit<-ppoutfit[,(8*N+1):(9*N)]
  his_predmix_fit<-ppoutfit[,(9*N+1):(10*N)]
  fatal_pred_fit<-ppoutfit[,(10*N+1):(11*N)]
  severe_predf1_fit<-ppoutfit[,(13*N+1):(14*N)]
  severe_predv1_fit<-ppoutfit[,(14*N+1):(15*N)]
  severe_predmix_fit<-ppoutfit[,(15*N+1):(16*N)]
  
  return(list(vmw_predf1_fit=vmw_predf1_fit,
              vmw_predv1_fit=vmw_predv1_fit,
              vmw_predmix_fit=vmw_predmix_fit,
              his_predf1_fit=his_predf1_fit, 
              his_predv1_fit=his_predv1_fit, 
              his_predmix_fit=his_predmix_fit,
              fatal_pred_fit=fatal_pred_fit,
              severe_predf1_fit=severe_predf1_fit,
              severe_predv1_fit=severe_predv1_fit, 
              severe_predmix_fit=severe_predmix_fit
  ))
}



