# ************************************************************************************* #
# Function to calculate inputs to transition rates
# ************************************************************************************* #
# define the error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## create an input for MM_Malrates function
MM_Inputs<-function(parmal, maldata = NULL, climatedata=NULL){

  if(!is.null(maldata) && !is.null(climatedata)){

    # ************************************************************************************* #
    # import data
    # ************************************************************************************* #

    #alldata = loadWorkbook(system.file("extdata", "Cambodia_PfPv_provincial_data.xlsx", package="MEEM"))
    #climate = loadWorkbook(system.file("extdata", "climate.xlsx", package="MEEM"))

    if(file.exists(maldata)){
      alldata = loadWorkbook(maldata);
    }else{
      stop(paste("I can't find your data file ",maldata,"."));
    }
    if(file.exists(climatedata)){
      climate = loadWorkbook(climatedata);
    }else{
      stop(paste("I can't find your data file ",climatedata,"."));
    }

    # population and villages
    pvxy = readWorksheet(alldata, sheet="Sheet1")

    N<-nrow(pvxy)   # number of patches
    B<-23   # number of variables per patch
    A<-261  # number of transitions per patch
    V<-N*B # total number of variables
    L<-N*A #total number of transitions
    startyear=1995 # starting year of simulation
    tyears<-20 # total years of simulation
    dtout<-1/12 # output timestep
    tsteps<-round(tyears/dtout) # number of time steps
    time<-startyear+seq(0,tyears,dtout) # time vector


    pvxy<-pvxy[1:N,]
    yrprim<-pvxy[,6] #year primaquine adopted
    g6pDd<-pvxy[,7] # G6PDd def proportion


    # VMW villages over time
    vmwdat = readWorksheet(alldata, sheet="Sheet2")
    vmw_vil<-vmwdat[,(1:N)+2]
    vmw_time<-vmwdat[,1]+vmwdat[,2]/12 # end of the month
    cov_vmw_dat<-matrix(0,nrow=length(vmw_time),ncol=N)
    for (n in 1:N){
      cov_vmw_dat[,n]<-100*vmw_vil[,n]/pvxy[n,3]
    }

    # ITN distributions over time
    itndat = readWorksheet(alldata, sheet="Sheet3")
    itn_dis<-itndat[,(1:N)+2]
    itn_time<-itndat[,1]+itndat[,2]/12 # end of the month
    cov_itn_dat<-matrix(0,nrow=length(itn_time),ncol=N)
    for (n in 1:N){
      cov_itn_dat[,n]<-100*itn_dis[,n]/pvxy[n,2]
    }

    # VMW cases over time
    vmw_casesf = readWorksheet(alldata, sheet="Sheet8")
    vmw_incf<-vmw_casesf[,(1:N)+2]
    vmw_casesv = readWorksheet(alldata, sheet="Sheet10")
    vmw_incv<-vmw_casesv[,(1:N)+2]
    vmw_casesmix = readWorksheet(alldata, sheet="Sheet12")
    vmw_incmix<-vmw_casesmix[,(1:N)+2]
    vmw_inc_t<-vmw_casesf[,1]+vmw_casesf[,2]/12 # end of the month


    # HIS cases over time
    his_casesf = readWorksheet(alldata, sheet="Sheet9")
    his_incf<-his_casesf[,(1:N)+2]
    his_casesv = readWorksheet(alldata, sheet="Sheet11")
    his_incv<-his_casesv[,(1:N)+2]
    his_casesmix = readWorksheet(alldata, sheet="Sheet13")
    his_incmix<-his_casesmix[,(1:N)+2]
    his_inc_t<-his_casesf[,1]+his_casesf[,2]/12 # end of the month

    # severe cases over time
    severe_cases = readWorksheet(alldata, sheet="Sheet6")
    severe_inc<-severe_cases[,(1:N)+2]
    severe_inc_t<-severe_cases[,1]+severe_cases[,2]/12 # end of the month

    # fatal cases over time
    fatal_cases = readWorksheet(alldata, sheet="Sheet7")
    fatal_inc<-fatal_cases[,(1:N)+2]
    fatal_inc_t<-fatal_cases[,1]+fatal_cases[,2]/12 # end of the month


    # Best El Nino index over time
    elnino_index = readWorksheet(climate, sheet="Sheet1")
    elnino<-elnino_index[,7]
    elnino_t<-elnino_index[,1]+elnino_index[,2]/12 # end of the month #b sure that the start years match


    # ************************************************************************************* #
    # define variables
    # ************************************************************************************* #
    # FALCIPARUM
    # 1=S: uninfected non-immune
    # 2=In: infected submicro
    # 3=Ia: infected asymp
    # 4=Ic: infected clinical
    # 5=Is: infected severe
    # 6=To: treated unrecorded
    # 7=Tv: treated VMW
    # 8=Th: treated HIS
    # 9=R: uninfected immune
    # 10=H: uninfected and HRP2 positive

    #VIVAX
    # 11=S: uninfected non-immune
    # 12=In: infected submicro
    # 13=Ia: infected asymp
    # 14=Ic: infected clinical
    # 15=Is: infected severe
    # 16=To: treated unrecorded
    # 17=Tv: treated VMW without G6PDd
    # 18=Th: treated HIS without G6PDd
    # 19=R: uninfected immune with no hypnozoites
    # 20=L: uninfected immune with hypnozoites
    # 21=Togd: treated unrecorded with G6PDd
    # 22=Tvgd: treated VMW with G6PDd
    # 23=Thgd: treated HIS with G6PDd

    falpop<-1:10
    vivpop<-11:23

    # ************************************************************************************* #
    # define indices
    # ************************************************************************************* #
    varind<-matrix(0,nrow=B,ncol=N)
    traind<-matrix(0,nrow=A,ncol=N)
    for (n in 1:N){
      for (b in 1:B){
        varind[b,n]<-(n-1)*B+b
      }
      for (a in 1:A){
        traind[a,n]<-(n-1)*A+a
      }
    }

    #####
  }

  if((is.null(maldata) && is.null(climatedata)) && !is.data.frame(pvxy)){
    stop("Please load all parameters required by MEEM from MM_Init() if you would like to use our default data files.")
  }

  # sensitivities
  sens_c_micro<-1-0.5*(1+erf((parmal$dl_micro-parmal$mn_c)/(parmal$sd_c*(2^0.5))))
  sens_c_RDT<-1-0.5*(1+erf((parmal$dl_RDT-parmal$mn_c)/(parmal$sd_c*(2^0.5))))
  vsens_c_micro<-1-0.5*(1+erf((parmal$vdl_micro-parmal$vmn_c)/(parmal$vsd_c*(2^0.5))))
  vsens_c_RDT<-1-0.5*(1+erf((parmal$vdl_RDT-parmal$vmn_c)/(parmal$vsd_c*(2^0.5))))

  sens_vmw<-sens_c_RDT  # test used by VMW
  sens_his<-sens_c_micro # test used by HIS
  vsens_vmw<-vsens_c_RDT  # test used by VMW
  vsens_his<-vsens_c_micro # test used by HIS
  sens_oth<-1 # test used by other

  # detection limits
  dl_0<-parmal$dl_qPCR
  vdl_0<-parmal$vdl_qPCR

  # durations of infection
  nun<-1/((1/parmal$nua)*(parmal$mn_n-dl_0)/(parmal$mn_a-parmal$mn_n)) #(r_n)
  vnun<-1/((1/parmal$vnua)*(parmal$vmn_n-vdl_0)/(parmal$vmn_a-parmal$vmn_n)) #(r_n)

  # ************************************************************************************* #
  # set up historical coverage over time using data and parameter values
  # ************************************************************************************* #
  # VMW coverage
  vmw_time<-vmwdat[,1]+vmwdat[,2]/12 # end of the month
  cov_vmw<-matrix(0,nrow=length(vmw_time),ncol=N)
  for (n in 1:N){
    cov_vmw[,n]<-vmw_vil[,n]/pvxy[n,3]
  }

  if (vmw_time[1]>startyear){
    vmw_time<-c(startyear,vmw_time[1]-1/12,vmw_time)
    cov_vmw<-rbind(matrix(0,nrow=2,ncol=N),cov_vmw)
  }

  fvmw<-length(vmw_time)

  if (vmw_time[fvmw]<(startyear+tyears)){
    vmw_time<-c(vmw_time,vmw_time[fvmw]+1/12,(startyear+tyears)+1/12)
    cov_vmw<-rbind(cov_vmw,cov_vmw[fvmw,],cov_vmw[fvmw,])
  }

  vmw_eff<-matrix(0,nrow=length(vmw_time),ncol=N)
  vmw_dec<-c(parmal$stableveff+exp(log(0.02)/5),parmal$stableveff+exp(log(0.02)/5*(seq(1,12,1/12)))) # exponentially decreasing function for effectiveness (0.85 to 0.4 within 5 years)

  cov_st<-c() #cover start
  for (n in 1:N){
    st<- cbind(tapply(cov_vmw[,n], vmw_time,  function(a) which(a>0)[1]))
    cov_st<-c(cov_st,sum(is.na(st)))
  }

  for (n in 1:N){
    vmw_eff[cov_st[n]:length(vmw_time),n]<-vmw_dec[1:(length(vmw_time)-cov_st[n]+1)]
  }


  #itn coverage
  cov_itn<-matrix(0,nrow=length(itn_time),ncol=N)

  for (n in 1:N){
    cov_itn[,n]<-itn_dis[,n]/pvxy[n,2]
  }

  fitn<-length(itn_time)

  if (itn_time[fitn]<(startyear+tyears)){
    itn_time<-c(itn_time,itn_time[fitn]+1/12,(startyear+tyears)+1/12)
    cov_itn<-rbind(cov_itn,cov_itn[fitn,],cov_itn[fitn,])
  }

  c_itn<-matrix(0,nrow=length(itn_time),ncol=N)
  c_itn[1,]<-cov_itn[1,]
  for (tt in 2:length(itn_time))
  {
    dt_itn<-itn_time[tt]-itn_time[tt-1]
    c_itn[tt,]<-cov_itn[tt,]+0.5*c_itn[tt-1,]*exp(-dt_itn/parmal$hl_net)
  }

  c_itn<-parmal$eff_itn*c_itn

  if (itn_time[1]>startyear){
    itn_time<-c(startyear,itn_time[1]-1/12,itn_time)
    c_itn<-rbind(matrix(0,nrow=2,ncol=N),c_itn)
  }


  #set up El nino effect
  mineln<-rep(1,N)
  eln_inp<-runmed(elnino_index[,7], parmal$elneff)
  eln_t<-elnino_index[,1]+elnino_index[,2]/12 # end of the month #b sure that the start years match
  # for (n in 1:N) {
  #   mineln[n]<-abs(min(eln_inp[(494-length(time)+1):494]*parmal$amp[n]*cos(2*pi*(time-parmal$phi[n]))))
  # }


  # ************************************************************************************* #
  # set up spatial connectivity
  # ************************************************************************************* #
  dist<-matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    for (j in 1:N){
      dist[i,j]<-0.001*((pvxy[i,4]-pvxy[j,4])^2+(pvxy[i,5]-pvxy[j,5])^2)^0.5
    }
  }

  # matrix of connectivities
  connect<-matrix(0,nrow=N,ncol=N)
  het<-10^(parmal$indexhet)
  # distance only
  for (i in 1:N){
    for (j in 1:N){
      connect[i,j]<-(1/(1+het[i]*dist[i,j]))/(sum(1/(1+het[i]*dist[i,])))
    }
  }
  #weight distance probabilities by population
  connect2<-connect
  for (i in 1:N){
    for (j in 1:N){
      connect2[i,j]<-pvxy[j,2]*connect[i,j]/sum(pvxy[,2])
    }
  }
  #scaling pop-dis weights back to sum to 1 to indicate probability of movement
  connect3<-connect2
  for (i in 1:N){
    for (j in 1:N){
      connect3[i,j]<-connect2[i,j]/sum(connect2[i,])
    }
  }

  return(list(nun=nun,vnun=vnun, connect=connect, connect3=connect3,cov_vmw=cov_vmw,vmw_eff=vmw_eff, vmw_time=vmw_time, c_itn=c_itn, itn_time=itn_time, sens_vmw=sens_vmw, sens_his=sens_his, sens_oth=sens_oth, eln_inp=eln_inp, eln_t=eln_t, mineln=mineln))
}
#######################################################################



