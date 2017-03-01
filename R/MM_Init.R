## parameter list created by Sheetal
MM_SheetalPars <- function(){
    N <- 24
    g6pDd <- rep(0.143,N)
# ************************************************************************************* #
# Set the parameters (vivax-specific parameters prefixed with v)
# ************************************************************************************* #
    return(list(
      # Common parameters to both species
      Pmaxv=40000*seq_len(N)/seq_len(N), # the maximum populations of mosquitoes in each patch
      Pmaxf=40000*seq_len(N)/seq_len(N), # the maximum populations of mosquitoes in each patch
      amp=rep(1,N), # amplitude of seasonal forcing
      phi=rep(0.75,N), # phase angle of seasonal forcing
      indexhet=rep(0.5,N), # parameter for connectivity about 1
      propgd=g6pDd, #Proportion of patch that is G6PDdef
      delta_m=365.25/14, # death rate of mosquitoes
      bites=365.25/3, # biting rate
      prob_h=0.5, # probability that a bite will result in infection (mosquito to human)
      demog=1/50, # birth/death rate (mu)
      eff_itn=0.15, #efficacy of bednets
      stableveff=0.4, # steady state of vmw_effective
      hl_net=1.5, # half-life of bednets
      eff_vmw=0.6, # efficacy of VMW
      eff_his=0.25, # efficacy of HIS
      eff_oth=0.1, # efficacy of other systems
      elneff=21, # smoothing effect of el nino patterns.
      pgd=0.3, #probability of clinical if G6PDdef
      sensgd=0.97, #Sensitivity of test to detect G6PDdef
      mask=0.05, #probability of vivax mix cases masked as falciparum
      # Humans
      #falciparum parameters
      gamma_m=365.25/10, # latent period
      prob_m=0.5, # probability that a bite will result in infection (human to mosquito)
      ps=0.9, # probability of clinical if non-immune
      psn=0.1, # probability of sub-patent given not clin, non-immune
      pr=0.1, # probability of clinical if immune
      prn=0.5, # probability of sub-patent given not clin, immune
      nuc=365.25/10, # recovery from clinical symptoms untreated (r_c)
      nua=365.25/130, # recovery from asym untreated (r_a)
      nus=365.25/5, # recovery from severe (r_s)
      psev=0.03, # probability of severe disease given clinical
      pmort=0.2, # proportional of all untreated severe cases that die (theta1)
      pmortt=0.15, # proportional of all treated severe cases that die (theta2)
      rep_fat=0.6, # proportion of fatalities reported
      tausev=0.8, # probability that a severe infection is treated
      gamma_h=365.25/21, # latent period in humans
      zeta_a=12.6/27, # relative infectiousness of asym to clinical
      zeta_n=3.9/27, # relative infectiousness of submicro to clinical
      chi=365.25/28, # rate of loss of detectable HRP2 by RDT
      # this will change if the RDT detection limit changes chi=chi_old*(mn_c-dl_RDTold)/(mn_c-dl_RDTnew)
      omega=1, # loss of immunity
      #control
      nut=365.25/3, # recovery rate under treatment (r_t)
      nuq=365.25/6, # recovery rate under quinine treatment (r_q)
      ptf=0.05, #Probability of treatment failure
      ptfc=0.75, # Probability of being clinical after treatment failure
      ptftr=0.27, #Probability of seeking trt if clinical, after treatment failure

      # diagnostics
      dl_RDT=log10(200), # standard RDT detection limit
      dl_micro=log10(100), # micro detection limit
      dl_qPCR=log10(2), # standard PCR detection limit
      # dl_inf=log10(26), # detection limit for infectiousness
      mn_n=log10(5), # mean parasiteamia for sub micro
      mn_a=log10(1000), # mean parasiteamia for asym
      mn_c=log10(25000), # mean parasiteamia for clinical
      mn_s=log10(350000), # mean parasiteamia for severe
      sd_n=0.75, # sd parasiteamia for sub micro
      sd_a=1.5, # sd parasiteamia for asym
      sd_c=1.3, # sd parasiteamia for clinical
      sd_s=0.26, # sd parasiteamia for severe

      #vivax parameters
      vgamma_m=365.25/12, # latent period
      vprob_h=0.23, # probability that a bite will result in infection (human to mosquito)
      vps=0.9, # probability of clinical if non-immune
      vpsn=0.1, # probability of sub-patent given not clin, non-immune
      vpr=0.1, # probability of clinical if immune
      vprn=0.17, # probability of sub-patent given not clin, immune
      vnuc=365.25/20, # recovery from clinical symptoms untreated (r_c)
      vnua=365.25/365.25, # recovery from asym untreated (r_a)
      vnus=365.25/3, # recovery from severe (r_s)
      vpsev=0.03, # probability of severe disease given clinical
      vpmort=0.2, # proportional of all untreated severe cases that die (theta)
      vpmortt=0.02, # proportional of all treated severe cases that die (theta)
      vtausev=0.8, # probability that a severe infection is treated
      vgamma_h=365.25/17, # latent period in humans
      vzeta_a=1, # relative infectiousness of asym to clinical
      vzeta_n=1, # relative infectiousness of submicro to clinical
      vomega=1, # loss of immunity
      vprel=0.25, # probability of relapse
      vincprel=0.30, #Probabily of relapse due to triggering
      vrel=365.25/100, # rate of relapse
      vph=0.68 , # probability of recovering with hypnozoites under ACT
      vphprim=0.13, # probability of recovering with hypnozoites under primaquine
      vkappa=365.25/400 , #  hypnozoite death rate
      vnut=365.25/3, # recovery rate under treatment (r_t)
      vnuq=365.25/6, # recovery rate under quinine treatment (r_q)
      vnup=365/14, # recovery rate under primaquine
      vptf=0.05, #Probability of treatment failure on FLT
      vptfp=1-(0.9*0.85), #Probability of treatment failure on Primaquine( clinical failure and adherance)
      vptfc=0.5, # Probability of being clinical after treatment failure
      vptftr=0.27, #Probability of seeking trt if clinical, after treatment failure

      vdl_RDT=log10(200), # standard RDT detection limit
      vdl_micro=log10(100), # micro detection limit
      vdl_qPCR=log10(2), # standard PCR detection limit
      vmn_n=log10(5), # mean parasiteamia for sub micro
      vmn_a=log10(750), # mean parasiteamia for asym
      vmn_c=log10(5000), # mean parasiteamia for clinical
      vmn_s=log10(20000), # mean parasiteamia for severe
      vsd_n=0.75, # sd parasiteamia for sub micro
      vsd_a=1.5, # sd parasiteamia for asym
      vsd_c=1.3, # sd parasiteamia for clinical
      vsd_s=log10(9900), # sd parasiteamia for severe

      t1=1,                  # Entanglement 1 - dual treatment switch
      t2=1                 # Entanglement 2 - triggering relapse from Pf infection switch
    ))
}

## temporary function for creating all parameters needed by Sheetal code
MM_Init<-function(parmal=NULL, maldata=NULL,climatedata=NULL){
  if(is.null(maldata) && is.null(climatedata) && is.null(parmal)){
    fpath <- system.file("extdata", "init.R", package="MEEM")
    cat("Creating some parameters required by MEEM...\n")
    source(fpath)
  }

  if(!is.null(maldata) && !is.null(climatedata) && !is.null(parmal)){
      return(MM_Inputs(parmal, maldata = maldata, climatedata = climatedata))
  }
}


