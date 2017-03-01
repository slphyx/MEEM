
# for removing NA from a data frame/list
removeNA<-function(inputDF){
  tmpDF <- inputDF

  for(n in names(inputDF)){
    txt<-paste0("tmpDF$",n,"<-","tmpDF$",n,"[!is.na(tmpDF$",n,")]")
    eval(parse(text = txt))
  }
  return(tmpDF)
}


MM_readpars <- function(parfile){

  if(!file.exists(parfile)){
    stop("Your parameter file does not exist!")
  }

  tmp <- read.csv(parfile)
  tmp <- removeNA(as.list(tmp))

  return(list(
    # Common parameters to both species
    Pmaxv=tmp$Pmaxv,
    Pmaxf=tmp$Pmaxf,
    amp=tmp$amp,
    phi=tmp$phi,
    indexhet=tmp$indexhet,
    propgd=tmp$propgd,
    delta_m=tmp$delta_m,
    bites=tmp$bites,
    prob_h=tmp$prob_h,
    demog=tmp$demog,
    eff_itn=tmp$eff_itn,
    stableveff=tmp$stableveff,
    hl_net=tmp$hl_net,
    eff_vmw=tmp$eff_vmw,
    eff_his=tmp$eff_his,
    eff_oth=tmp$eff_oth,
    elneff=tmp$elneff,
    pgd=tmp$pgd,
    sensgd=tmp$sensgd,
    mask=tmp$mask,
    # Humans
    #falciparum parameters
    gamma_m=tmp$gamma_m,
    prob_m= tmp$prob_m,
    ps=tmp$ps,
    psn=tmp$psn,
    pr=tmp$pr,
    prn=tmp$prn,
    nuc=tmp$nuc,
    nua=tmp$nua,
    nus=tmp$nus,
    psev=tmp$psev,
    pmort=tmp$pmort,
    pmortt=tmp$pmortt,
    rep_fat=tmp$rep_fat,
    tausev=tmp$tausev,
    gamma_h=tmp$gamma_h,
    zeta_a=tmp$zeta_a,
    zeta_n=tmp$zeta_n,
    chi=tmp$chi,
    # this will change if the RDT detection limit changes chi=chi_old*(mn_c-dl_RDTold)/(mn_c-dl_RDTnew)
    omega=tmp$omega, # loss of immunity
    #control
    nut=tmp$nut, # recovery rate under treatment (r_t)
    nuq=tmp$nuq, # recovery rate under quinine treatment (r_q)
    ptf=tmp$ptf, #Probability of treatment failure
    ptfc=tmp$ptfc, # Probability of being clinical after treatment failure
    ptftr=tmp$ptftr, #Probability of seeking trt if clinical, after treatment failure

    # diagnostics
    dl_RDT=tmp$dl_RDT, # standard RDT detection limit
    dl_micro=tmp$dl_micro, # micro detection limit
    dl_qPCR=tmp$dl_qPCR, # standard PCR detection limit
    # dl_inf=log10(26), # detection limit for infectiousness
    mn_n=tmp$mn_n, # mean parasiteamia for sub micro
    mn_a=tmp$mn_a, # mean parasiteamia for asym
    mn_c=tmp$mn_c, # mean parasiteamia for clinical
    mn_s=tmp$mn_s, # mean parasiteamia for severe
    sd_n=tmp$sd_n, # sd parasiteamia for sub micro
    sd_a=tmp$sd_a, # sd parasiteamia for asym
    sd_c=tmp$sd_c, # sd parasiteamia for clinical
    sd_s=tmp$sd_s, # sd parasiteamia for severe

    #vivax parameters
    vgamma_m=tmp$vgamma_m, # latent period
    vprob_h=tmp$vprob_h, # probability that a bite will result in infection (human to mosquito)
    vps=tmp$vps, # probability of clinical if non-immune
    vpsn=tmp$vpsn, # probability of sub-patent given not clin, non-immune
    vpr=tmp$vpr, # probability of clinical if immune
    vprn=tmp$vprn, # probability of sub-patent given not clin, immune
    vnuc=tmp$vnuc, # recovery from clinical symptoms untreated (r_c)
    vnua=tmp$vnua, # recovery from asym untreated (r_a)
    vnus=tmp$vnus, # recovery from severe (r_s)
    vpsev=tmp$vpsev, # probability of severe disease given clinical
    vpmort=tmp$vpmort, # proportional of all untreated severe cases that die (theta)
    vpmortt=tmp$vpmortt, # proportional of all treated severe cases that die (theta)
    vtausev=tmp$vtausev, # probability that a severe infection is treated
    vgamma_h=tmp$vgamma_h, # latent period in humans
    vzeta_a=tmp$vzeta_a, # relative infectiousness of asym to clinical
    vzeta_n=tmp$vzeta_n, # relative infectiousness of submicro to clinical
    vomega=tmp$vomega, # loss of immunity
    vprel=tmp$vprel, # probability of relapse
    vincprel=tmp$vincprel, #Probabily of relapse due to triggering
    vrel=tmp$vrel, # rate of relapse
    vph=tmp$vph, # probability of recovering with hypnozoites under ACT
    vphprim=tmp$vphprim, # probability of recovering with hypnozoites under primaquine
    vkappa=tmp$vkappa , #  hypnozoite death rate
    vnut=tmp$vnut, # recovery rate under treatment (r_t)
    vnuq=tmp$vnuq, # recovery rate under quinine treatment (r_q)
    vnup=tmp$vnup, # recovery rate under primaquine
    vptf=tmp$vptf, #Probability of treatment failure on FLT
    vptfp=tmp$vptfp, #Probability of treatment failure on Primaquine( clinical failure and adherance)
    vptfc=tmp$vptfc, # Probability of being clinical after treatment failure
    vptftr=tmp$vptftr, #Probability of seeking trt if clinical, after treatment failure

    vdl_RDT=tmp$vdl_RDT, # standard RDT detection limit
    vdl_micro=tmp$vdl_micro, # micro detection limit
    vdl_qPCR=tmp$vdl_qPCR, # standard PCR detection limit
    vmn_n=tmp$vmn_n, # mean parasiteamia for sub micro
    vmn_a=tmp$vmn_a, # mean parasiteamia for asym
    vmn_c=tmp$vmn_c, # mean parasiteamia for clinical
    vmn_s=tmp$vmn_s, # mean parasiteamia for severe
    vsd_n=tmp$vsd_n, # sd parasiteamia for sub micro
    vsd_a=tmp$vsd_a, # sd parasiteamia for asym
    vsd_c=tmp$vsd_c, # sd parasiteamia for clinical
    vsd_s=tmp$vsd_s, # sd parasiteamia for severe

    t1=tmp$t1,                  # Entanglement 1 - dual treatment switch
    t2=tmp$t2                 # Entanglement 2 - triggering relapse from Pf infection switch
  ))


}
