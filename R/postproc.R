

# POST PROCESSING function
MM_Postproc <- function(parpro,out,tran, npatch) {
  
  N<-npatch
  
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
  
  
  
  # sensitivities
  sens_n_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_n_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_n_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_a_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_a_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_a_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_c_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_c_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_c_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_s_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_s_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_s_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_H_micro<-0
  sens_H_RDT<-1
  sens_H_qPCR<-0
  vsens_n_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_n_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_n_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_a_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_a_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_a_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_c_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_c_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_c_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_s_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))
  vsens_s_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))
  vsens_s_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))

  sens_vmw<-sens_c_RDT  # test used by VMW
  sens_his<-sens_c_micro # test used by HIS
  vsens_vmw<-vsens_c_RDT  # test used by VMW
  vsens_his<-vsens_c_micro # test used by HIS
  sens_oth<-1 # test used by other


  # ************************************************************************************* #
  # for outputting the  time series for each patch
  # ************************************************************************************* #

  # VMW outputs
  vmw_predv<-vmw_predf<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    vmw_predf[,n]<-rowSums(tran[,c(traind[c(16,29,34,40,47,56),n])])/12
    vmw_predv[,n]<-rowSums(tran[,c(traind[c(78,81,105,108,112,115,122,125,132,135,141,144,153,162,212,215,222,231,240,249,258),n])])/12
  }

  vmw_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    vmw_predmix[,n]<-vmw_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + vmw_predf[,n]*(0.5*out[,1+(13+(n-1)*B)]+out[,1+(14+(n-1)*B)]+out[,1+(15+(n-1)*B)]+out[,1+(16+(n-1)*B)]+out[,1+(18+(n-1)*B)]+out[,1+(21+(n-1)*B)]+out[,1+(23+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }

  vmw_predv1<-vmw_predf1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    vmw_predf1[,n]<-vmw_predf[,n]-vmw_predmix[,n]
    vmw_predv1[,n]<-vmw_predv[,n]-vmw_predmix[,n]
  }

  #HIS outputs
  his_predv<-his_predf<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predf[,n]<-rowSums(tran[,c(traind[c(17,30,35,41,48,53,59),n])])/12
    his_predv[,n]<-rowSums(tran[,c(traind[c(79,82,106,109,113,116,123,126,133,136,142,145,150,156,159,165,213,216,219,225,228,234,237,243,246,252,255,261),n])])/12
  }

  his_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predmix[,n]<-his_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + his_predf[,n]*(0.5*out[,1+(13+(n-1)*B)]+out[,1+(14+(n-1)*B)]+out[,1+(15+(n-1)*B)]+out[,1+(16+(n-1)*B)]+out[,1+(17+(n-1)*B)]+out[,1+(21+(n-1)*B)]+out[,1+(22+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }

  his_predv1<-his_predf1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predf1[,n]<-his_predf[,n]-his_predmix[,n]
    his_predv1[,n]<-his_predv[,n]-his_predmix[,n]
  }

  # Treated cases (His+other)(Uncomplicated cases accessing Health system (incidence +relapses+failures))
  trt_predf<-trt_predv<-trt_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    trt_predf[,n]<-rowSums(tran[,c(traind[c(15:17,28:30,33:35,39:41,46:48,53,56,59),n])])/12
    trt_predv[,n]<-rowSums(tran[,c(traind[c(80:82,104:109,111:116,121:126,131:136,140:145,150,153,156,159,162,165,211:216, 219,222, 225,228,231,234,237,240,243,246,249,252,255,258,261),n])])/12
  }

  trt_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    trt_predmix[,n]<-trt_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + trt_predf[,n]*(0.5*out[,1+(13+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }

  trt_predf1<- trt_predf-trt_predmix
  trt_predv1<- trt_predv-trt_predmix

  # Reported severe cases (successfully treated +deaths (trt+untrt))
  severe_predf<-severe_predv<-severe_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    severe_predf[,n]<-(tran[,traind[25,n]]+tran[,traind[50,n]])/12
    severe_predv[,n]<-(tran[,traind[94,n]]+tran[,traind[101,n]]+tran[,traind[147,n]])/12
  }

  severe_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    severe_predmix[,n]<-severe_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + severe_predf[,n]*(0.5*out[,1+(13+(n-1)*B)]+out[,1+(14+(n-1)*B)]+out[,1+(16+(n-1)*B)]+out[,1+(17+(n-1)*B)]+out[,1+(18+(n-1)*B)]+out[,1+(21+(n-1)*B)]+out[,1+(22+(n-1)*B)]+out[,1+(23+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }
  severe_predf1=severe_predf - severe_predmix
  severe_predv1=severe_predv - severe_predmix

  # fatalities
  fatal_pred<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    fatal_pred[,n]<-(tran[,traind[50,n]]+tran[,traind[147,n]])/12
  }

  # true clinical burden (all clinical cases+failures - uncomp+severe; trt+untreated)
  totclininc_predf<-totclininc_predv<-totclininc_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totclininc_predf[,n]<-rowSums(tran[,c(traind[c(14:17,27:35,38:41,45:48,52,53,55,56,58,59),n])])/12
    totclininc_predv[,n]<-rowSums(tran[,c(traind[c(76:82,103:117,120:126,130:136,139:145,149,150,152,153,155,156,158,159,161,162,164,165,210:216,218,219,221,222,224,225,227,228,230,231,233,234,236,237,239,240,242,243,245,246,248,249,251,252,254,255,257,258,260,261),n])])/12
  }

  totclininc_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totclininc_predmix[,n]<-totclininc_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + totclininc_predf[,n]*(0.5*out[,1+(13+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }
  totclininc_predf1=totclininc_predf- totclininc_predmix
  totclininc_predv1=totclininc_predv- totclininc_predmix

  # true incidence (all uncomplicated +severe+failures)
  totinc_predf<-totinc_predv<-totinc_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totinc_predf[,n]<-rowSums(tran[,c(traind[c(12:17,26:41,43:48,51:59),n])])/12
    totinc_predv[,n]<-rowSums(tran[,c(traind[c(74:82,102:126,128:145,148:165,208:261),n])])/12
  }

  totinc_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totinc_predmix[,n]<-totinc_predv[,n]*(1-(out[,1+(1+(n-1)*B)]+out[,1+(9+(n-1)*B)]+out[,1+(10+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]))
  }

  totinc_predf1=totinc_predf - totinc_predmix
  totinc_predv1=totinc_predv - totinc_predmix

  # ************************************************************************************* #
  # for predicting prevalence with different tests
  # ************************************************************************************* #
  # true prevalence
  totalinf_pred<-rowSums(out[,c(varind[2:5,],varind[12:15,])+1])

  prevalence_predf<-prevalence_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    prevalence_predf[,n]<-100*rowSums(out[,c(varind[2:5,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prevalence_predv[,n]<-100*rowSums(out[,c(varind[12:15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
  }
  prev_micro_predf<-prev_micro_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  prev_RDT_predf<- prev_RDT_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  prev_qPCR_predf<-prev_qPCR_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    # prevalence by Micro
    prev_micro_predf[,n]<-100*(sens_n_micro*out[,c(varind[2,n])+1]+sens_a_micro*out[,c(varind[3,n])+1]+sens_c_micro*out[,c(varind[4,n])+1]+sens_s_micro*out[,c(varind[5,n])+1]+sens_H_micro*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_micro_predv[,n]<-100*(vsens_n_micro*out[,c(varind[12,n])+1]+vsens_a_micro*out[,c(varind[13,n])+1]+vsens_c_micro*out[,c(varind[14,n])+1]+vsens_s_micro*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])

    # prevalence by RDT
    prev_RDT_predf[,n]<-100*(sens_n_RDT*out[,c(varind[2,n])+1]+sens_a_RDT*out[,c(varind[3,n])+1]+sens_c_RDT*out[,c(varind[4,n])+1]+sens_s_RDT*out[,c(varind[5,n])+1]+sens_H_RDT*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_RDT_predv[,n]<-100*(vsens_n_RDT*out[,c(varind[12,n])+1]+vsens_a_RDT*out[,c(varind[13,n])+1]+vsens_c_RDT*out[,c(varind[14,n])+1]+vsens_s_RDT*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])

    # prevalence by qPCR
    prev_qPCR_predf[,n]<-100*(sens_n_qPCR*out[,c(varind[2,n])+1]+sens_a_qPCR*out[,c(varind[3,n])+1]+sens_c_qPCR*out[,c(varind[4,n])+1]+sens_s_qPCR*out[,c(varind[5,n])+1]+sens_H_qPCR*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_qPCR_predv[,n]<-100*(vsens_n_qPCR*out[,c(varind[12,n])+1]+vsens_a_qPCR*out[,c(varind[13,n])+1]+vsens_c_qPCR*out[,c(varind[14,n])+1]+vsens_s_qPCR*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
  }

  return(cbind(vmw_predf,    #1
               vmw_predv,     #2
               vmw_predf1,  #3
               vmw_predv1,   #4
               vmw_predmix,    #5
               his_predf,    #6
               his_predv,    #7
               his_predf1,  #8
               his_predv1,   #9
               his_predmix,    #10
               fatal_pred,    #11
               severe_predf, #12
               severe_predv, #13
               severe_predf1, #14
               severe_predv1, #15
               severe_predmix, #16
               totclininc_predf,    #17
               totclininc_predv,    #18
               totclininc_predf1,    #19
               totclininc_predv1,    #20
               totclininc_predmix,    #21
               prevalence_predf,    #22
               prevalence_predv,    #23
               prev_micro_predf,    #24
               prev_micro_predv,    #25
               prev_RDT_predf,    #26
               prev_RDT_predv,    #27
               prev_qPCR_predf,    #28
               prev_qPCR_predv,    #29
               trt_predf, #30
               trt_predv, #31
               trt_predf1, #32
               trt_predv1, #33
               trt_predmix, #34
               totinc_predf, #35
               totinc_predv, #36
               totinc_predf1, #37
               totinc_predv1, #38
               totinc_predmix, #39
               totalinf_pred
  ))

}

