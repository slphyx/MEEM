# ************************************************************************************* #
# Function to calculate transition rates, given variables and parameters
# ************************************************************************************* #
MM_Malrates <- function(x, input, parmal, t, ti) {
  ###should remove to somewhere!!
  N<-24   # number of patches
  B<-23   # number of variables per patch
  A<-261  # number of transitions per patch
  V<-N*B # total number of variables
  L<-N*A #total number of transitions
  startyear=1995 # starting year of simulation
  tyears<-20 # total years of simulation
  dtout<-1/12 # output timestep
  tsteps<-round(tyears/dtout) # number of time steps
  time<-startyear+seq(0,tyears,dtout) # time vector


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

  falpop<-1:10
  vivpop<-11:23
  ############################




  t_internal<-(ti-1)*dtout+t+startyear
  #Set up matrices
  seas<-c(rep(1,N)) # seasonality
  popf<-c(rep(0,N))    # population sizes
  popv<-c(rep(0,N))    # population sizes
  foif<-c(rep(0,N))     # forces of infection falciparum
  foiv<-c(rep(0,N))     # forces of infection vivax
  tranrate<-matrix(0,nrow=N,ncol=A)   # transition rate matrix
  itn<-c(rep(0,N))  #ITN coverage
  c_vmw<-c(rep(0,N)) #VMW coverage
  veff<-c(rep(0,N)) #VMW effect decline
  tau<-c(rep(0,N)) #Treatment probabilty
  tauo<-c(rep(0,N)) #Treatment with other
  tauh<-c(rep(0,N)) #Treatment with HIS
  tauv<-c(rep(0,N)) #Treatment with VMW

  eln<-approx(input$eln_t,input$eln_inp,xout=t_internal)$y
  for (n in 1:N){
    seas[n]<-1+eln*parmal$amp[n]*cos(2*pi*(t_internal-parmal$phi[n])) # seasonal forcing signal
    popf[n]<-sum(x[varind[falpop,n]]) # list of the variable indices for the human population
    popv[n]<-sum(x[varind[vivpop,n]]) # list of the variable indices for the human population
    itn[n]<-approx(input$itn_time,input$c_itn[,n],xout=t_internal)$y
    foif[n]<-sum(input$connect[n,]*seas[n]*((((1-itn[n])*parmal$bites)^2*parmal$prob_m*parmal$prob_h*(parmal$Pmaxf[n])/popf[n]*(parmal$zeta_n*x[varind[2,n]]+parmal$zeta_a*x[varind[3,n]]+x[varind[4,n]]+x[varind[5,n]])/sum(input$connect[n,]*popf))/((1-itn[n])*parmal$bites*parmal$prob_h*parmal$Pmaxf[n]/popf[n]+parmal$delta_m)*(parmal$gamma_m/(parmal$gamma_m+parmal$delta_m))))
    foiv[n]<-sum(input$connect[n,]*seas[n]*((((1-itn[n])*parmal$bites)^2*parmal$prob_m*parmal$vprob_h*(parmal$Pmaxv[n])/popv[n]*(parmal$vzeta_n*x[varind[12,n]]+parmal$vzeta_a*x[varind[13,n]]+x[varind[14,n]]+x[varind[15,n]])/sum(input$connect[n,]*popv))/((1-itn[n])*parmal$bites*parmal$vprob_h*parmal$Pmaxv[n]/popv[n]+parmal$delta_m)*(parmal$vgamma_m/(parmal$vgamma_m+parmal$delta_m))))

    foifa<-(1/(1/(foif[n])+1/parmal$gamma_h+1/parmal$gamma_m))
    foiva<-(1/(1/(foiv[n])+1/parmal$vgamma_h+1/parmal$vgamma_m))

    c_vmw[n]<-approx(input$vmw_time,input$cov_vmw[,n],xout=t_internal)$y
    veff[n]<-approx(input$vmw_time,input$vmw_eff[,n],xout=t_internal)$y
    tau[n]<-c_vmw[n]*input$sens_vmw*veff[n]+(1-c_vmw[n]*veff[n])*input$sens_his*parmal$eff_his+(1-c_vmw[n]*veff[n]-(1-c_vmw[n]*veff[n])*parmal$eff_his)*input$sens_oth*parmal$eff_oth #prob of treatment
    tauo[n]<-(1-c_vmw[n]*veff[n]-(1-c_vmw[n]*veff[n])*parmal$eff_his)*input$sens_oth*parmal$eff_oth #Prob treated by Other
    tauh[n]<-(1-c_vmw[n]*veff[n])*input$sens_his*parmal$eff_his #Prob treated by HIS
    tauv[n]<-c_vmw[n]*veff[n]*input$sens_vmw  #Prob treated by VMW

    t3=0
    if (t_internal>yrprim[n]){t3=1} #start radical cure from policy start date

    tranrate[n,]<-c(      #falciparum
      parmal$demog*popf[n]+(1-parmal$tausev)*parmal$pmort*parmal$nuq*x[varind[5,n]]+parmal$tausev*parmal$nuq*parmal$pmortt*x[varind[5,n]], # rate of birth 1
      parmal$demog*x[varind[1,n]],       # rate of death of S                                      2
      parmal$demog*x[varind[2,n]],       # rate of death of In                                     3
      parmal$demog*x[varind[3,n]],       # rate of death of Ia                                     4
      parmal$demog*x[varind[4,n]],       # rate of death of Ic                                     5
      parmal$demog*x[varind[5,n]],       # rate of death of Is                                     6
      parmal$demog*x[varind[6,n]],       # rate of death of To                                     7
      parmal$demog*x[varind[7,n]],       # rate of death of Tv                                     8
      parmal$demog*x[varind[8,n]],       # rate of death of Th                                     9
      parmal$demog*x[varind[9,n]],       # rate of death of R                                      10
      parmal$demog*x[varind[10,n]],       # rate of death of H                                     11
      parmal$psn*(1-parmal$ps)*foifa*x[varind[1,n]], #incidence S to In               12
      (1-parmal$psn)*(1-parmal$ps)*foifa*x[varind[1,n]],  #      incidence S to Ia    13
      (1-tau[n])*parmal$ps*foifa*x[varind[1,n]],          #      incidence S to Ic    14
      tauo[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to To    15
      tauv[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to Tv    16
      tauh[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to Th    17
      (1-parmal$tausev)*(1-parmal$pmort)*parmal$nuq*x[varind[5,n]],             #recovery Is to Ic    18
      (1-parmal$psev)*parmal$nuc*x[varind[4,n]],      # recovery Ic to Ia                         19
      parmal$nua*x[varind[3,n]],      # recovery Ia to In                                         20
      input$nun*x[varind[2,n]],      # recovery In to R                                           21
      (1-parmal$ptf)*parmal$nut*x[varind[6,n]],      # recovery To to H                                          22
      (1-parmal$ptf)*parmal$nut*x[varind[7,n]],      # recovery Tv to H                                          23
      (1-parmal$ptf)*parmal$nut*x[varind[8,n]],      # recovery Th to H                                          24
      parmal$tausev*parmal$nuq*(1-parmal$pmortt)*x[varind[5,n]],      # recovery Is to H                            25
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[2,n]], # incidence In to Ia          26
      parmal$pr*(1-tau[n])*foifa*x[varind[3,n]],         # incidence Ia to Ic          27
      parmal$pr*tauo[n]*foifa*x[varind[3,n]],         # incidence Ia to To             28
      parmal$pr*tauv[n]*foifa*x[varind[3,n]],         # incidence Ia to Tv             29
      parmal$pr*tauh[n]*foifa*x[varind[3,n]],         # incidence Ia to Th             30
      parmal$psev*parmal$nuc*x[varind[4,n]],         # incidence Ic to Is              31
      parmal$pr*(1-tau[n])*foifa*x[varind[2,n]],   # incidence In to Ic                32
      parmal$pr*tauo[n]*foifa*x[varind[2,n]],      # incidence In to To                33
      parmal$pr*tauv[n]*foifa*x[varind[2,n]],      # incidence In to Tv                34
      parmal$pr*tauh[n]*foifa*x[varind[2,n]],      # incidence In to Th                35
      parmal$prn*(1-parmal$pr)*foifa*x[varind[9,n]],     #       incidence R to In     36
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[9,n]],  #      incidence R to Ia     37
      (1-tau[n])*parmal$pr*foifa*x[varind[9,n]],          #      incidence R to Ic     38
      tauo[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to To     39
      tauv[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to Tv     40
      tauh[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to Th     41
      parmal$omega*x[varind[9,n]],                 #   loss of immunity from R to S    42
      parmal$prn*(1-parmal$pr)*foifa*x[varind[10,n]],     #       incidence H to In    43
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[10,n]],  #      incidence H to Ia    44
      (1-tau[n])*parmal$pr*foifa*x[varind[10,n]],          #      incidence H to Ic    45
      tauo[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to To    46
      tauv[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to Tv    47
      tauh[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to Th    48
      parmal$chi*x[varind[10,n]],                    #      loss HRP2 H to  R          49
      (1-parmal$tausev)*parmal$pmort*parmal$nuq*x[varind[5,n]]+parmal$tausev*parmal$nuq*parmal$pmortt*x[varind[5,n]],     # death untreated+treated Is  50
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Ia                 51
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Ic    52
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Th        53
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Ia                54
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Ic   55
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Tv       56
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[8,n]],          # failure Th to Ia             57
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[8,n]],  # failure Th to Ic        58
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[8,n]],    # failure Th to Th          59

      #vivax
      parmal$demog*popv[n]+(1-parmal$vtausev)*parmal$vpmort*parmal$vnuq*x[varind[15,n]]+parmal$vtausev*parmal$vpmortt*parmal$vnuq*x[varind[15,n]], # rate of birth 60
      parmal$demog*x[varind[11,n]],       # rate of death of S                                      61
      parmal$demog*x[varind[12,n]],       # rate of death of In                                     62
      parmal$demog*x[varind[13,n]],       # rate of death of Ia                                     63
      parmal$demog*x[varind[14,n]],       # rate of death of Ic                                     64
      parmal$demog*x[varind[15,n]],       # rate of death of Is                                     65
      parmal$demog*x[varind[16,n]],       # rate of death of To                                     66
      parmal$demog*x[varind[17,n]],       # rate of death of Tv                                     67
      parmal$demog*x[varind[18,n]],       # rate of death of Th                                     68
      parmal$demog*x[varind[19,n]],       # rate of death of R                                      69
      parmal$demog*x[varind[20,n]],       # rate of death of L                                      70
      parmal$demog*x[varind[21,n]],       # rate of death of Togd                                   71
      parmal$demog*x[varind[22,n]],       # rate of death of Tvgd                                   72
      parmal$demog*x[varind[23,n]],       # rate of death of Thgd                                   73
      parmal$vpsn*((1-parmal$vps)*(1-parmal$propgd[n])+(1-parmal$pgd)*parmal$propgd[n])*foiva*x[varind[11,n]], #incidence S to In          74
      (1-parmal$vpsn)*((1-parmal$vps)*(1-parmal$propgd[n])+(1-parmal$pgd)*parmal$propgd[n])*foiva*x[varind[11,n]],  # incidence S to Ia    75
      (1-tau[n])*(parmal$vps*(1-parmal$propgd[n])+parmal$pgd*parmal$propgd[n])*foiva*x[varind[11,n]],              #  incidence S to Ic    76
      (1-parmal$propgd[n])*parmal$vps*tauo[n]*foiva*x[varind[11,n]],                                     #         incidence S to To    77
      (1-parmal$propgd[n])*tauv[n]*parmal$vps*foiva*x[varind[11,n]],                                     #         incidence S to Tv    78
      (1-parmal$propgd[n])*tauh[n]*parmal$vps*foiva*x[varind[11,n]],                                     #         incidence S to Th    79
      parmal$propgd[n]*parmal$pgd*tauo[n]*foiva*x[varind[11,n]],                                       #         incidence S to Togd    80
      parmal$propgd[n]*parmal$pgd*tauv[n]*foiva*x[varind[11,n]],                                       #         incidence S to Tvgd    81
      parmal$propgd[n]*parmal$pgd*tauh[n]*foiva*x[varind[11,n]],                                       #         incidence S to Thgd    82
      (1-parmal$tausev)*(1-parmal$vpmort)*parmal$vnuq*x[varind[15,n]],            #          recovery Is to Ic    83
      (1-parmal$vpsev)*parmal$vnuc*x[varind[14,n]],                      # recovery Ic to Ia                      84
      parmal$vnua*x[varind[13,n]],                                      # recovery Ia to In                       85
      (1-parmal$vph)*input$vnun*x[varind[12,n]],                         # recovery In to R                       86
      parmal$vph*input$vnun*x[varind[12,n]],                             # recovery In to L                       87
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],        # recovery To to R         88
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],        # recovery Tv to R         89
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],        # recovery Th to R         90
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[21,n]],      # recovery Togd to R         91
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[22,n]],      # recovery Tvgd to R         92
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[23,n]],      # recovery Thgd to R         93
      (1-parmal$vph)*parmal$vtausev*(1-parmal$vpmortt)*parmal$vnuq*x[varind[15,n]],  # recovery Is to R    94
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],            # recovery To to L         95
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],           # recovery Tv to L          96
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],            # recovery Th to L         97
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[21,n]],          # recovery Togd to L         98
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[22,n]],        #  recovery Tvgd to L          99
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[23,n]],          # recovery Thgd to L         100
      parmal$vph*parmal$vtausev*(1-parmal$vpmortt)*parmal$vnuq*x[varind[15,n]], # recovery Is to L         101
      (1-parmal$vprn)*(1-parmal$vpr)*foiva*x[varind[12,n]],                   # incidence In to Ia         102
      parmal$vpr*(1-tau[n])*foiva*x[varind[12,n]],                     # incidence In to Ic                103
      (1-parmal$propgd[n])*parmal$vpr*tauo[n]*foiva*x[varind[12,n]],      # incidence In to To             104
      (1-parmal$propgd[n])*parmal$vpr*tauv[n]*foiva*x[varind[12,n]],      # incidence In to Tv             105
      (1-parmal$propgd[n])*parmal$vpr*tauh[n]*foiva*x[varind[12,n]],      # incidence In to Th             106
      parmal$propgd[n]*parmal$vpr*tauo[n]*foiva*x[varind[12,n]],        # incidence In to Togd             107
      parmal$propgd[n]*parmal$vpr*tauv[n]*foiva*x[varind[12,n]],        # incidence In to Tvgd             108
      parmal$propgd[n]*parmal$vpr*tauh[n]*foiva*x[varind[12,n]],        # incidence In to Thgd             109
      parmal$vpr*(1-tau[n])*foiva*x[varind[13,n]],                        # incidence Ia to Ic             110
      (1-parmal$propgd[n])*parmal$vpr*tauo[n]*foiva*x[varind[13,n]],         # incidence Ia to To          111
      (1-parmal$propgd[n])*parmal$vpr*tauv[n]*foiva*x[varind[13,n]],         # incidence Ia to Tv          112
      (1-parmal$propgd[n])*parmal$vpr*tauh[n]*foiva*x[varind[13,n]],         # incidence Ia to Th          113
      parmal$propgd[n]*parmal$vpr*tauo[n]*foiva*x[varind[13,n]],           # incidence Ia to Togd          114
      parmal$propgd[n]*parmal$vpr*tauv[n]*foiva*x[varind[13,n]],           # incidence Ia to Tvgd          115
      parmal$propgd[n]*parmal$vpr*tauh[n]*foiva*x[varind[13,n]],           # incidence Ia to Thgd          116
      parmal$vpsev*parmal$vnuc*x[varind[14,n]],                           # incidence Ic to Is             117
      parmal$vprn*(1-parmal$vpr)*foiva*x[varind[19,n]],                        #  incidence R to In        118
      (1-parmal$vprn)*(1-parmal$vpr)*foiva*x[varind[19,n]],                    # incidence R to Ia         119
      (1-tau[n])*parmal$vpr*foiva*x[varind[19,n]],                            #  incidence R to Ic         120
      (1-parmal$propgd[n])*tauo[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to To      121
      (1-parmal$propgd[n])*tauv[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to Tv      122
      (1-parmal$propgd[n])*tauh[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to Th      123
      parmal$propgd[n]*tauo[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Togd      124
      parmal$propgd[n]*tauv[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Tvgd      125
      parmal$propgd[n]*tauh[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Thgd      126
      parmal$vomega*x[varind[19,n]],                                   # loss of immunity from R to S      127
      parmal$vprn*(1-parmal$vpr)*(foiva)*x[varind[20,n]],                 #       incidence L to In        128
      (1-parmal$vprn)*(1-parmal$vpr)*(foiva)*x[varind[20,n]],              #      incidence L to Ia        129
      (1-tau[n])*parmal$vpr*(foiva)*x[varind[20,n]],                       #      incidence L to Ic        130
      (1-parmal$propgd[n])*tauo[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to To     131
      (1-parmal$propgd[n])*tauv[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to Tv     132
      (1-parmal$propgd[n])*tauh[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to Th     133
      parmal$propgd[n]*tauo[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Togd     134
      parmal$propgd[n]*tauv[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Tvgd     135
      parmal$propgd[n]*tauh[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Thgd     136
      parmal$vprn*(1-parmal$vpr)*parmal$vprel*parmal$vrel*x[varind[10,n]],                                #       relapse L to In   137
      (1- parmal$t2)*(1-parmal$vprn)*(1-parmal$vpr)*parmal$vprel*parmal$vrel*x[varind[10,n]],             #       relapse L to Ia   138
      (1- parmal$t2)*(1-tau[n])*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],                      #       relapse L to Ic   139
      (1- parmal$t2)*(1-parmal$propgd[n])*tauo[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],    #       relapse L to To   140
      (1- parmal$t2)*(1-parmal$propgd[n])*tauv[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],    #       relapse L to Tv   141
      (1- parmal$t2)*(1-parmal$propgd[n])*tauh[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],    #       relapse L to Th   142
      (1- parmal$t2)*parmal$propgd[n]*tauo[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],        #       relapse L to Togd 143
      (1- parmal$t2)*parmal$propgd[n]*tauv[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],        #       relapse L to Tvgd 144
      (1- parmal$t2)*parmal$propgd[n]*tauh[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[10,n]],        #       relapse L to Thgd 145
      parmal$vkappa*x[varind[10,n]],                      #        death of hypnozoites L to  S            146
      (1-parmal$vtausev)*parmal$vpmort*parmal$vnuq*x[varind[15,n]]+parmal$vtausev*parmal$vpmortt*parmal$vnuq*x[varind[15,n]],   # death untreated+treated Is  147

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[16,n]],                      # failed trt To to Ia         148
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],        # failed trt To to Ic         149
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],            # failed trt To to Th         150

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[17,n]],                      # failed trt Tv to Ia         151
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],        # failed trt Tv to Ic         152
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],            # failed trt Tv to Tv         153

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[18,n]],                      # failed trt Th to Ia         154
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],        # failed trt Th to Ic         155
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],            # failed trt Th to Th         156

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[21,n]],                    # failed trt Togd to Ia         157
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[21,n]],      # failed trt Togd to Ic         158
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[21,n]],          # failed trt Togd to Thgd       159

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[22,n]],                    # failed trt Tvgd to Ia         160
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[22,n]],      # failed trt Tvgd to Ic         161
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[22,n]],        # failed trt Tvgd to Tvgd         162

      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[23,n]],                    # failed trt Thgd to Ia         163
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[23,n]],      # failed trt Thgd to Ic         164
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[23,n]],        # failed trt Thgd to Thgd         165

      #Entanglements (dual treat)
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[2,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal In to H  166
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[3,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Ia to H  167
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[4,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Ic to H  168
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[5,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Is to H  169
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[12,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv In to R  170
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[13,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Ia to R  171
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[14,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Ic to R  172
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[15,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Is to R  173
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[12,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv In to L  174
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[13,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv Ia to L  175
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[14,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv Ic to L  176
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[15,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],           # dual trt Viv Is to L 177

      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[16,n]],                               # primaquine recovery To to R     178
      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[17,n]],                              #primaquine recovery Tv to R       179
      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[18,n]],                             #primaquine recovery Th to R        180
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[16,n]],                               # primaquine recovery To to L         181
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[17,n]],                              #primaquine recovery Tv to L           182
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[18,n]],                             #primaquine recovery Th to L            183

      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[21,n]],                  # test + act recovery Togd to R     184
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[22,n]],                 #test + act recovery Tvgd to R       185
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[23,n]],                #test + act recovery Thgd to R        186
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[21,n]],             # test + act recovery of GDef Togd to L      187
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[22,n]],            #test + act recovery of GDef Tvgd to L        188
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[23,n]],           #test + act recovery of GDef Thgd to L         189

      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[21,n]],        # test - primaquine recovery Togd to R     190
      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[22,n]],      # test - primaquine recovery Tvgd to R       191
      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[23,n]],      #test - primaquine recovery Thgd to R        192
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[21,n]],    # test - primaquine recovery of GDef Togd to L     193
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[22,n]],   #test - primaquine recovery of GDef Tvgd to L       194
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[23,n]],  #test - primaquine recovery of GDef Thgd to L        195

      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],                             # masked act recovery To to R     196
      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],                            #masked act recovery Tv to R       197
      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],                           #masked act recovery Th to R        198
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],                             # masked act recovery To to L         199
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],                            #masked act recovery Tv to L           200
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],                           #masked act recovery Th to L            201
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[21,n]],                           # masked act recovery Togd to R     202
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[22,n]],                          #masked act recovery Tvgd to R       203
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[23,n]],                         #masked act recovery Thgd to R        204
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[21,n]],                       # masked act recovery of GDef Togd to L     205
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[22,n]],                     #masked act recovery of GDef Tvgd to L        206
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[23,n]],                   #masked act recovery of GDef Thgd to L          207

      parmal$t2*parmal$vprn*(1-parmal$vpr)*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                       #       relapse L to In +triggering  208
      parmal$t2*(1-parmal$vprn)*(1-parmal$vpr)*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                   #      relapse L to Ia  +triggering  209
      parmal$t2*(1-tau[n])*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                           #      relapse L to Ic+triggering     210
      parmal$t2*(1-parmal$propgd[n])*tauo[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],     #          relapse L to To +triggering    211
      parmal$t2*(1-parmal$propgd[n])*tauv[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],   #           relapse L to Tv+triggering      212
      parmal$t2*(1-parmal$propgd[n])*tauh[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],     #          relapse L to Th  +triggering   213
      parmal$t2*parmal$propgd[n]*tauo[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],       #          relapse L to Togd +triggering    214
      parmal$t2*parmal$propgd[n]*tauv[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],      #          relapse L to Tvgd+triggering      215
      parmal$t2*parmal$propgd[n]*tauh[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],       #           relapse L to Thgd  +triggering  216

      #Failed treatments

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[16,n]],                            # primaquine failed trt To to Ia     217
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[16,n]],              # primaquine failed trt To to Ic     218
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[16,n]],                  # primaquine failed trt To to Th     219

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[17,n]],                           #primaquine failed trt Tv to Ia       220
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[17,n]],             #primaquine failed trt Tv to Ic       221
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[17,n]],                 #primaquine failed trt Tv to Tv       222

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[18,n]],                           #primaquine failed trt Th to Ia       223
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[18,n]],             #primaquine failed trt Th to Ic       224
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[18,n]],                 #primaquine failed trt Th to Th       225

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[16,n]],                                 # masked act failed trt To to Ia     226
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],                   # masked act failed trt To to Ic     227
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],                       # masked act failed trt To to Th     228

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[17,n]],                                 #masked act failed trt Tv to Ia      229
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],                   #masked act failed trt Tv to Ic      230
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],                       #masked act failed trt Tv to Tv      231

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[18,n]],                                 #masked act failed trt Th to Ia      232
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],                   #masked act failed trt Th to Ic      233
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],                       #masked act failed trt Th to Th      234

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[21,n]],                                 # masked act failed trt Togd to Ia   235
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[21,n]],                   # masked act failed trt Togd to Ic   236
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[21,n]],                     # masked act failed trt Togd to Thgd   237

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[22,n]],                                #masked act failed trt Tvgd to Ia     238
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[22,n]],                  #masked act failed trt Tvgd to Ic     239
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[22,n]],                    #masked act failed trt Tvgd to Tvgd     240

      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[23,n]],                               #masked act failed trt Thgd to Ia      241
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[23,n]],                 #masked act failed trt Thgd to Ic      242
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[23,n]],                   #masked act failed trt Thgd to Thgd      243

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],                  # test + act failed trt Togd to Ia     244
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],    # test + act failed trt Togd to Ic     245
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],        # test + act failed trt Togd to Thgd   246

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],                  #test + act failed trt Tvgd to Ia      247
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],    #test + act failed trt Tvgd to Ic      248
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],        #test + act failed trt Tvgd to Tvgd    249

      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],                  #test + act failed trt Thgd to Ia      250
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],    #test + act failed trt Thgd to Ic      251
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],        #test + act failed trt Thgd to Thgd    252

      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]],               # test - primaquine failed trt Togd to Ia       253
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]], # test - primaquine failed trt Togd to Ic       254
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]],     # test - primaquine failed trt Togd to Thgd     255

      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]],               # test - primaquine failed trt Tvgd to Ia       256
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]], # test - primaquine failed trt Tvgd to Ic       257
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]],     # test - primaquine failed trt Tvgd to Tvgd     258

      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]],               #test - primaquine failed trt Thgd to Ia        259
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]], #test - primaquine failed trt Thgd to Ic        260
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]]      #test - primaquine failed trt Thgd to Thgd      261

    )
  }
  return(c(t(tranrate)))
}

