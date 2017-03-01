
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


# ************************************************************************************* #
# import data; All data were removed from the Github version
# ************************************************************************************* #

alldata = loadWorkbook(system.file("extdata", "Cambodia_PfPv_provincial_data.xlsx", package="MEEM"))
climate = loadWorkbook(system.file("extdata", "climate.xlsx", package="MEEM"))

# population and villages
pvxy = readWorksheet(alldata, sheet="Sheet1")
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

# ************************************************************************************* #
# define transitions
# ************************************************************************************* #
# first transition is given without index
transitions =ssa.maketrans(V,rbind(varind[4,1], +1)) # birth mos patch 1
for (n in 1:N){
  # Falciparum indep flows (1:50)
  transitions[traind[1,n]]<-ssa.maketrans(V,rbind(varind[1,n],0,varind[1,n], +1)) # birth humans
  transitions[traind[2,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[1,n],0)) # death S=1
  transitions[traind[3,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1,varind[1,n],0)) # death In=2
  transitions[traind[4,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[1,n],0)) # death Ia=3
  transitions[traind[5,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1,varind[1,n],0)) # death Ic=4
  transitions[traind[6,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1,varind[1,n],0)) # death Is=5
  transitions[traind[7,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1,varind[1,n],0)) # death To=6
  transitions[traind[8,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1,varind[1,n],0)) # death Tv=7
  transitions[traind[9,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1,varind[1,n],0)) # death Th=8
  transitions[traind[10,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1,varind[1,n],0)) # death R=9
  transitions[traind[11,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1,varind[1,n],0)) # death H=10
  transitions[traind[12,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[2,n],+1)) # incidence S to In
  transitions[traind[13,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[3,n],+1)) # incidence S to Ia
  transitions[traind[14,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[4,n],+1)) # incidence S to Ic
  transitions[traind[15,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[6,n],+1)) # incidence S to To
  transitions[traind[16,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[7,n],+1)) # incidence S to Tv
  transitions[traind[17,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[8,n],+1)) # incidence S to Th
  transitions[traind[18,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[4,n], +1)) # recovery Is to Ic
  transitions[traind[19,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[3,n], +1)) # recovery Ic to Ia
  transitions[traind[20,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[2,n], +1)) # recovery Ia to In
  transitions[traind[21,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[9,n], +1)) # recovery In to R
  transitions[traind[22,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[10,n], +1)) # recovery To to H
  transitions[traind[23,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[10,n], +1)) # recovery Tv to H
  transitions[traind[24,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[10,n], +1)) # recovery Th to H
  transitions[traind[25,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[10,n], +1)) # recovery Is to H
  transitions[traind[26,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[3,n], +1)) # incidence In to Ia
  transitions[traind[27,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[4,n], +1)) # incidence Ia to Ic
  transitions[traind[28,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[6,n], +1)) # incidence Ia to To
  transitions[traind[29,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[7,n], +1)) # incidence Ia to Tv
  transitions[traind[30,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[8,n], +1)) # incidence Ia to Th
  transitions[traind[31,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[5,n], +1)) # incidence Ic to Is
  transitions[traind[32,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[4,n], +1)) # incidence In to Ic
  transitions[traind[33,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[6,n], +1)) # incidence In to To
  transitions[traind[34,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[7,n], +1)) # incidence In to Tv
  transitions[traind[35,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[8,n], +1)) # incidence In to Th
  transitions[traind[36,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[2,n], +1)) # incidence R to In
  transitions[traind[37,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[3,n], +1)) # incidence R to Ia
  transitions[traind[38,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[4,n], +1)) # incidence R to Ic
  transitions[traind[39,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[6,n], +1)) # incidence R to To
  transitions[traind[40,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[7,n], +1)) # incidence R to Tv
  transitions[traind[41,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[8,n], +1)) # incidence R to Th
  transitions[traind[42,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[1,n], +1)) # loss imm R to S
  transitions[traind[43,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[2,n], +1)) # incidence H to In
  transitions[traind[44,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[3,n], +1)) # incidence H to Ia
  transitions[traind[45,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[4,n], +1)) # incidence H to Ic
  transitions[traind[46,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[6,n], +1)) # incidence H to To
  transitions[traind[47,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[7,n], +1)) # incidence H to Tv
  transitions[traind[48,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[8,n], +1)) # incidence H to Th
  transitions[traind[49,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[9,n], +1)) # loss HRP2 H to  R
  transitions[traind[50,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[1,n],0)) # death Is fatal malaria
  transitions[traind[51,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[3,n],+1)) # failed trt To to Ia
  transitions[traind[52,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[4,n],+1)) # failed trt To to Ic
  transitions[traind[53,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[8,n],+1)) # failed trt To to Th
  transitions[traind[54,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[3,n],+1)) # failed trt Tv to Ia
  transitions[traind[55,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[4,n],+1)) # failed trt Tv to Ic
  transitions[traind[56,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[7,n],+1)) # failed trt Tv to Th
  transitions[traind[57,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[3,n],+1)) # failed trt Th to Ia
  transitions[traind[58,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[4,n],+1)) # failed trt Th to Ic
  transitions[traind[59,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[8,n],+1)) # failed trt Th to Th

  #Vivax indep flows (51-138)
  transitions[traind[60,n]]<-ssa.maketrans(V,rbind(varind[11,n],0,varind[11,n], +1)) # birth humans
  transitions[traind[61,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[11,n],0)) # death S=11
  transitions[traind[62,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1,varind[11,n],0)) # death In=12
  transitions[traind[63,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1,varind[11,n],0)) # death Ia=13
  transitions[traind[64,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1,varind[11,n],0)) # death Ic=14
  transitions[traind[65,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1,varind[11,n],0)) # death Is=15
  transitions[traind[66,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1,varind[11,n],0)) # death To=16
  transitions[traind[67,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1,varind[11,n],0)) # death Tv=17
  transitions[traind[68,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1,varind[11,n],0)) # death Th=18
  transitions[traind[69,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1,varind[11,n],0)) # death R=19
  transitions[traind[70,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1,varind[11,n],0)) # death L=20
  transitions[traind[71,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1,varind[11,n],0)) # death Togd=21
  transitions[traind[72,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1,varind[11,n],0)) # death Tvgd=22
  transitions[traind[73,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1,varind[11,n],0)) # death Thgd=23
  transitions[traind[74,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[12,n],+1)) # incidence S to In
  transitions[traind[75,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[13,n],+1)) # incidence S to Ia
  transitions[traind[76,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[14,n],+1)) # incidence S to Ic
  transitions[traind[77,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[16,n],+1)) # incidence S to To
  transitions[traind[78,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[17,n],+1)) # incidence S to Tv
  transitions[traind[79,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[18,n],+1)) # incidence S to Th
  transitions[traind[80,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[21,n],+1)) # incidence S to Togd
  transitions[traind[81,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[22,n],+1)) # incidence S to Tvgd
  transitions[traind[82,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[23,n],+1)) # incidence S to Thgd
  transitions[traind[83,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[14,n],+1)) # recovery Is to Ic
  transitions[traind[84,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[13,n],+1)) # recovery Ic to Ia
  transitions[traind[85,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[12,n],+1)) # recovery Ia to In
  transitions[traind[86,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[19,n],+1)) # recovery In to R
  transitions[traind[87,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[20,n],+1)) # recovery In to L
  transitions[traind[88,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n], +1)) # recovery To to R
  transitions[traind[89,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n], +1)) # recovery Tv to R
  transitions[traind[90,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n], +1)) # recovery Th to R
  transitions[traind[91,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n], +1)) # recovery Togd to R
  transitions[traind[92,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n], +1)) # recovery Tvgd to R
  transitions[traind[93,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n], +1)) # recovery Thgd to R
  transitions[traind[94,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[19,n], +1)) # recovery Is to R
  transitions[traind[95,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n], +1)) # recovery To to L
  transitions[traind[96,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n], +1)) # recovery Tv to L
  transitions[traind[97,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n], +1)) # recovery Th to L
  transitions[traind[98,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n], +1)) # recovery Togd to L
  transitions[traind[99,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n], +1)) # recovery Tvgd to L
  transitions[traind[100,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n], +1)) # recovery Thgd to L
  transitions[traind[101,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[20,n], +1)) # recovery Is to L
  transitions[traind[102,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[13,n], +1)) # incidence In to Ia
  transitions[traind[103,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[14,n], +1)) # incidence In to Ic
  transitions[traind[104,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[16,n], +1)) # incidence In to To
  transitions[traind[105,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[17,n], +1)) # incidence In to Tv
  transitions[traind[106,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[18,n], +1)) # incidence In to Th
  transitions[traind[107,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[21,n], +1)) # incidence In to Togd
  transitions[traind[108,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[22,n], +1)) # incidence In to Tvgd
  transitions[traind[109,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[23,n], +1)) # incidence In to Thgd
  transitions[traind[110,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[14,n], +1)) # incidence Ia to Ic
  transitions[traind[111,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[16,n], +1)) # incidence Ia to To
  transitions[traind[112,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[17,n], +1)) # incidence Ia to Tv
  transitions[traind[113,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[18,n], +1)) # incidence Ia to Th
  transitions[traind[114,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[21,n], +1)) # incidence Ia to Togd
  transitions[traind[115,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[22,n], +1)) # incidence Ia to Tvgd
  transitions[traind[116,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[23,n], +1)) # incidence Ia to Thgd
  transitions[traind[117,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[15,n], +1)) # incidence Ic to Is
  transitions[traind[118,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[12,n], +1)) # incidence R to In
  transitions[traind[119,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[13,n], +1)) # incidence R to Ia
  transitions[traind[120,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[14,n], +1)) # incidence R to Ic
  transitions[traind[121,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[16,n], +1)) # incidence R to To
  transitions[traind[122,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[17,n], +1)) # incidence R to Tv
  transitions[traind[123,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[18,n], +1)) # incidence R to Th
  transitions[traind[124,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[21,n], +1)) # incidence R to Togd
  transitions[traind[125,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[22,n], +1)) # incidence R to Tvgd
  transitions[traind[126,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[23,n], +1)) # incidence R to Thgd
  transitions[traind[127,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[11,n], +1)) # loss imm R to S
  transitions[traind[128,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # incidence L to In
  transitions[traind[129,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # incidence L to Ia
  transitions[traind[130,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # incidence L to Ic
  transitions[traind[131,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # incidence L to To
  transitions[traind[132,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # incidence L to Tv
  transitions[traind[133,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # incidence L to Th
  transitions[traind[134,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # incidence L to Togd
  transitions[traind[135,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # incidence L to Tvgd
  transitions[traind[136,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # incidence L to Thgd
  transitions[traind[137,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # relapse L to In
  transitions[traind[138,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # relapse L to Ia
  transitions[traind[139,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # relapse L to Ic
  transitions[traind[140,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # relapse L to To
  transitions[traind[141,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # relapse L to Tv
  transitions[traind[142,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # relapse L to Th
  transitions[traind[143,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # relapse L to Togd
  transitions[traind[144,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # relapse L to Tvgd
  transitions[traind[145,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # relapse L to Thgd
  transitions[traind[146,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[11,n], +1)) # death of hypnozoites L to S
  transitions[traind[147,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[11,n],0)) # death Is fatal malaria

  transitions[traind[148,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # failed trt To to Ia
  transitions[traind[149,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # failed trt To to Ic
  transitions[traind[150,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # failed trt To to Th

  transitions[traind[151,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # failed trt Tv to Ia
  transitions[traind[152,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # failed trt Tv to Ic
  transitions[traind[153,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # failed trt Tv to Tv

  transitions[traind[154,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # failed trt Th to Ia
  transitions[traind[155,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # failed trt Th to Ic
  transitions[traind[156,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # failed trt Th to Th

  transitions[traind[157,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # failed trt Togd to Ia
  transitions[traind[158,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # failed trt Togd to Ic
  transitions[traind[159,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # failed trt Togd to Thgd

  transitions[traind[160,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # failed trt Tvgd to Ia
  transitions[traind[161,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # failed trt Tvgd to Ic
  transitions[traind[162,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # failed trt Tvgd to Tvgd

  transitions[traind[163,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # failed trt Thgd to Ia
  transitions[traind[164,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # failed trt Thgd to Ic
  transitions[traind[165,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # failed trt Thgd to Thgd

  #Entanglements (166-)
  transitions[traind[166,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[10,n],+1)) # dual trt Fal In to H
  transitions[traind[167,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[10,n],+1)) # dual trt Fal Ia to H
  transitions[traind[168,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[8,n],+1)) # dual trt Fal Ic to H
  transitions[traind[169,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[6,n],+1)) # dual trt Fal Is to H
  transitions[traind[170,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[19,n],+1)) # dual trt Viv In to R
  transitions[traind[171,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[19,n],+1)) # dual trt Viv Ia to R
  transitions[traind[172,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[19,n],+1)) # dual trt Viv Ic to R
  transitions[traind[173,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[19,n],+1)) # dual trt Viv Is to R
  transitions[traind[174,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[20,n],+1)) # dual trt Viv In to L
  transitions[traind[175,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[20,n],+1)) # dual trt Viv Ia to L
  transitions[traind[176,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[20,n],+1)) # dual trt Viv Ic to L
  transitions[traind[177,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[20,n],+1)) # dual trt Viv Is to L
  transitions[traind[178,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n],+1)) # primaquine recovery To to R
  transitions[traind[179,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n],+1)) # primaquine recovery Tv to R
  transitions[traind[180,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n],+1)) # primaquine recovery Th to R

  transitions[traind[181,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n],+1)) # primaquine recovery To to L
  transitions[traind[182,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n],+1)) # primaquine recovery Tv to L
  transitions[traind[183,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n],+1)) # primaquine recovery Th to L
  transitions[traind[184,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # primaquine recovery Togd to R
  transitions[traind[185,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # primaquine recovery Tvgd to R
  transitions[traind[186,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # primaquine recovery Thgd to R
  transitions[traind[187,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # primaquine recovery Togd for Gdef to L
  transitions[traind[188,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Tvgd to L
  transitions[traind[189,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Thgd to L
  transitions[traind[190,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # primaquine recovery Togd to R
  transitions[traind[191,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # primaquine recovery Tvgd to R
  transitions[traind[192,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # primaquine recovery Thgd to R
  transitions[traind[193,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # primaquine recovery Togd for Gdef to L
  transitions[traind[194,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Tvgd to L
  transitions[traind[195,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Thgd to L
  transitions[traind[196,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n],+1)) # masked ACT recovery To to R
  transitions[traind[197,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n],+1)) # masked ACT recovery Tv to R
  transitions[traind[198,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n],+1)) # masked ACT recovery Th to R
  transitions[traind[199,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n],+1)) # masked ACT recovery To to L
  transitions[traind[200,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n],+1)) # masked ACT recovery Tv to L
  transitions[traind[201,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n],+1)) # masked ACT recovery Th to L
  transitions[traind[202,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # masked ACT recovery Togd to R
  transitions[traind[203,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # masked ACT recovery Tvgd to R
  transitions[traind[204,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # masked ACT recovery Thgd to R
  transitions[traind[205,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # masked ACT recovery Togd for Gdef to L
  transitions[traind[206,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # masked ACT recovery for Gdef Tvgd to L
  transitions[traind[207,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # masked ACT recovery for Gdef Thgd to L

  transitions[traind[208,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # relapse L to In + triggering
  transitions[traind[209,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # relapse L to Ia+ triggering
  transitions[traind[210,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # relapse L to Ic+ triggering
  transitions[traind[211,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # relapse L to To+ triggering
  transitions[traind[212,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # relapse L to Tv+ triggering
  transitions[traind[213,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # relapse L to Th+ triggering
  transitions[traind[214,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # relapse L to Togd+ triggering
  transitions[traind[215,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # relapse L to Tvgd+ triggering
  transitions[traind[216,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # relapse L to Thgd+ triggering

  #failed treatments
  transitions[traind[217,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # primaquine failed trt To to Ia
  transitions[traind[218,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # primaquine failed trt To to Ic
  transitions[traind[219,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # primaquine failed trt To to Th

  transitions[traind[220,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # primaquine failed trt Tv to Ia
  transitions[traind[221,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # primaquine failed trt Tv to Ic
  transitions[traind[222,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # primaquine failed trt Tv to Tv

  transitions[traind[223,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # primaquine failed trt Th to Ia
  transitions[traind[224,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # primaquine failed trt Th to Ic
  transitions[traind[225,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # primaquine failed trt Th to Th

  transitions[traind[226,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # masked act failed trt To to Ia
  transitions[traind[227,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # masked act failed trt To to Ic
  transitions[traind[228,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # masked act failed trt To to Th

  transitions[traind[229,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # masked act failed trt Tv to Ia
  transitions[traind[230,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # masked act failed trt Tv to Ic
  transitions[traind[231,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # masked act failed trt Tv to Tv

  transitions[traind[232,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # masked act failed trt Th to Ia
  transitions[traind[233,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # masked act failed trt Th to Ic
  transitions[traind[234,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # masked act failed trt Th to Th

  transitions[traind[235,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # masked failed trt Togd to Ia
  transitions[traind[236,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # masked failed trt Togd to Ic
  transitions[traind[237,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # masked failed trt Togd to Thgd

  transitions[traind[238,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # masked failed trt Tvgd to Ia
  transitions[traind[239,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # masked failed trt Tvgd to Ic
  transitions[traind[240,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # masked failed trt Tvgd to Tvgd

  transitions[traind[241,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # masked failed trt Thgd to Ia
  transitions[traind[242,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # masked failed trt Thgd to Ic
  transitions[traind[243,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # masked failed trt Thgd to Thgd

  transitions[traind[244,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # test + ACT failed trt Togd to Ia
  transitions[traind[245,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # test + ACT failed trt Togd to Ic
  transitions[traind[246,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # test + ACT failed trt Togd to Thgd

  transitions[traind[247,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # test + ACT failed trt Tvgd to Ia
  transitions[traind[248,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # test + ACT failed trt Tvgd to Ic
  transitions[traind[249,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # test + ACT failed trt Tvgd to Tvgd

  transitions[traind[250,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # test + ACT failed trt Thgd to Ia
  transitions[traind[251,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # test + ACT failed trt Thgd to Ic
  transitions[traind[252,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # test + ACT failed trt Thgd to Thgd

  transitions[traind[253,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # test - primaquine failed trt Togd to Ia
  transitions[traind[254,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # test - primaquine failed trt Togd to Ic
  transitions[traind[255,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # test - primaquine failed trt Togd to Thgd

  transitions[traind[256,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # test - primaquine failed trt Tvgd to Ia
  transitions[traind[257,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # test - primaquine failed trt Tvgd to Ic
  transitions[traind[258,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # test - primaquine failed trt Tvgd to Tvgd

  transitions[traind[259,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # test - primaquine failed trt Thgd to Ia
  transitions[traind[260,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # test - primaquine failed trt Thgd to Ic
  transitions[traind[261,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # test - primaquine failed trt Thgd to Thgd

}
#Alternate formulation of transitions matrix (used in epimodel function)
transitions2<-NULL
for (i in 1: length(transitions)){
  transitions2<-rbind(transitions2,cbind(as.integer(names(transitions[[i]]))[1],as.integer(names(transitions[[i]]))[2], transitions[[i]][1], transitions[[i]][2]))
}
row.names(transitions2)<-NULL
transitionsiu1<-transitions2[,1]
transitionsiu2<-transitions2[,2]
transitionsiv1<-transitions2[,3]
transitionsiv2<-transitions2[,4]



# ************************************************************************************* #
# Set the parameters (vivax-specific parameters prefixed with v)
# ************************************************************************************* #
pars = list(
  # Common parameters to both species
  Pmaxv=40000*seq_len(N)/seq_len(N), # the maximum populations of mosquitoes in each patch
  Pmaxf=40000*seq_len(N)/seq_len(N), # the maximum populations of mosquitoes in each patch
  delta_m=365.25/14, # death rate of mosquitoes
  amp=rep(1,N), # amplitude of seasonal forcing
  phi=rep(0.75,N), # phase angle of seasonal forcing
  bites=365.25/3, # biting rate
  prob_h=0.5, # probability that a bite will result in infection (mosquito to human)
  indexhet=rep(0.5,N), # parameter for connectivity about 1
  demog=1/50, # birth/death rate (mu)
  eff_itn=0.15, #efficacy of bednets
  stableveff=0.4, # steady state of vmw_effective
  hl_net=1.5, # half-life of bednets
  eff_vmw=0.6, # efficacy of VMW
  eff_his=0.25, # efficacy of HIS
  eff_oth=0.1, # efficacy of other systems
  elneff=21, # smoothing effect of el nino patterns.
  propgd=g6pDd, #Proportion of patch that is G6PDdef
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

);
