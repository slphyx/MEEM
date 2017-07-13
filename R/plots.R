
###
# taken from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}


Replace0<-function(obdat, patchID, val,type="vmw"){
  #obdat is a data frame from MM_Plot (vmwdata,hisdata)
  #vmwdata<-data.frame(time=vmw_inc_t,vmw_incf=vmw_incf,
  #  vmw_incv=vmw_incv, vmw_incmix=vmw_incmix)
  # hisdata<-data.frame(his_inc_t=his_inc_t,his_incf=his_incf,
  #  his_incv=his_incv, his_incmix=his_incmix)
  #
  # val - value to replace 0
  # type - vmw or his

  tmpobdat<-NULL

  if(type=="vmw"){
    tmpobdat<-obdat
    eval(parse(text = paste0("tmpobdat$vmw_incf.X",patchID,"[tmpobdat$vmw_incf.X",patchID,"==0]<-val")))
    eval(parse(text = paste0("tmpobdat$vmw_incv.X",patchID,"[tmpobdat$vmw_incv.X",patchID,"==0]<-val")))
    eval(parse(text = paste0("tmpobdat$vmw_incmix.X",patchID,"[tmpobdat$vmw_incmix.X",patchID,"==0]<-val")))

    # tmpobdat[tmpobdat$vmw_incv[,patchID]==0,]<-val
    # tmpobdat[tmpobdat$vmw_incmix[,patchID]==0,]<-val
  }

  if(type=="his"){
    tmpobdat<-obdat
    eval(parse(text = paste0("tmpobdat$his_incf.X",patchID,"[tmpobdat$his_incf.X",patchID,"==0]<-val")))
    eval(parse(text = paste0("tmpobdat$his_incv.X",patchID,"[tmpobdat$his_incv.X",patchID,"==0]<-val")))
    eval(parse(text = paste0("tmpobdat$his_incmix.X",patchID,"[tmpobdat$his_incmix.X",patchID,"==0]<-val")))

    # tmpobdat[tmpobdat$his_incf[,patchID]==0,]<-val
    # tmpobdat[tmpobdat$his_incv[,patchID]==0,]<-val
    # tmpobdat[tmpobdat$his_incmix[,patchID]==0,]<-val
    #
  }

  if(is.null(tmpobdat)){
    stop("Please check your input.")
  }else{
    tmpobdat
  }
}


MM_Plot<-function(object, maldata, patchID=1, log10scale=FALSE){
  
  ## data
  if(file.exists(maldata)){
    alldata = loadWorkbook(maldata);
  }else{
    stop(paste("I can't find your data file ",maldata,"."));
  }
  
  if(is.null(object)){
    stop("Input object is NULL.")
  }
  
  
  # VMW cases over time
  vmw_casesf = readWorksheet(alldata, sheet="Sheet8")
  N<-ncol(vmw_casesf)-2
  
  if(patchID > N){
    stop("Max patch id : ",N)
  }
  cat("patch id: ",patchID,"/",N)

  vmw_incf <- vmw_casesf[,(1:N)+2]
  vmw_casesv = readWorksheet(alldata, sheet="Sheet10")
  vmw_incv <- vmw_casesv[,(1:N)+2]
  vmw_casesmix = readWorksheet(alldata, sheet="Sheet12")
  vmw_incmix <- vmw_casesmix[,(1:N)+2]
  vmw_inc_t <- vmw_casesf[,1]+vmw_casesf[,2]/12 # end of the month
  
  vmwdata <- data.frame(time=vmw_inc_t,vmw_incf=vmw_incf,vmw_incv=vmw_incv, vmw_incmix=vmw_incmix)
  #write.csv(vmwdata,file = "vmwdata.csv")
  
  # HIS cases over time
  his_casesf = readWorksheet(alldata, sheet="Sheet9")
  his_incf<-his_casesf[,(1:N)+2]
  his_casesv = readWorksheet(alldata, sheet="Sheet11")
  his_incv<-his_casesv[,(1:N)+2]
  his_casesmix = readWorksheet(alldata, sheet="Sheet13")
  his_incmix <- his_casesmix[,(1:N)+2]
  his_inc_t <- his_casesf[,1]+his_casesf[,2]/12 # end of the month
  
  hisdata <- data.frame(his_inc_t=his_inc_t,his_incf=his_incf,his_incv=his_incv, his_incmix=his_incmix)
  #write.csv(hisdata,file = "hisdata.csv")
  
  
  popdat <- readWorksheet(alldata,sheet="Sheet1")
  popdatscale <- 1/popdat[1:nrow(popdat),2]
  
  
  #define time vector
  startyear=1995 # starting year of simulation
  tyears<-20 # total years of simulation
  dtout<-1/12 # output timestep
  time<-startyear+seq(0,tyears,dtout)
  
  #object is the output from MM_RunMod
  vmw_predf1_fit<-object$vmw_predf1_fit
  vmw_predv1_fit<-object$vmw_predv1_fit
  vmw_predmix_fit<-object$vmw_predmix_fit
  his_predf1_fit<-object$his_predf1_fit
  his_predv1_fit<-object$his_predv1_fit 
  his_predmix_fit<-object$his_predmix_fit
  fatal_pred_fit<-object$fatal_pred_fit
  severe_predf1_fit<-object$severe_predf1_fit
  severe_predv1_fit<-object$severe_predv1_fit 
  severe_predmix_fit<-object$severe_predmix_fit
  if(!is.null("ppout")){
    ppout<-object$"ppout"
  }

  vmwmodel <- data.frame(time=time, vmw_predf1_fit=vmw_predf1_fit,
                       vmw_predv1_fit=vmw_predv1_fit, vmw_predmix_fit=vmw_predmix_fit)
  hismodel <- data.frame(time=time, his_predf1_fit=his_predf1_fit,
                       his_predv1_fit=his_predv1_fit, his_predmix_fit=his_predmix_fit)
 
  if(log10scale){
     vmwdata<-Replace0(obdat = vmwdata,patchID = patchID, val = popdatscale[patchID], type = "vmw")
     hisdata<-Replace0(obdat = hisdata,patchID = patchID, val = popdatscale[patchID], type = "his")
  }
  #write.csv(hisdata,file = "log10hisdata.csv")
  
  graph_vmw <- ggplot(data=vmwmodel,aes(x=time))+
    geom_line(data=vmwmodel,aes(y=vmw_predf1_fit[,patchID],color="model Pf"))+
    geom_point(data=vmwdata,aes(x=vmw_inc_t,y=vmw_incf[,patchID], color="data Pf"),shape=4)+
    geom_line(data=vmwdata,aes(x=vmw_inc_t, y=vmw_incf[,patchID], color = "data Pf"))+

    geom_point(data=vmwdata,aes(x=vmw_inc_t,y=vmw_incv[,patchID], color="data Pv"))+
    geom_line(data=vmwdata,aes(x=vmw_inc_t,y=vmw_incv[,patchID], color="data Pv"))+
    geom_line(data=vmwmodel,aes(x=time,y=vmw_predv1_fit[,patchID],color="model Pv"))+
     
    geom_point(data=vmwdata,aes(x=vmw_inc_t,y=vmw_incmix[,patchID], color="data mix"))+
    geom_line(data=vmwdata,aes(x=vmw_inc_t,y=vmw_incmix[,patchID], color="data mix"))+
    geom_line(data=vmwmodel, aes(x=time,y=vmw_predmix_fit[,patchID], color="model mix"))+
     
    ggtitle(label = paste0("VMW     Patch:",as.character(patchID)))+
    scale_colour_manual(name="",breaks = c("data Pf", "model Pf", "data Pv", "model Pv", "data mix", "model mix"),
                          values = c("data Pf" = "red", "model Pf"="red", "data Pv"="blue", "model Pv"="blue", 
                                     "data mix"="green", "model mix"="green"))+
       xlab("Year") +
       ylab("Incidence")+theme_bw()
    
  
  graph_his <- ggplot(data=hismodel,aes(x=time))+
    geom_line(data=hismodel,aes(y=his_predf1_fit[,patchID],color="model Pf"))+
    geom_point(data=hisdata,aes(x=his_inc_t,y=his_incf[,patchID], color="data Pf"),shape=3)+
    geom_line(data=hisdata,aes(x=his_inc_t, y=his_incf[,patchID], color = "data Pf"))+
    
    geom_point(data=hisdata,aes(x=his_inc_t,y=his_incv[,patchID], color="data Pv"),shape=4)+
    geom_line(data=hisdata,aes(x=his_inc_t,y=his_incv[,patchID], color="data Pv"))+
    geom_line(data=hismodel,aes(x=time,y=his_predv1_fit[,patchID],color="model Pv"))+
    
    geom_point(data=hisdata,aes(x=his_inc_t,y=his_incmix[,patchID], color="data mix"),shape=8)+
    geom_line(data=hisdata,aes(x=his_inc_t,y=his_incmix[,patchID], color="data mix"))+
    geom_line(data=hismodel, aes(x=time,y=his_predmix_fit[,patchID], color="model mix"))+
    
    ggtitle(label = paste0("HIS     Patch:",as.character(patchID)))+
    scale_color_manual(name="", breaks = c("data Pf", "model Pf", "data Pv", "model Pv", "data mix", "model mix"),
                       values = c("data Pf" = 'red', "model Pf"='red', "data Pv"='blue', "model Pv"='blue', 
                                  "data mix"='green', "model mix"='green'))+
    xlab("Year") +
    ylab("Incidence")+theme_bw()
    
  if(log10scale){
    grid_arrange_shared_legend(graph_vmw+scale_y_log10("Incidence", limits=c(1,max(vmwmodel))),
                               graph_his+scale_y_log10("Incidence", limits=c(1,max(hismodel))),ncol=2, nrow = 1)
  }else{  
    grid_arrange_shared_legend(graph_vmw,graph_his,ncol=2, nrow = 1)
  }
}


