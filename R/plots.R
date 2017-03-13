
MM_Plot<-function(object, ToPlot="vmw", patchID=1){
  
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
  
  if(ToPlot=="vmw"){
    data<-data.frame(time=time,vmw_predf1_fit=vmw_predf1_fit)
      
      graph<-ggplot(data=data,aes(x=time))+
        geom_line(aes(y=vmw_predf1_fit[,patchID],color="predf1"))+
        geom_line(aes(y=vmw_predv1_fit[,patchID],color="predv1"))+
        ggtitle(label = paste0("Patch:",as.character(patchID)))+
        scale_color_manual("", breaks = c("predf1", "predv1"),
          values = c('predf1' = 'red', 'predv1' = "blue"))+
        xlab("year") +
        ylab("numbers XXXX")

  }
  graph
}