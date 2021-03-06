
MM_EpiModel<-function(t,state, parode,input) # returns a list, need this for deSolve run
{
  with(as.list(c(state, parode)),
       {

         #   # ************************************************************************************* #
         #   # define variables
         #   # ************************************************************************************* #

         Z<-state
          
         # rates of change
         ti<-1
         V<-input$V
         transitions2<-input$transitions2
         
         transitionsiu1<-transitions2[,1]
         transitionsiu2<-transitions2[,2]
         transitionsiv1<-transitions2[,3]
         transitionsiv2<-transitions2[,4]
         L<-input$L
         N<-input$N
         
         transit<-MM_Malrates(Z[1:V],input,parode,Z[V+1],ti)

         if (sum(is.na(transit))>0)  {
           stop("transit NA   ",Z[V+1], "                                      ",
                as.data.frame(transit))
         }

         eq<-rep(0.0, V)

         eq<-EQ(L, N, eq, transit,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2)

         eq[V+1]<-1

         dZ<-eq

         # return the rate of change
         list(c(dZ))
       }
  )
  # end with(as.list ...
}



