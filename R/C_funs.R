# Call EQ function 
EQ<-function(L, N, oldeq, transit,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2){
  len<-length(oldeq)
  .C("EQ",
     as.integer(L),
     as.integer(N),
     as.double(oldeq),
     as.double(transit),
     as.integer(transitionsiu1),
     as.integer(transitionsiu2),
     as.integer(transitionsiv1),
     as.integer(transitionsiv2),
     as.double(vector("double",length = len))
  )[[9]]

}

