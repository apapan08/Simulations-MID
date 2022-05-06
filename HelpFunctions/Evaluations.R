
# 1. scaled Hausdorff distance
# n = length of signal
# chps = true changepoints
# est = estimated changepoints
scaled_ha<-function(n,chps,est){
  ns <- max(diff(c(0,chps,n)))
  diff <- abs(matrix(est,nrow=length(chps),ncol=length(est),byr=T)-matrix(chps,nrow=length(chps),ncol=length(est),byr=F))
  dh <- max(apply(diff,1,min),apply(diff,2,min))/ns
  return(dh)
}

#scaled_ha(1000,c(100,200,300,400),c(200,300))


# 2. Adjusted rank index
# n = length of timeserie
# changepoints = true changepoints 
label_timeserie<-function(n,changepoints){ # function to be used for ARI
  df = data.table::data.table(n=c(1:n))
  df[,class := findInterval(n,changepoints,left.open = TRUE)]
  return(df)
}

# Computes ARI
# n = length of timeserie
# true = true changepoints
# est = estimated changepoints

ARI_compute<- function(n,true,est){
  df_true=label_timeserie(n,true)
  df_est=label_timeserie(n,est)
  value=ARI(df_true$class,df_est$class)
  return(value)
}
# example
#ARI_compute(100,c(10,20),c(11,21))

