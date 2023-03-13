# KCPrs default -----------------------------------------------------------



kcprs_function_default <- function(data,wsize = 20,nperm =  100,Kmax = 100){
  if ((nrow(data))/2 < (Kmax)){Kmax <- nrow(data)/2 - 1}
  results <- kcpRS(
    data = data
    ,RS_name = "Mean"
    ,RS_fun = runMean
    ,wsize = wsize
    ,nperm = nperm
    ,Kmax = Kmax
  )
  
  allchangepoints <- results$changePoints
  
  return(allchangepoints)
}


kcprs_function_scree <- function(data,wsize = 20,nperm =  100,Kmax = 100){
  if ((nrow(data))/2 < (Kmax)){Kmax <- nrow(data)/2 - 1}
  results <- kcpRS(
    data = data
    ,RS_name = "Mean"
    ,RS_fun = runMean
    ,wsize = wsize
    ,nperm = nperm
    ,Kmax = Kmax
  )
  
  OptimalPoints <- results$changePoints_scree_test
  CPsGivenK <- results$CPs_given_K
  Changepoints <- CPsGivenK[CPsGivenK$k == OptimalPoints,startsWith(colnames(CPsGivenK),prefix = "CP")]
  Changepoints <- as.numeric(Changepoints[Changepoints != 0])
  return(Changepoints)
}
