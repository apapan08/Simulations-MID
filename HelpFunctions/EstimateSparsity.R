# Functions for estimating the sparsity of a multivariate time serie signal (S1 scenario)

# 1. Function to estimate changepoints from triplets matrix

EstimateTriplets <- function(X,A,threshold){
  ChagePointsPerRow <- foreach(row = 1:nrow(A),.combine='c')%do%
    {
      triplet = A[row,]
    
      CusumValue = IDetect:::cusum_one(X,s = triplet[1],e = triplet[3],b = triplet[2])
      return(CusumValue)
    }
  #for linear check mad(diff(diff(X)))/sqrt(6)
  sigmaestimated = mad(diff(X))/sqrt(2)
  thresholdnew = threshold * sigmaestimated
  VectorPositionChangePoints <- which(ChagePointsPerRow>thresholdnew)
  FoundChangePoints <- A[VectorPositionChangePoints,2]
  return(FoundChangePoints)
  
}

EstimatedSparsity <- function(X,cpt){
  LengthTimeserie <- nrow(X)
  EstimatedChangepoints <- cpt
  if (length(EstimatedChangepoints)==0){
    ComponentsWithChp <- 0
  }else{
    triplets <- embed(c(1,EstimatedChangepoints,LengthTimeserie),3)[ ,3:1] 
    toAddMatrix <- matrix(rep(0,length(EstimatedChangepoints)*3),nrow = length(EstimatedChangepoints))
    toAddMatrix[2:nrow(toAddMatrix)] <- 1
    triplets <- triplets + toAddMatrix
    # we add one to changepoints
    threshold = 1.1*sqrt(2*log(LengthTimeserie)) # 1.1 is the default value for IDetect
    ChangepointComponents = apply(
      X, MARGIN = 2,
      FUN = EstimateTriplets,A = triplets,
      threshold = threshold )
    ComponentsWithChp = max(table(unlist(ChangepointComponents)))
  }
  sparsity <- ComponentsWithChp/ncol(X)
  return(sparsity)
}




