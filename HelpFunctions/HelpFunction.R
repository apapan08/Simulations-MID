get.signal <- function(model){
  
  if(model$cpt.type == "pcwsConstMean"){
    
    signal <- rep(0, model$n)
    
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    signal[segments[1,1]:segments[1,2]] <- model$start[1]
    
    for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
    
    
  }else if(model$cpt.type == "pcwsLinContMean"){
    if(length(model$cpt)==0){
      signal <- rep(0, model$n)
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      return(signal)}
    else{
      signal <- rep(0, model$n)
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      for(j in 2:nrow(segments)) {
        
        slope <- slope +  model$jump.size[j-1]
        
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      
      
    } 
    return(signal[1:model$n])
  }}



sim.model <- function(model, sigma=1){
  get.signal(model) + sigma * rnorm(model$n)} 


all.slopechanges.are.cpts <- function(x) {
  
  diff.x <- abs(diff(diff(x)))
  
  cpts <- which(diff.x > 0)
  no.of.cpt <- length(cpts)
  est <- x
  
  
  list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
  
  
}


##function to create signals 

create_matrix<-function(n,d,sd){
  model.sim = list(name = "wave", cpt.type = "pcwsLinContMean", cpt = c() ,
                    jump.size =c() , n = n, start = c(1, 1))
  timeserie = matrix(rnorm(d*n,sd = sd ),nrow=n)
  timeserie = timeserie + get.signal(model.sim) 
  return(timeserie)
}
