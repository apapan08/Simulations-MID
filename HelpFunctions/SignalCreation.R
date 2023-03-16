# Function for creating linear singal
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


# function that creates a random sample from 1 to k without replacement
rdu <- function(n,k) sample(1:k,n,replace=FALSE)
#EXAMPLE

#rdu(7,16)

##function to create signals 

create_matrix<-function(n,d,sd){
  model.sim = list(name = "wave", cpt.type = "pcwsLinContMean", cpt = c() ,
                    jump.size =c() , n = n, start = c(1, 1))
  timeserie = matrix(rnorm(d*n,sd = sd ),nrow=n)
  timeserie = timeserie + get.signal(model.sim) 
  return(timeserie)
}



#creates the underlying signal of length 'n' at positions 'chps'. jump between 
##changepoints is 'jump_vector'

create_signal=function(n,chps,jump_vector){
  #unique to prevent if a changepoint=n
  differences=diff(c(chps,n))
  signal_values=cumsum(jump_vector)
  output=rep(0,chps[1])
  for (i in 1:length(differences)){
    output=append(output,rep(signal_values[i],differences[i]))  
  }
  return(output)
}


#example 
#it will create signal with changepoints at 30 and 50
#create_signal(100,c(30,50),c(1,2))


##this function will create a vector of length n_chps in which its elements will have distance
##greater than `distance parameter`. Also will not output elements from 0-15 and from [n-15,-n]
place_chps<-function(len,n_chps,distance){
  min_dist <- distance
  a <- seq(1,len)
  picked <- integer(n_chps+2)
  copy <- a
  for (i in 1:(n_chps+2)) {
    if (length(copy)==0){
      result=sort(picked[picked!=0 & picked <(len-15) & picked > 15])
      min_value=min(n_chps,length(result))
      return(sort(sample(result,min_value)))
    }
    picked[i] <- sample(copy, 1)
    copy <- copy[abs(copy - picked[i]) >= min_dist]
    #print(copy)
  }
  result=sort(picked[picked!=0 & picked <(len-15) & picked > 15])
  min_value=min(n_chps,length(result))
  return(sort(sample(result,min_value)))
}


#EXAMPLE
#place_chps(1000,10,50)


#randomly changes sign of a vector
change_sign<-function(x){
  value<-sample(c(-1,1), size=length(x), replace=TRUE) * abs(x)
  return(value)
}


##EXAMPLE
#change_sign(c(1,2,4,3))

# 'random_matrix' creates a matrix of dimension 'd' and Length 'T'
# other parameters : number of changepoints , minimum distance between changepoints
# parameters of the uniform distribution that will be used for the random jump
random_matrix <- function(d,n,number_of_changepoints,sparsity,s = 2, 
                          noise_distr = "Gaussian",DoF = 8,
                          uniform_lower = -sqrt(3),uniform_upper = sqrt(3),
                          SettingSpatial = "setting2"){
  #sparsity(in how many columns there will it be at least one changepoint)
  ks=round(sparsity*d)
  components_for_changepoints<-matrix(NA,nrow = number_of_changepoints,ncol = d)
  #one change point will appear in those ks columns
  sparsity_max<-rdu(ks,d)
  #changepoint that will be placed in ks components
  max_chps=sample(c(1:number_of_changepoints),1)
  components_for_changepoints[max_chps,sparsity_max]=1
  for (i in c(1:number_of_changepoints)[-max_chps]){
    x=rdu(sample(c(1:ceiling(ks/16)),1),d)
    components_for_changepoints[i,x]=1
  }
  components_for_changepoints[is.na(components_for_changepoints)]<-0
  timeserie=list()
  #find the positions of the changepoints
  changepoints_positions <- seq(floor(n/(number_of_changepoints+1)), n - floor(n/(number_of_changepoints+1)), floor(n/(number_of_changepoints+1)))
  ind <- rowSums(components_for_changepoints)
  for (i in 1:d){
    x=which(components_for_changepoints[,i]==1)
    if (length(x)==0){
      timeserie=list.append(timeserie,rep(0,n))
    }else{
      changepoints<-changepoints_positions[x]
      jumps=change_sign(runif(length(changepoints),s,s))/sqrt(ind[x])
      timeserie=list.append(timeserie,create_signal(n,changepoints,jumps))
    }
    
    
  }
  noiseless_signal <- timeserie
  #adds noise to each one of the components
  if (noise_distr == "Gaussian"){
    timeserie = t(do.call(rbind,timeserie))+ matrix(rnorm(d*n),nrow = n)
  } 
  if (noise_distr == "Student"){
    timeserie = t(do.call(rbind,timeserie))+ sqrt((DoF-2)/DoF) * matrix(rt(d*n,df = DoF),nrow = n)
  }
  if (noise_distr == "Unif"){
    timeserie = t(do.call(rbind,timeserie))+ matrix(runif(d*n,min = uniform_lower,max = uniform_upper),nrow = n)
  } 
  
  if (noise_distr == "Spatial"){
    if (SettingSpatial == "setting1"){
      SigmaMatrix <- matrix(NA, d,d)
      for (i in 1: d){
        for (j in 1:d){
          SigmaMatrix[i,j] = 2^{-abs(i-j)}
        }
      }
    }else{
      SigmaMatrix <- matrix(1, d,d)
    }
    TempMatrix <- rmvnorm(n = n, mean = rep(0,d), sigma= SigmaMatrix)
    timeserie = t(do.call(rbind,timeserie)) + TempMatrix
  }
  return(list("ts"=timeserie,"chps"=changepoints_positions,"comp"=components_for_changepoints, "signal" = noiseless_signal))
}

# EXAMPLE
# d=10
# n=1000
# x=random_matrix(d,n,number_of_changepoints = 3,sparsity = 0.2,distance = 50)
# 
# A=x$ts
# chps=x$chps
# components=x$comp



# Linear signal 


##'random_matrix' creates a matrix of dimension 'd' and Length 'T'
random_matrix_linear <- function(d,n,number_of_changepoints,sparsity,distance,a_uniform=1,b_uniform=4){
  #sparsity
  ks=round(sparsity*d)
  components_for_changepoints<-matrix(NA,nrow = number_of_changepoints,ncol = d)
  #one change point will appear in those ks columns
  sparsity_max<-rdu(ks,d)
  #changepoint that will be placed in ks components
  max_chps=sample(c(1:number_of_changepoints),1)
  components_for_changepoints[max_chps,sparsity_max]=1
  for (i in c(1:number_of_changepoints)[-max_chps]){
    x=rdu(sample(c(1:ks),1),d)
    components_for_changepoints[i,x]=1
  }
  
  
  #define the "no-changepoints" model
  model.sim_no <- list(cpt.type = "pcwsLinContMean", cpt = c() ,
                       jump.size =c() , n = n, start = c(1, 1))
  
  
  timeserie=list()
  #find the positions of the changepoints
  changepoints_positions = place_chps(n,number_of_changepoints,distance)
  for (i in 1:d){
    x=which(components_for_changepoints[,i]==1)
    if (length(x)==0){
      timeserie=list.append(timeserie,get.signal(model.sim_no))
    }else{
      changepoints<-changepoints_positions[x]
      #jump here is the change in slope 
      jumps=change_sign(runif(length(changepoints),a_uniform,b_uniform))
      model.sim_chp <- list(cpt.type = "pcwsLinContMean", cpt = changepoints ,
                            jump.size = jumps , n = n, start = c(1, 1))
      timeserie=list.append(timeserie,get.signal(model.sim_chp))
    }
    
    
  }
  timeserie=t(do.call(rbind,timeserie))+matrix(rnorm(d*n),nrow = n)
  components_for_changepoints[is.na(components_for_changepoints)]<-0
  
  return(list("ts"=timeserie,"chps"=changepoints_positions,"comp"=components_for_changepoints))
}

# EXAMPLE
# d=10
# n=1000
# x=random_matrix_linear(d,n,number_of_changepoints = 3,sparsity = 0.2,distance = 50,a_uniform = 0.3,b_uniform = 0.6)
# 
# A=x$ts
# chps=x$chps
# components=x$comp

