library(data.table)
library(rlist)
library(foreach)

source("algorithms/forsimulations/additive_thr_L2.R")
source("algorithms/forsimulations/additive_thr_Linf.R")

# 1. Set the parameters 
parameters <- list (
  dimension = seq (1,50,1), # dimension grid
  thresholds = seq(0.5,2,0.05),# threshold grid
  m = 500, # number of simulations per setup  
  n = 700 , # length of each timeserie
  a = 0.95 # percentile 
) 


# 2. Create signals
allsignals <- list()

for (d in parameters$dimension){
  #create list for each dimension
  allsignals[[paste0(as.character(d))]]=list()
  for (iteration in 1:parameters$m)
    allsignals[[paste0(as.character(d))]][[paste0(as.character(iteration))]] =  matrix(rnorm(d*parameters$n),nrow = parameters$n)
}



# 3. Simulation 

allmethods <- foreach(method = c("additive_thr_L2","additive_thr_Linf")) %do%
  {
    alldimension <- foreach(d = parameters$dimension) %do% 
      {
        allthresholds <- foreach(t = parameters$thresholds) %do% 
          {
            alliteration <- foreach(i = 1:parameters$m) %do% 
              {
                temp_signal = allsignals[[paste0(as.character(d))]][[paste0(as.character(i))]]
                estimation = eval(parse(text = paste0(method,'(temp_signal,thr_const=t)')))
                dt_temp  =  data.table(
                  dimension = d,
                  iteration = i,
                  threshold = t,
                  method = method,
                  estimated = length(estimation)
                )
                return(dt_temp)
              }
            alliteration <- rbindlist(alliteration)
            return(alliteration)
          }
        allthresholds <- rbindlist(allthresholds)
        return(allthresholds)
      }
    alldimension <- rbindlist(alldimension)
    return(alldimension)
  }

allmethods <- rbindlist(allmethods)


allmethods[,IsZero := ifelse(estimated == 0,1,0)]

allmethodsgrouped <- allmethods[,.(Total = sum(IsZero)) , by = .(method,dimension,threshold)]



allmethodsgrouped[,diff := abs(parameters$a*parameters$m-Total)]
BestThresholds <- allmethodsgrouped[,.SD[diff == min(diff)],keyby = .(method,dimension)]

l2best <- BestThresholds[method == "additive_thr_L2",.(dimension,threshold)]
linfbest <- BestThresholds[method == "additive_thr_Linf",.(dimension,threshold)]
fwrite(l2best,"Thresholds/constant/l2_best.csv")
fwrite(linfbest,"Thresholds/constant/linf_best.csv")





