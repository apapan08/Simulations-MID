library(data.table)
library(rlist)
library(foreach)

source("algorithms/additive_thr_L2_trend.R")
source("algorithms/additive_thr_Linf_trend.R")
source("HelpFunctions/SignalCreation.R")


# 1. Set the parameters 
parameters <- list (
  dimension = seq (1,50,1), # dimension grid
  thresholds = seq(0.8,1.4,0.5),# threshold grid
  m = 500, # number of simulations per setup  
  n = 700 , # length of each timeserie
  a = 0.95 # percentile 
) 


#This chunk of code will create a list of 50 elements (as many as the dimensions)
#Each element will have 500 objects (as many as the number of signals per dimension)
#Each Object is a timeseries of a specific lenght
#For exammple Allsignals[[3]][[5]] is a 3 dimension signal. It is the 5th iteration 
AllSignals <- foreach(i=1:length(parameters$dimension))%do%
  {
    TempDimension = parameters$dimension[i]
    AllIteration <- foreach(j=1:parameters$m)%do%
      {
        TempMatrix = create_matrix(n = parameters$n, d = TempDimension,sd = 1)
      }
    return(AllIteration)
    
  }


# Simulation 

allmethods <- foreach(method = c("additive_thr_L2_trend","additive_thr_Linf_trend")) %do%
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

l2best <- BestThresholds[method == "additive_thr_L2_trend",.(dimension,threshold)]
linfbest <- BestThresholds[method == "additive_thr_Linf_trend",.(dimension,threshold)]
fwrite(l2best,"Thresholds/linear/l2_best.csv")
fwrite(linfbest,"Thresholds/linear/linf_best.csv")

