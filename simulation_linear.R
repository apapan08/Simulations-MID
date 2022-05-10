#Script to evaluate additive_linf_trend and additive_l2_trend

library(data.table)
library(foreach)
library(MixGHD)#for ARI


#best thresholds (Find through best_thresholds.R)
l2_dimension<-fread("Thresholds/linear/l2_best.csv")
linf_dimension<-fread("Thresholds/linear/linf_best.csv")

#Import Signals 
all_signal <- readRDS("Signals/signals_linear.rds")
#Import Evaluation Metric
source("HelpFunctions/Evaluations.R")


#import methods 

source("algorithms/additive_thr_L2_trend.R")
source("algorithms/additive_thr_Linf_trend.R")


parameters = list(
  iterations=100,#same as the number of signals generated for
  #each simulation setup ( m in Create_signals_Linear.R)
  algorithms_to_run = c("additive_thr_L2_trend","additive_thr_Linf_trend")
)


#import all scenarios
all_combos<-fread("Signals/all_combos_linear.csv")


Allresults <- foreach(i=1:nrow(all_combos))%do%
  {
    

    signal_name <- paste0(
      as.character(all_combos[i,dimension]),
      "_",as.character(all_combos[i,number_of_changepoints]),
      "_",as.character(all_combos[i,sparsity]))
    
    TempDimension <- all_combos[i,dimension]
    if (TempDimension >= 50 ){
      TempDimension <-  50
    }
    #Theshold to use 
    TempThresholdL2 <- l2_dimension[dimension == TempDimension,threshold]
    TempThresholdLinf <- linf_dimension[dimension == TempDimension,threshold]
    #all simulated signals of the current setup
    signal <- all_signal[[signal_name]]
    LengthTimeserie <- all_combos[i,LengthTimeserie]
    AllIteration <- foreach(iteration = 1:parameters$iterations)%do% 
      {
        signal_on_test=signal[[as.character(iteration)]]$signal
        chp_signal_on_test=signal[[as.character(iteration)]]$chps
        ResultsAllMethods <- foreach(method = 1:length(parameters$algorithms_to_run)) %do% 
          {
            TempMethod <- parameters$algorithms_to_run[method]
            if (TempMethod == "additive_thr_L2_trend"){
              TimeAlgo <- system.time(result <- additive_thr_L2_trend(signal_on_test,thr_const = TempThresholdL2))[[3]]
            }else{
              TimeAlgo <- system.time(result <- additive_thr_Linf_trend(signal_on_test,thr_const = TempThresholdLinf))[[3]]
            }
            dt_results=data.table(
              setup = signal_name,
              iteration = iteration,
              method = TempMethod,
              ARI = ARI_compute(LengthTimeserie,chp_signal_on_test,result),
              true_chps = length(chp_signal_on_test),
              est_chps = length(result),
              sh = scaled_ha(LengthTimeserie,chp_signal_on_test,result),
              time_algo=TimeAlgo)
            return(dt_results)
          }
        ResultsAllMethods <- rbindlist(ResultsAllMethods)
        return(ResultsAllMethods)
      }
    AllIteration = rbindlist(AllIteration)
    return(AllIteration)
    gc()
  }
Allresults <- rbindlist(Allresults)


dir.create(outpath <- paste0("results","/"),showWarnings = F)

saveRDS(Allresults,file=paste0(outpath,"/","results_linear.rds"))
