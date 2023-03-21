#Script to evaluate additive_linf_trend and additive_l2_trend

library(data.table)
library(foreach)
library(MixGHD)#for ARI


#best thresholds (Find through best_thresholds.R)
l2_dimension_lin<-fread("Thresholds/linear/l2_best.csv")
linf_dimension_lin<-fread("Thresholds/linear/linf_best.csv")

#Import Signals 
all_signal <- readRDS("Signals/signals_linear.rds")
#Import Evaluation Metrics + EstimatedSparsity
source("HelpFunctions/Evaluations.R")
#import methods + estimated sparsity functions
source("HelpFunctions/EstimateSparsityLinear.R")

#import methods 

# List all the R scripts in the folder
r_files <- list.files(path = "algorithms/Linear/",pattern = "*.R")

# Loop through the list of files and source each R script
for (file in r_files) {
  source(paste0("algorithms/Linear/",file))
}


parameters = list(
  iterations=100,#same as the number of signals generated for
  #each simulation setup ( m in Create_signals_Linear.R)
  algorithms_to_run = c("additive_thr_L2_trend","additive_thr_Linf_trend","MIDopt_lin")
)

Execution_Table <- data.table(
  methods = c(
    "additive_thr_L2_trend",
    "additive_thr_Linf_trend",
    "MIDopt_lin"
  ),
  Execution = c(
    "system.time(result<-additive_thr_L2_trend(signal_on_test))[[3]]",
    "system.time(result<-additive_thr_Linf_trend(signal_on_test))[[3]]",
    "system.time(result<-MIDopt_lin(signal_on_test))[[3]]"
  )
)

# Filter only algorithms of interest 

Execution_Table <- Execution_Table[methods %in% parameters$algorithms_to_run]

#import all scenarios
all_combos<-fread("Signals/all_combos_linear.csv")


Allresults<-foreach(i=1:nrow(all_combos))%do%
  {
    
    
    signal_name <- paste0(
      as.character(all_combos[i,dimension]),
      "_",as.character(all_combos[i,number_of_changepoints]),
      "_",as.character(all_combos[i,sparsity]))
    #all simulated signals of the current setup
    signal <- all_signal[[signal_name]]
    LengthTimeserie <- all_combos[i,LengthTimeserie]
    AllIteration <- foreach(iteration = 1:parameters$iterations)%do% 
      {
        signal_on_test=signal[[as.character(iteration)]]$signal
        chp_signal_on_test=signal[[as.character(iteration)]]$chps
        ResultsAllMethods <- foreach(method = 1:nrow(Execution_Table)) %do% 
          {
            exe <- Execution_Table[method,Execution]
            TimeAlgo <- eval(parse(text = exe))
            dt_results=data.table(
              setup = signal_name,
              iteration = iteration,
              method = Execution_Table[method,methods],
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
