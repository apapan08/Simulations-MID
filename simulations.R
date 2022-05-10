# Script to compare the following algorithms/methods

#A)MID_OPT
#B)additive_thr_L2_perm
#C)additive_thr_Linf_perm
#D)InspectChangepoint
#E)DC
#F)SBS

library(data.table)
library(foreach)
library(MixGHD)#for ARI
library(hdbinseg)#for dcbs.alg
library(InspectChangepoint)#for inspect



#best thresholds (Find through best_thresholds.R)
l2_dimension<-fread("Thresholds/constant/l2_best.csv")
linf_dimension<-fread("Thresholds/constant/linf_best.csv")

#Import Signals 
# signals created through Create_signals.R
all_signal <- readRDS("Signals/signal.rds")

# Import Evaluation Metric functions
source("HelpFunctions/Evaluations.R")
 
#import methods + estimated sparsity functions
source("algorithms/additive_thr_Linf.R")
source("algorithms/additive_thr_L2.R")
source("HelpFunctions/EstimateSparsity.R")
source("algorithms/MIDopt.R")
source("algorithms/additive_permutations.R")


parameters <-  list(
  iterations=100,#same as the number of signals generated for each simulation setup
  algorithms_to_run = c("MIDopt","dc",
                        "sbs","inspect",
                        "additive_thr_Linf_perm","additive_thr_L2_perm")
)

Execution_Table <- data.table(
  methods = c("MIDopt",
              "additive_thr_L2_perm",
              "additive_thr_Linf_perm",
              "dc","sbs","inspect"),
  Execution = c("system.time(result<-MIDopt(signal_on_test))[[3]]",
                "system.time(result<-additive_thr_L2_perm(signal_on_test,m = 1000,a1 = 0.001,points = 20))[[3]]",
                "system.time(result<-additive_thr_Linf_perm(signal_on_test,m = 1000,a1 = 0.001,points = 20))[[3]]",
                "system.time(result<-dcbs.alg(t(signal_on_test),cp.type = 1, phi = -1)$ecp)[[3]]",
                "system.time(result<-sbs.alg(t(signal_on_test),cp.type = 1)$ecp)[[3]]",
                "system.time(result<-InspectChangepoint::inspect(t(signal_on_test))$changepoints[,1])[[3]]")
)

# Filter only algorithms of interest 

Execution_Table <- Execution_Table[methods %in% parameters$algorithms_to_run]


#import all scenarios (exported from Create_signals.R)
all_combos <- fread("Signals/all_combos.csv")

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

dir.create(outpath <- paste0("results","/"))

saveRDS(Allresults,file=paste0(outpath,"/","results.rds"))
