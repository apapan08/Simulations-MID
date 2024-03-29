#Code for creating linear signals for simulations
library(rlist)
library(data.table)


source("HelpFunctions/SignalCreation.R")
#parameters
dimension=c(10,30,100)
number_of_changepoints<-c(3,20,50)
sparsity <- c(0.2,0.5,0.8)

#How many simulations to create 
m <- 100
#minimum distance between 2 changepoints
distance <- 20

simulation_setup <- as.data.table(tidyr::crossing(dimension,number_of_changepoints,sparsity))

#Depending on the number of changepoints we fix the length of the timeserie
simulation_setup[,LengthTimeserie := 1500]

fwrite(simulation_setup,"Signals/all_combos_linear.csv")

#list objest. Its elements is of the form dimension_number_of_chps_sparsity
timeserie_list=list()

for (i in 1:nrow(simulation_setup)){
  dime=simulation_setup[i,dimension]
  n_changepoints=simulation_setup[i,number_of_changepoints]
  sparsity=simulation_setup[i,sparsity]
  length_of_timeserie=simulation_setup[i,LengthTimeserie]
  timeserie_list[[paste0(
    as.character(dime),"_",as.character(n_changepoints),"_",as.character(sparsity)
  )]]=list()
  for (j in 1:m){
    matrix_temp=random_matrix_linear(d=dime,n=length_of_timeserie,number_of_changepoints = n_changepoints,
                              sparsity = sparsity,distance = distance,a_uniform = 0.3,b_uniform = 0.6)
    timeserie_list[[paste0(
      as.character(dime),"_",as.character(n_changepoints),"_",as.character(sparsity)
    )]][[as.character(j)]][["signal"]]=matrix_temp$ts
    timeserie_list[[paste0(
      as.character(dime),"_",as.character(n_changepoints),"_",as.character(sparsity)
    )]][[as.character(j)]][["chps"]]=matrix_temp$chps
    timeserie_list[[paste0(
      as.character(dime),"_",as.character(n_changepoints),"_",as.character(sparsity)
    )]][[as.character(j)]][["components"]]=matrix_temp$comp
  }
}

saveRDS(timeserie_list,file = "Signals/signals_linear.rds")


