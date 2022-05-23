## Simulations-MID : Instructions for replicating the simulations 
To replicate the simulations, as presented in Section 4 follow these steps : 
1. Clone the repo 
```sh
   git clone https://github.com/.git](https://github.com/apapan08/Simulations-MID.git
   ```
2. Open `MIDsimulations.Rproj` and run the `best_thresholds.R`. This R script will create the thresholds (see below) that will be used for the simulation (Section 3.2). 
Outputs : 
* ./Thresholds/constant/l2_best.csv
* ./Thresholds/constant/l2_inf.csv

3. Run Create_Signals.R. This R script will create the signals used for the simulation (Section 4)
Output: 
* ./Signals/signal.rds
* ./Signals/all_combos.csv

4. Run simulations.R 
Output : 
* ./results/results.rds

To replicate the simulations of Scenario (S2) as presented in Table 5, the required procedure is similar to the above. 
1. Run best_thresholds_linear.R
2. Run Create_signals_Linear.R
3. Run simulation_Linear.R

