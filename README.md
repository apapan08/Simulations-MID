## Simulations-MID : Instructions for replicating the simulations 
To replicate the simulations, as presented in Section 4 follow these steps : 
1. Clone the repo 
```sh
   git clone https://github.com/.git](https://github.com/apapan08/Simulations-MID.git
   ```
2. Open `MIDsimulations.Rproj` and run the `best_thresholds.R`. This R script will create the thresholds (see below) that will be used for the simulation (Section 3.3). 

* ./Thresholds/constant/l2_best.csv
* ./Thresholds/constant/l2_inf.csv

3. Execute `Create_Signals.R` to generate the signals required for the simulation (Section 4). 

* ./Signals/signal.rds
* ./Signals/all_combos.csv 

You can customize the script below to create signals corresponding to sections 4.2 (Spatially Dependent Data) and 4.3 (Non-Gaussian Noise) by adjusting the following parameters:

* `noise_distr`: Determines the noise distribution. Acceptable values include "Gaussian", "Student", "Uniform", and "Spatial".

* `DoF`: Applicable when `noise_distr` is set to "Student". Adjusts the degrees of freedom

* `SettingSpatial`: Applicable when `noise_distr` is set to "Spatial". Alters the settings as described in the paper.

* `uniform_lower` and `uniform_upper`: Applicable when `noise_distr` is set to "Uniform". Defines the lower and upper bounds of the uniform distribution.


```r
matrix_temp = random_matrix(
      d = dime,
      n = length_of_timeserie,
      number_of_changepoints = n_changepoints,
      sparsity = sparsity,
      s = 2.5,
      noise_distr = "Gaussian",
      DoF = 8,
      uniform_lower = -sqrt(3),
      uniform_upper = sqrt(3),
      SettingSpatial = "setting2"
    )
```



4. Run `simulations.R` to output the results in an rds format

* ./results/results.rds

Note that depending on the algorithms you want to run each time you can customize the following part of the code 

```r
parameters <-  list(
  iterations = 100,
  #same as the number of signals generated for each simulation setup
  algorithms_to_run = c(
    "MIDopt",
    "dc",
    "sbs",
    "inspect",
    "additive_thr_Linf_perm",
    "additive_thr_L2_perm",
    "KCP_Grid",
    "KCP_Scree",
    "ht_MIDopt"
  )
)
```



To replicate the simulations of Scenario (S2) as presented in Table 5, the required procedure is similar to the above. 
1. Run `best_thresholds_linear.R`
2. Run `Create_signals_Linear.R`
3. Run `simulation_Linear.R`

