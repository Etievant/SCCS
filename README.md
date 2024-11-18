# SCCS

Code for replicating the simulations in "A Causal Framework for the Self-Controlled Case Series Design" by Etievant, Gail and Follmann (2024).

### Required packages 

```
tidyr, parallel, ggplot2.
```

### Scripts

* Script `Simul_Scenario1.R` replicates the simulation for scenario 1 proposed by Etievant, Gail and Follmann in Section 8.1 and Supplementary Section S4.2.

* Script `Simul_Scenario2.R` replicates the simulation for scenario 2 in Section 8.1. 

* Script `Simul_ATS.R` replicates the simulation in Section 8.2. 

* Script `Simul_Detection.R` replicates the simulation in Supplementary Section S4.1.


### Instructions to run each script

* Save the chosen script(s).

* Open and run the whole script(s).

* The results of the simulations are saved Rdata files and in figures. For example, when running script `Simul_Scenario1.R`, file contrasts_example1.pdf will give the same results as displayed in Figure 3 of Section 8.1. Folder SimulResults in the repository contains all the results files.
