## CHS-speciation-variation-in-RI
### matlab code for simulations

###### 1. put all files in the working folder of matlab
   
###### 2. input_params.m is the file where you can adjust the parameter values accordingly
   
###### 3. after adjust the input parameteres in input_params.m, adjust the line 10 in the file main_hybridiz_model_parallel. Here you need to change the number in the line "p = parpool(23); " % here, 23 stands for the number of cores you can use for parallel running. the number of cores you put = number of your available cores-1
   
###### 4.run code 
```
main_hybridiz_model_parallel(seed,kappa)
```
######## seed is the seed number you assigned for this run, kappa is the temporal correlation parameter. kappa=0 means no temporal correlation, kappa=0.6 means temporal correlation.

###### 5. results.summary.m is the file for summarizing the results of matlab simulation runs, including for creating the heatmaps. 
   


