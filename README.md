# BayesFit
Fitting program for refining analytical models with small-angle scattering data using Bayesian refinement.  

### Cite
For description of the theory as perspectives of the program, see: 

Andreas Haahr Larsen, Lise Arleth, Steen Hansen, 2018, Analysis of small-angle scattering data using model fitting and Bayesian regularization, J. Appl. Cryst. 51, 1151-1161.

## Install and run the program

### Installation
Download bayesfit.f

### Compiling the program

>> gfortran -m32 -O2 bayefit.f -o bayesfit

### Running the program

>> bayesfit input_file.d    

### The input file  
Line 1: filename with data (headerlines are ignored)  
Line 2: model name (describes the function for fitting data)  
Line 3: log(alpha_min), log(alpha_max), number of alpha values, max number of steps in function integrals, max iterations*   
Line 4: q_min, q_max (to be included in the fit)  
Line 5 and following lines: parameter 1 prior value, parameter 1 prior width, 1 or 0 (fit or fix)  

*log(alpha_min) and log(alpha_max) and number of alpha values, are not used, as the optimal value of alpha is found automatically  

### Dependencies  
- fortran  
- gfortran (for compilation)  

## Output

### Terminal Output
while running the following is printet:  
Column 1: chi^2_r  
Column 2: lambda (used in the Levenberg-Marquardt algorithm)  
Column 3: log(alpha)  
Column 4: Q  
Column 5: dQ (change in Q from previous values - used in convergence criteria)  
Column 6: S (the prior value)  

For each value of alpha, the best fit is found and the following information is printet  
end lambda, iterations
refined parameter values and uncertainties  
prior parameter values and prior uncertainties  
final chi^2_r, -log(evidence), and Ng

At the end the following is printed:
optimal alpha and evidence at that value
computation time  

### Output files

##### parameters_filename
Column 1: log(alpha)  
Column 2: evidence  
Column 3: chi^2_r  
Column 4: N_g (number of good parameters)  
Column 5: alpha*S (S is the prior)  
Column 6: sum (the Occam term)  
Column 7: chi^2 (chi^2 = chi^2_r * mtot)  
Column 8: Iterations (before convergence)  
Column 9: mtot (number of datapoints)  
Column 10 and subsequent: refined parameter values  

The last line is the most probable fit  

##### start_filename 
Intensity calculated with the prior values, which are also the initial values in the fitting algorithm  

##### fit_filename 
Most probable fit after model refinement for all alpha-values. The last is the most probable fit.

##### global_par.d
Gives some parameters from the model, depending on the chosen model

##### errors.d
As parameters_<filename> but with errors on parameters instead of refined values  

#### Qcheck.d
Used for debugging

## Models 

Three models have been implemented so far

#### Nanodisc
Name of model (for input file): nanodisc  
Name of function (in BayesFit.f): fct_nano  
Description: Elliptical nanodiscs with a purification tag. The disc is build up of cylinder form factors and the tag is modelled as a random coil.  

#### Micelle
Name of model (for input file): micelle  
Name of function (in BayesFit.f): fct_mic  
Description: Core-shell detergent micelle.  

#### Test model of spheres
Name of model (for input file): test  
Name of function (in BayesFit.f): fct_test  
Description: Simple sphere with background.  

## Versions  

#### version 1
- Release August 2018  

## License
BayesFit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.          

BayesFit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).  

## Citing the program  
If you use the BayesFit software in your work, please cite:                       

Andreas Haahr Larsen, Lise Arleth, and Steen Hansen (2018). J Appl Cryst, 51, 1151-1161.  


## Acknowledgements
To CoNeXT and University of Copenhagen for co-funding the project.   

## Expected future development  
- include SANS contrasts.  
- include resolution effects.  
- option for inclusion of several datasets.  
- restructure, so models are independent modules.  
