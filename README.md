# BayesFit
Program for fitting analytical models to small-angle scattering data using Bayesian refinement.  

*obs: the code is not well documented, so if you wish to use BayesApp, but your model is not implemented, do not hesitate to contact Andreas Larsen for questions or a collaborative project, andreas(dot)larsen(a)bioch(dot)ox(dot)ac(dot)uk. 

## Citing the program  
If you BayesFit in your work, please cite:                       

Andreas Haahr Larsen, Lise Arleth, and Steen Hansen (2018). Analysis of small-angle scattering data using model fitting and Bayesian regularization. J Appl Cryst, 51, 1151-1161.  

## Install and run the program

### Operating systems
Only tested on Linux Ubuntu        

### Installation
Download bayesfit.f and compile it:

>> gfortran -m64 -O2 bayesfit.f -o bayesfit

#### Dependencies  
- fortran  
- gfortran (for compilation) 

### Running the program

>> bayesfit input_file.d    

### The input file must contain
Line 1: filename with data (headerlines are ignored)  
Line 2: model name (describes the function for fitting data, se below)  
Line 3: number of steps in function integrals, max iterations in minimization routine at each alpha (latter fixed to 100 in the last minimization step)
Line 4: q_min, q_max (to be included in the fit. to include everything, type, e.g., 0 and 1)  
Line 5 and following lines: parameter 1 prior value, parameter 1 prior width, 0 or 1 or 2 (0:fix 1:fit 2:fit with positive constraint)  

#### Example input file (for the for coreshell model): 

Isim.dat   
coreshell   
30 7   
0.0 1.0   
0.808 0.08 2   
4.95e-3 5e-4 1   
30.3 3 2    
49.5 5 2   
101 10 2   
3.96 0.4 2   

## Output

### Terminal Output
Column 1: chi^2/M, where M is the number of datapoints
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

##### parameters.d

prior and refined parameter values, as well as Chi-square, reduced Chi-square, final alpha, number of good parameters (Ng), number of fitted datapoints (M), evidence, value of constraint (S), and value of alpha*S.     

##### prior.d 
Intensity calculated with the prior values, which are also the initial values in the fitting algorithm.  

##### fit.d
Most probable fit after model refinement for all alpha-values.    

#### data.d     
Data used in the fitting process (no headerlines

## Models 
Four models have been implemented so far in BayesFit.   
The models can be run with the python3 script run_examples.py. Edit the line "i = 0" to select what example model to run, and run the script: 

>> python run_examples.py    

#### Nanodisc
Name of model (for input file): nanodisc  
Name of function (in BayesFit.f): fct_nano  
Description: Elliptical nanodiscs with a purification tag. The disc is build up of cylinder form factors and the tag is modelled as a random coil.  
Reference: Andreas Haahr Larsen, Lise Arleth, and Steen Hansen (2018). Analysis of small-angle scattering data using model fitting and Bayesian regularization. J Appl Cryst, 51, 1151-1161.     
Test data: Dataset5.rad    
Test input file: input_nano.d    
parameters:   
- Background [1/cm]       
- Concentration [uM]   
- Volume, lipid [A^3]   
- Volume, lipid tail [A^3]    
- Correction of volume, protein   
- Number of lipids per disc  
- Thickness of protein belt [A]   
- Roughness [A^2]     
- Area per lipid headgroup [A^2]    
- Ellipticity of lipid core     
- Number of waters per lipid  
- Radius of gyration of tag [A]    

#### Micelle
Name of model (for input file): micelle   
Name of function (in BayesFit.f): fct_mic   
Description: Core-shell detergent micelle, DDM micelles, concentration 30 mM.    
Reference: Andreas Haahr Larsen, Lise Arleth, and Steen Hansen (2018). Analysis of small-angle scattering data using model fitting and Bayesian regularization. J Appl Cryst, 51, 1151-1161.     
Test data: RB_DDM30mM.dat    
Test input file: input_micelle.d    
parameters:  
- Concentration [mM]   
- Background [1/cm]   
- Aggregation number, N   
- Contrast shell [e/A^3]*  
- Contrast core [e/A^3 ]*
- Roughness, sigma [A]    
- Ellipticity, epsilon   
e is the scattering lenght density of an electron    

#### CoreShell
Name of model (for input file): coreshell    
Name of cunction (in BayesFit.f): fct_coreshell    
Description: Idealised spherical core-shell particle.      
Reference:   
Test data: Isim.dat       
Test input file: input_coreshell_good_prior.d, input_coreshell_poor_prior.d    
parameters:
- I(0) [1/cm]    
- Constant backbround [1/cm]   
- Radius of core [A]    
- Radius of shell [A]   
- Excess scattering length density (contrast) core [arb. units]    
- Excess scattering length density (contrast) shell [arb. units]   

#### Test model of spheres
Name of model (for input file): test  
Name of function (in BayesFit.f): fct_test  
Description: Simple sphere with background.     

## Versions  

#### version 1
- First release, August 2018  

#### version 2
- Release September 2020    
- Cleaned up output    
- Added model CoreShell    

## License
BayesFit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.          

BayesFit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).  

## Acknowledgements
CoNeXT, the Carlsberg Foundation, and University of Copenhagen for co-funding the project.   

## Suggestions for future development  
- include resolution effects.  
- option for inclusion of several datasets.  
- restructure, so models are independent modules.  
