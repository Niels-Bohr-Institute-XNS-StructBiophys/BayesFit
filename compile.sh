rm bayesfit
#gfortran -m64 -O3 bayesfit.f -o bayesfit
#gfortran bayesfit.f -march=native -O2 -o bayesfit
gfortran bayesfit.f -march=native -O3 -o bayesfit
