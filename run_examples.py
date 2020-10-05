import numpy as np
import matplotlib.pyplot as plt
import os

# choose model.  
#  0: nanodisc, 
#  1: micelle, 
#  2: coreshell with good prior, 
#  3: coreshell with poor prior.

i = 0

filename = ['Dataset5.rad','RBS_DDM30mM.dat','Isim.dat','Isim.dat']
headerlines = [0,3,2,2]
model = ['Nanodisc','MicelleDDM','CoreShellGoodPrior','CoreShellPoorPrior']
inputfile = ['nano','micelle','coreshell_good_prior','coreshell_poor_prior']

os.system('rm bayesfit')
os.system('gfortran -m64 -O2 bayesfit.f -o bayesfit')
os.system('./bayesfit input_%s.d' % inputfile[i])

q,I,dI = np.genfromtxt(filename[i],skip_header=headerlines[i],usecols=[0,1,2],unpack=True)
qprior,Iprior = np.genfromtxt('prior.d',skip_header=0,usecols=[0,1],unpack=True)
qfit,Ifit = np.genfromtxt('fit.d',skip_header=1,usecols=[0,1],unpack=True)

plt.errorbar(q,I,yerr=dI,marker='.',linestyle='none',color='red',zorder=0,label='data')
plt.plot(qprior,Iprior,color='grey',zorder=1,label='prior')
plt.plot(qfit,Ifit,color='black',zorder=2,label='fit')

plt.yscale('log')
plt.xscale('log')

plt.xlabel('q')
plt.ylabel('I(q)')
plt.legend()

plt.title(model[i])
plt.savefig(model[i])


plt.show()
