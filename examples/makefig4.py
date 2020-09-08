#This program prints out data to *.dat files, for Fig. 4 in our paper
import numpy as np
from Zstats import DeltaP

#X-axis

n_array = np.arange(0,31,1)

#Fig 4: Discovery case [left panel]

#Set signal and background means
s=5
b=5

#Y-axis

#Calculate probabilities *DeltaP(n,m,tau,s,b)* as a function of event count *n* in the signal region, for a fixed bhat = m/tau
#For more information about DeltaP, use the Python help function *help(DeltaP)*

#tau=Infinity
temp1 = [DeltaP(n,np.inf,np.inf,s,b) for n in n_array]

#tau=3
temp2 = [DeltaP(n,b*3,3,s,b) for n in n_array]

#tau=1
temp3 = [DeltaP(n,b*1,1,s,b) for n in n_array]

#Printing data to a *.dat file
np.savetxt('fig4_disc.dat',np.transpose([n_array, temp1, temp2, temp3]),delimiter='\t',header='n \t tau=Infinity \t tau=3 \t tau=1',comments='#Fig4: Probabilities DeltaP(n,m,tau,s) for discovery, when s=%s, and b=%s\n#' %(s,b))

#Fig 4: Exclusion case [right panel]

#Set signal and background means
s=0
b=10

#Y-axis

#Calculate probabilities *DeltaP(n,m,tau,s,b)* as a function of event count *n* in the signal region, for a fixed bhat = m/tau

#tau=Infinity
temp1 = [DeltaP(n,np.inf,np.inf,s,b) for n in n_array]

#tau=3
temp2 = [DeltaP(n,b*3,3,s,b) for n in n_array]

#tau=1
temp3 = [DeltaP(n,b*1,1,s,b) for n in n_array]

#Printing data to a *.dat file
np.savetxt('fig4_excl.dat',np.transpose([n_array, temp1, temp2, temp3]),delimiter='\t',header='n \t tau=Infinity \t tau=3 \t tau=1',comments='#Fig4: Probabilities DeltaP(n,m,tau,0) for exclusion, when s=%s, and b=%s\n#' %(s,b))
