#This program prints out data to *.dat files, for Fig. 7 in our paper
import numpy as np
from Zstats import Zexcl

#WARNING: Each computation, particularly when *asimov_only* is set to False, takes a lot more time when background uncertainty is non-zero

#And, the computation time increases as s/bhat (*sbyb*) gets smaller, and also when *s* gets larger. Not recommended to run this script on a single cpu if generating data when background uncertainty is non-zero.

#For faster computation, significantly reduce the number of points in *s_array* and/or run multiple batch jobs on a compute cluster.

#X-axis

s_array = np.append(np.arange(1.00,2.00,0.01),np.arange(2.00,10.02,0.02))
s_array = np.append(s_array,np.arange(10.05,100.05,0.05))

#Fig 7: Known background case [left panel]

#Set fractional uncertainty in the background
dbbyb=0

#Y-axis

#Calculate P(Zexcl > 1.645) for exclusion case
#*Zcriteria* is set to 1.645 by default for function *Zexcl*

#s/bhat=0.1
temp1 = [Zexcl(s,s/0.1,dbbyb*s/0.1,asimov_only=False)[5] for s in s_array]

#s/bhat=1
temp2 = [Zexcl(s,s/1,dbbyb*s/1,asimov_only=False)[5] for s in s_array]

#s/bhat=10
temp3 = [Zexcl(s,s/10,dbbyb*s/10,asimov_only=False)[5] for s in s_array]

#Printing data to a *.dat file
np.savetxt('fig7_excl_dbbyb%s.dat' %(dbbyb),np.transpose([s_array, temp1, temp2, temp3]),delimiter='\t',header='s \t s/bhat=0.1 \t s/bhat=1 \t s/bhat=10',comments='#Fig7: The probability of obtaining a significance Zexcl > 1.645 in a large number of pseudo-experiments generated for the exclusion case, when dbhat/bhat=%s\n#' %(dbbyb))

#Fig 7: Uncertain background case, with dbhat/bhat=0.5 [right panel]

#Set fractional uncertainty in the background
dbbyb=0.5

#Y-axis

#Calculate P(Zexcl > 1.645) for exclusion case

#s/bhat=0.1
temp1 = [Zexcl(s,s/0.1,dbbyb*s/0.1,asimov_only=False)[5] for s in s_array]

#s/bhat=1
temp2 = [Zexcl(s,s/1,dbbyb*s/1,asimov_only=False)[5] for s in s_array]

#s/bhat=10
temp3 = [Zexcl(s,s/10,dbbyb*s/10,asimov_only=False)[5] for s in s_array]

#Printing data to a *.dat file
np.savetxt('fig7_excl_dbbyb%s.dat' %(dbbyb),np.transpose([s_array, temp1, temp2, temp3]),delimiter='\t',header='s \t s/bhat=0.1 \t s/bhat=1 \t s/bhat=10',comments='#Fig7: The probability of obtaining a significance Zexcl > 1.645 in a large number of pseudo-experiments generated for the exclusion case, when dbhat/bhat=%s\n#' %(dbbyb))
