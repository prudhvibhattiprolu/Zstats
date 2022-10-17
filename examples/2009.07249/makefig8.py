#This program prints out data to *.dat files, for Fig. 8 in our paper
import numpy as np
from Zstats import ZDiscExp

#WARNING: Each computation, particularly when *asimov_only* is set to False, takes a lot more time when background uncertainty is non-zero

#And, the computation time increases as s/bhat (*sbyb*) gets smaller, and also when *s* gets larger. Not recommended to run this script on a single cpu if generating data when background uncertainty is non-zero.

#For faster computation, significantly reduce the number of points in *s_array* and/or run multiple batch jobs on a compute cluster.

#X-axis

s_array = np.append(np.arange(0.10,2.00,0.01),np.arange(2.00,10.02,0.02))
s_array = np.append(s_array,np.arange(10.05,100.05,0.05))

#Fig 8: Known background case [left panel]

#Set fractional uncertainty in the background
dbbyb=0

#Y-axis

#Calculate P(Zdisc > 5.0) for discovery case
#*Zcriteria* is set to 5.0 by default for function *ZDiscExp*

#s/bhat=2
temp1 = [ZDiscExp(s,s/2,dbbyb*s/2,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=5
temp2 = [ZDiscExp(s,s/5,dbbyb*s/5,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=10
temp3 = [ZDiscExp(s,s/10,dbbyb*s/10,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=50
temp4 = [ZDiscExp(s,s/50,dbbyb*s/50,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#Printing data to a *.dat file
np.savetxt('fig8_disc_dbbyb%s.dat' %(dbbyb),np.transpose([s_array, temp1, temp2, temp3, temp4]),delimiter='\t',header='s \t s/bhat=2 \t s/bhat=5 \t s/bhat=10 \t s/bhat=50',comments='#Fig8: The probability of obtaining a significance Zdisc > 5 in a large number of pseudo-experiments generated for the discovery case, when dbhat/bhat=%s\n#' %(dbbyb))

#Fig 8: Uncertain background case, with dbhat/bhat=0.5 [right panel]

#Set fractional uncertainty in the background
dbbyb=0.5

#Y-axis

#Calculate P(Zdisc > 5.0) for discovery case

#s/bhat=2
temp1 = [ZDiscExp(s,s/2,dbbyb*s/2,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=5
temp2 = [ZDiscExp(s,s/5,dbbyb*s/5,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=10
temp3 = [ZDiscExp(s,s/10,dbbyb*s/10,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#s/bhat=50
temp4 = [ZDiscExp(s,s/50,dbbyb*s/50,asimov_only=False,CLDiscbool=False)[5] for s in s_array]

#Printing data to a *.dat file
np.savetxt('fig8_disc_dbbyb%s.dat' %(dbbyb),np.transpose([s_array, temp1, temp2, temp3, temp4]),delimiter='\t',header='s \t s/bhat=2 \t s/bhat=5 \t s/bhat=10 \t s/bhat=50',comments='#Fig8: The probability of obtaining a significance Zdisc > 5 in a large number of pseudo-experiments generated for the discovery case, when dbhat/bhat=%s\n#' %(dbbyb))
