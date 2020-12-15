#This program prints out data to *.dat files, for Fig. 3 in our paper
import numpy as np
from Zstats import Zdisc, Zexcl

#WARNING: Each computation, particularly when *asimov_only* is set to False, takes a lot more time when background uncertainty is non-zero

#And, the computation time increases as s/bhat (*sbyb*) gets smaller, and also when *s* gets larger. Not recommended to run this script on a single cpu as is.

#For faster computation, significantly reduce the number of points in *s_array* and/or run multiple batch jobs on a compute cluster. 

#X-axis 

s_array = np.append(np.arange(1.00,2.00,0.01),np.arange(2.00,10.02,0.02))
s_array = np.append(s_array,np.arange(10.05,100.05,0.05))

#Fig 3: Discovery case [left panel]

#Y-axis

#set fractional background uncertainty
dbbyb=0.2

for sbyb in [2, 10, 100]:
    #Calculate Asimov Z, Mean Z, Median Z for fixed s/bhat
    temp = np.transpose([Zdisc(s,s/sbyb,dbbyb*s/sbyb,asimov_only=False) for s in s_array])

    #Printing data to a *.dat file
    np.savetxt('fig3_disc_sbyb%s.dat' %(sbyb),np.transpose([s_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig3: The exact Asimov, mean, and median expected significances for discovery, when s/bhat=%s, and dbhat/bhat=%s\n#' %(sbyb,dbbyb))

#Fig 3: Exclusion case [right panel]

#Y-axis

for sbyb in [5, 0.5]:
    #Calculate Asimov Z, Mean Z, Median Z for fixed s/bhat
    temp = np.transpose([Zexcl(s,s/sbyb,dbbyb*s/sbyb,asimov_only=False) for s in s_array])

    #Printing data to a *.dat file
    np.savetxt('fig3_excl_sbyb%s.dat' %(sbyb),np.transpose([s_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig3: The exact Asimov, mean, and median expected significances for exclusion, when s/bhat=%s, and dbhat/bhat=%s\n#' %(sbyb,dbbyb))
