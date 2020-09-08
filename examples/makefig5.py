#This program prints out data to *.dat files, for Fig. 5 in our paper
import numpy as np
from Zstats import Zdisc, Zexcl

#X-axis 

dbbyb_array=np.arange(0,0.5+0.005,0.005)

#Fig 2: Discovery case [left panel]

#Y-axis

#Set signal and background means
s=24
b=10

#Calculate Asimov Z, Mean Z, Asimov Z
temp = np.transpose([Zdisc(s,b,dbbyb*b,asimov_only=False) for dbbyb in dbbyb_array])

#Printing data to a *.dat file
np.savetxt('fig5_disc.dat',np.transpose([dbbyb_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig5: The exact Asimov, mean, and median expected significances for discovery, when s=%s, and bhat=%s\n#' %(s,b))

#Fig 2: Exclusion case [right panel]

#Y-axis

#Set signal and background means
s=12
b=20

#Calculate Asimov Z, Mean Z, Asimov Z
temp = np.transpose([Zexcl(s,b,dbbyb*b,asimov_only=False) for dbbyb in dbbyb_array])

#Printing data to a *.dat file
np.savetxt('fig5_excl.dat',np.transpose([dbbyb_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig5: The exact Asimov, mean, and median expected significances for exclusion, when s=%s, and bhat=%s\n#' %(s,b))
