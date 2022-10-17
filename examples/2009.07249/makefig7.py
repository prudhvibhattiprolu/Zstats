#This program prints out data to *.dat files, for Fig. 7 in our paper
import numpy as np
from Zstats import ZDiscExp, ZExclExp

#X-axis 

dbbyb_array=np.arange(0,0.5+0.005,0.005)

#Fig 7: Discovery case [left panel]

#Y-axis

#Set signal and background means
s=24
b=10

#Calculate Asimov Z, Mean Z, Median Z
temp = np.transpose([ZDiscExp(s,b,dbbyb*b,asimov_only=False, CLDiscbool=False) for dbbyb in dbbyb_array])

#Printing data to a *.dat file
np.savetxt('fig7_disc.dat',np.transpose([dbbyb_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig7: The exact Asimov, mean, and median expected significances for discovery, when s=%s, and bhat=%s\n#' %(s,b))

#Fig 7: Exclusion case [right panel]

#Y-axis

#Set signal and background means
s=12
b=20

#Calculate Asimov Z, Mean Z, Median Z
temp = np.transpose([ZExclExp(s,b,dbbyb*b,asimov_only=False, CLExclbool=False) for dbbyb in dbbyb_array])

#Printing data to a *.dat file
np.savetxt('fig7_excl.dat',np.transpose([dbbyb_array, temp[0], temp[1], temp[3]]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z',comments='#Fig7: The exact Asimov, mean, and median expected significances for exclusion, when s=%s, and bhat=%s\n#' %(s,b))
