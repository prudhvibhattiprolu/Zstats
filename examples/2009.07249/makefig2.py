#This program prints out data to *.dat files, for Fig. 2 in our paper
import numpy as np
from Zstats import ZDiscExp, ZExclExp, ZDiscExpCCGV, ZExclExpKM

#Fig 2: Discovery case [left panel]

#X-axis 

s_arraya = np.arange(0,2.5+0.001,0.001)

#Y-axis

#set background mean
ba=10**(-6)

#Computing Asimov Z, Mean Z, and Median Z
temp1 = np.transpose([ZDiscExp(s, ba, asimov_only=0, CLDiscbool=False) for s in s_arraya])

#Computing CCGV Z
temp2 = [ZDiscExpCCGV(s, ba) for s in s_arraya]

#Printing data to a *.dat file
np.savetxt('fig2_disc_b0.dat',np.transpose([s_arraya, temp1[0], temp1[1], temp1[3], temp2]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z \t CCGV Z',comments='#Fig2: Expected significances for discovery for an extremely small background mean b=%s\n#' %(ba))

#Fig 2: Exclusion case [right panel]

#X-axis 

s_arrayb = np.arange(0,12+0.01,0.01)

#Y-axis

#set background mean
bb=0

#Computing Asimov Z, Mean Z, and Median Z
temp1 = np.transpose([ZExclExp(s, bb, asimov_only=0, CLExclbool=False) for s in s_arrayb])

#Computing KM Z
temp2 = [ZExclExpKM(s, bb) for s in s_arrayb]

#Printing data to a *.dat file
np.savetxt('fig2_excl_b0.dat',np.transpose([s_arrayb, temp1[0], temp1[1], temp1[3], temp2]),delimiter='\t',header='s \t Asimov Z \t Mean Z \t Median Z \t KM Z',comments='#Fig2: Expected significances for exclusion for zero background mean\n#')
