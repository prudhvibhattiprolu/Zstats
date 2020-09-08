#This program prints out data to *.dat files, for Fig. 1 in our paper
import numpy as np
from Zstats import Zdisc, Zexcl, ZdiscAsimovCCGV, ZexclAsimovKM

#Fig 1: Discovery case [left panel]

#X-axis

b_arraya = np.append(np.arange(0.1,2.00,0.01),np.arange(2.00,10.02,0.02))
b_arraya = np.append(b_arraya,np.arange(10.05,30.05,0.05))

#Y-axis

for s in [3, 6, 12]:
    #Computing Asimov Z, Mean Z, and Median Z
    temp1 = np.transpose([Zdisc(s,b,asimov_only=False) for b in b_arraya])

    #Computing CCGV Z
    temp2 = [ZdiscAsimovCCGV(s,b) for b in b_arraya]

    #Printing data to a *.dat file
    np.savetxt('fig1_disc_s%s.dat' %(s),np.transpose([b_arraya, temp1[0], temp1[1], temp1[3], temp2]),delimiter='\t',header='b \t Asimov Z \t Mean Z \t Median Z \t CCGV Z',comments='#Fig1: Expected significances for discovery, when s=%s, and background is known exactly\n#' %(s))

#Fig 1: Exclusion case [right panel]

#X-axis

b_arrayb = np.append(np.arange(0.1,2.00,0.01),np.arange(2.00,10.02,0.02))
b_arrayb = np.append(b_arrayb,np.arange(10.05,100.05,0.05))

#Y-axis

for s in [3, 6, 12]:
    #Computing Asimov Z, Mean Z, and Median Z
    temp1 = np.transpose([Zexcl(s,b,asimov_only=False) for b in b_arrayb])

    #Computing KM Z
    temp2 = [ZexclAsimovKM(s,b) for b in b_arrayb]

    #Printing data to a *.dat file
    np.savetxt('fig1_excl_s%s.dat' %(s),np.transpose([b_arrayb, temp1[0], temp1[1], temp1[3], temp2]),delimiter='\t',header='b \t Asimov Z \t Mean Z \t Median Z \t KM Z',comments='#Fig1: Expected significances for exclusion, when s=%s, and background is known exactly\n#' %(s))
