#This program prints out data to *.dat files, for Fig. 10 in our paper
import numpy as np
from Zstats import Zdisc, Zexcl

#Fig 10: Discovery case [left panel]

#X-axis

Z_array = np.arange(0,5.005,0.005)

#Y-axis

#Set signal and background means
for (s, b) in [(1, 10), (10, 5)]:
    #Calculate P(Zdisc > Z)
    #Known background case
    temp1 = [Zdisc(s, b, asimov_only=False, Zcriteria=Z)[5] for Z in Z_array]

    #Unknown background case, with dbhat/bhat=0.5
    temp2 = [Zdisc(s, b, 0.5*b, asimov_only=False, Zcriteria=Z)[5] for Z in Z_array]

    #Printing data to a *.dat file
    np.savetxt('fig10_disc_s%s_b%s.dat' %(s,b),np.transpose([Z_array, temp1, temp2]),delimiter='\t',header='Z \t dbhat/bhat=0 \t dbhat/bhat=0.5',comments='#Fig10: The probability of obtaining a discovery significance Zdisc > Z in a large number of pseudo-experiments, when s=%s, and b=%s\n#' %(s,b))

#Fig 10: Exclusion case [right panel]

#X-axis

Z_array = np.arange(0,1.8+0.005,0.005)

#Y-axis

#Set signal and background means
for (s, b) in [(1, 25), (5, 5)]:
    #Calculate P(Zexcl > Z)
    #Known background case
    temp1 = [Zexcl(s, b, asimov_only=False, Zcriteria=Z)[5] for Z in Z_array]

    #Unknown background case, with dbhat/bhat=0.5
    temp2 = [Zexcl(s, b, 0.5*b, asimov_only=False, Zcriteria=Z)[5] for Z in Z_array]

    #Printing data to a *.dat file
    np.savetxt('fig10_excl_s%s_b%s.dat' %(s,b),np.transpose([Z_array, temp1, temp2]),delimiter='\t',header='Z \t dbhat/bhat=0 \t dbhat/bhat=0.5',comments='#Fig10: The probability of obtaining an exclusion significance Zexcl > Z in a large number of pseudo-experiments, when s=%s, and b=%s\n#' %(s,b))
