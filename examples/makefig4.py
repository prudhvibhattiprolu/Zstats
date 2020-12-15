#This program prints out data to *.dat files, for Fig. 4 in our paper
import numpy as np
from Zstats import Zdisc, Zexcl

#WARNING: Each computation, particularly when *asimov_only* is set to False, takes a lot more time when background uncertainty is non-zero

#And, the computation time increases as s/bhat (*sbyb*) gets smaller, and also when *s* gets larger.

#For faster computation, significantly reduce the number of points in *s_array* and/or run multiple batch jobs on a compute cluster.

#X-axis

s_array = np.append(np.arange(1.00,2.00,0.01),np.arange(2.00,10.02,0.02))
s_array = np.append(s_array,np.arange(10.05,100.05,0.05))

#Fig 4: Discovery case [left panel]

#Y-axis

for sbyb in [1, 10, 100]:
    #Calculate Asimov Z for fixed s/bhat

    #fractional background uncertainty is 0
    temp1 = np.transpose([Zdisc(s,s/sbyb) for s in s_array])

    #fractional background uncertainty is 0.2
    temp2 = np.transpose([Zdisc(s,s/sbyb,0.2*s/sbyb) for s in s_array])

    #fractional background uncertainty is 0.5
    temp3 = np.transpose([Zdisc(s,s/sbyb,0.5*s/sbyb) for s in s_array])

    #Printing data to a *.dat file
    np.savetxt('fig4_disc_sbyb%s.dat' %(sbyb),np.transpose([s_array, temp1, temp2, temp3]),delimiter='\t',header='s \t dbhat/bhat=0 \t dbhat/bhat=0.2 \t dbhat/bhat=0.5',comments='#Fig4: The exact Asimov expected significance for discovery, when s/bhat=%s\n#' %(sbyb))

#Fig 4: Exclusion case [right panel]

#Y-axis

for sbyb in [5, 0.5]:
    #Calculate Asimov Z for fixed s/bhat

    #fractional background uncertainty is 0
    temp1 = np.transpose([Zexcl(s,s/sbyb) for s in s_array])

    #fractional background uncertainty is 0.2
    temp2 = np.transpose([Zexcl(s,s/sbyb,0.2*s/sbyb) for s in s_array])

    #fractional background uncertainty is 0.5
    temp3 = np.transpose([Zexcl(s,s/sbyb,0.5*s/sbyb) for s in s_array])

    #Printing data to a *.dat file
    np.savetxt('fig4_excl_sbyb%s.dat' %(sbyb),np.transpose([s_array, temp1, temp2, temp3]),delimiter='\t',header='s \t dbhat/bhat=0 \t dbhat/bhat=0.2 \t dbhat/bhat=0.5',comments='#Fig4: The exact Asimov expected significance for exclusion, when s/bhat=%s\n#' %(sbyb))
