#This program prints out data to *.dat files, for Fig. 5 in our paper
import numpy as np
from Zstats import ZDiscExp, ZExclExp

#Fig 5: Discovery case [left panel]

#X-axis

b_arraya = np.append(np.arange(0.1,2.00,0.01),np.arange(2.00,10.02,0.02))
b_arraya = np.append(b_arraya,np.arange(10.05,30.05,0.05))

#Y-axis

#set signal mean
s=6

#Calculate Z(16%, 50%, and 84% quantiles of the number of pseudo-experiments n)
temp1 = [ZDiscExp(s, b, asimov_only=False, quantile=0.16, CLDiscbool=False)[3] for b in b_arraya]
temp2 = [ZDiscExp(s, b, asimov_only=False, quantile=0.50, CLDiscbool=False)[3] for b in b_arraya]
temp3 = [ZDiscExp(s, b, asimov_only=False, quantile=0.84, CLDiscbool=False)[3] for b in b_arraya]

#Printing data to a *.dat file
np.savetxt('fig5_quantiles_disc_s%s.dat' %(s),np.transpose([b_arraya, temp1, temp2, temp3]),delimiter='\t',header='b \t 16% quantile \t 50% quantile (Median Z) \t 84% quantile',comments='#Fig5: The 16%%, 50%%, and 84%% quantile lines for the number of pseudo-experiments n for the discovery case when s=%s\n#' %(s))

#Fig 5: Exclusion case [right panel]

#X-axis

b_arrayb = np.append(np.arange(0.1,2.00,0.01),np.arange(2.00,10.02,0.02))
b_arrayb = np.append(b_arrayb,np.arange(10.05,100.05,0.05))

#Y-axis

#set signal mean
s=6

#Calculate Z(16%, 50%, and 84% quantiles of the number of pseudo-experiments n)
temp1 = [ZExclExp(s, b, asimov_only=False, quantile=0.16, CLExclbool=False)[3] for b in b_arrayb]
temp2 = [ZExclExp(s, b, asimov_only=False, quantile=0.50, CLExclbool=False)[3] for b in b_arrayb]
temp3 = [ZExclExp(s, b, asimov_only=False, quantile=0.84, CLExclbool=False)[3] for b in b_arrayb]

#Printing data to a *.dat file
np.savetxt('fig5_quantiles_excl_s%s.dat' %(s),np.transpose([b_arrayb, temp1, temp2, temp3]),delimiter='\t',header='b \t 16% quantile \t 50% quantile (Median Z) \t 84% quantile',comments='#Fig5: The 16%%, 50%%, and 84%% quantile lines for the number of pseudo-experiments n for the exclusion case when s=%s\n#' %(s))
