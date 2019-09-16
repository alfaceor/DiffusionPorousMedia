#!/usr/bin/python

# Objective: Analyze the case where 1> pBB > pBF, where there is confiment and high mobility inside a surface of the spheres.

import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
###### BEGIN PLOT DECORATION VARIABLES
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES

strF2B = "1.0" 
strB2F = "1.0" #"0.000001"
#strB2BList= ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "0.2", "0.5", "1.0"]
strB2BList= ["1.0"]
strHradius = "20"
strBradius = "21"
flnPrefix="rw3D_RegularPorousDepletion"
strNtrials = str(10**5)

fig, ax = plt.subplots(1,1, figsize=(8,6))
axR2 = ax
left, bottom, width, height = [0.6, 0.25, 0.3, 0.25]

axR2.set_xlabel("$t$", fontsize='30', rotation=0)
axR2.set_ylabel("$ \\langle \\Delta r^2 \\rangle $", fontsize='30')

axR2.set_xscale("log")
axR2.set_yscale("log")

#axR2.set_xlim(1e2,1e6)
#axR2.set_ylim(1e1,1e6)

axR2.set_title("$ p_{F \\to B} = "+strF2B+" $, $ p_{B \\to F} = "+strB2F+" $")


#pB2B = map(float,strB2B)
DiffCoeffLong=np.zeros(len(strB2BList))

colvals = np.abs(np.log( map(float, strB2BList)))+3.0
colvals = colvals/colvals.max()

for b2b in range(len(strB2BList)):
  #strB2F = strB2FList[nparam]
  strB2B = strB2BList[b2b]
  dataFile = flnPrefix+"_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHradius+"__Brad_"+strBradius+"__nt_"+strNtrials+".dat"
  print dataFile 
  data = np.loadtxt(dataFile)
  ttime = data[:,0]
  DelX2 = data[:,7]
  axR2. plot(ttime, DelX2, 'o', linewidth=2, color=plt.cm.Purples(colvals[b2b]), label="$ "+strB2B+" $")

axR2. plot(data[:,0], data[:,0], linewidth=2, color='k' )
#axIns.plot(data[:,0], 1.*data[:,0]+1 )

handles, labels = axR2.get_legend_handles_labels()
axR2.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.02,0.5), loc="lower left", title="$p_{B \\to B}$", prop={'size':12}, frameon=False)
fig.tight_layout()
outFilename = flnPrefix+"_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+"VAR"+"__Hrad_"+strHradius+"__Brad_"+strBradius+".pdf" 
print outFilename
plt.savefig(outFilename)


plt.show()



