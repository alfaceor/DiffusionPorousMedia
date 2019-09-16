#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 

#strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ]
strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01"]
#strB2FList = strB2FList[::-1]

strHrad = "10"
strBrad = "11"
Directory = "RandomIniCond_02/"
flnPrefix = Directory + "rw3D_RegularPorousDepletion"
strNtrials = str(10**5)

################## FIGURES CONFIGURATION ##################
###### BEGIN PLOT DECORATION VARIABLES
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES

def lineStyle_axDiff(strB2B):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2B       = float(strB2B)
  lblB2B     = latex_float02( pB2B )
  etiqueta   = "$"+lblB2B+" $"
  expBB      = np.abs(np.log10( pB2B ))
  colval     = (expBB+4.0)/10.0
  colmap     = plt.cm.Blues
  colorcito  = colmap(colval)
  
  case = expBB
  if case == 0 :
    marcador = 'o'
  elif case == 1 :
    marcador = '*'
  elif case == 2 :
    marcador = 's'
  else:
    marcador = (int(case), 0, 0)
  
  estilo = {'linestyle'   : tipoLinea, # (0, (5, 10)), 
        'color'           : colorcito,
        'linewidth'       : anchoLinea,
        'label'           : etiqueta,
        'alpha'           : alfa,
        'marker'          : marcador,
        'markersize'      : 9,
        'markeredgewidth' : 2,
        'markeredgecolor' : colorcito
        } 
  return estilo #'solid'


#######################  Diffusion vs pBF  ############################
# TODO: Plot a figure and save the Difusion coefficient in organized files.
figScaDiff, axScaDiff = plt.subplots(1,1, figsize=(8,6))
#axR2 = ax
axScaDiff.set_xlim((1e-5,1e7))
axScaDiff.set_ylim((1e-1, 1e6))

axScaDiff.tick_params(axis='x',which='minor', top='off')
axScaDiff.tick_params(axis='x',which='major', top='on')
axScaDiff.tick_params(axis='x',which='minor', bottom='off')
axScaDiff.tick_params(axis='x',which='major', bottom='on')

axScaDiff.tick_params(axis='y',which='minor', right='off')
axScaDiff.tick_params(axis='y',which='major', right='on')
axScaDiff.tick_params(axis='y',which='minor', left='off')
axScaDiff.tick_params(axis='y',which='major', left='on')


axScaDiff.set_xlabel("$p_S t$", fontsize='30', rotation=0)
axScaDiff.set_ylabel("$ r_2 / (p_S \\times t)$", fontsize='30')

axScaDiff.set_xscale("log")
axScaDiff.set_yscale("log")

axScaDiff.set_xlabel("$p_{S}/p_{D}$", fontsize='30')
axScaDiff.set_ylabel("$D/p_{D}$", fontsize='30')

#### Inset Figure
left, bottom, width, height = [0.25, 0.5, 0.42, 0.42]
axDiff = figScaDiff.add_axes([left, bottom, width, height])
#axDiff.set_title("$ R_H = "+strHrad+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiff.set_xscale("log")
axDiff.set_yscale("log")
axDiff.set_xlabel("$p_{S}$") #, fontsize='30')
axDiff.set_ylabel("$D$", rotation=0) #, fontsize='30')
#axDiff.set_ylim((0,0.2))

axDiff.tick_params(axis='x', labelsize=10)
axDiff.tick_params(axis='y', labelsize=10)  


axDiff.tick_params(axis='both', which='major', left='on', right='on', top='on', bottom='on')  
axDiff.tick_params(axis='both', which='minor', left='off', right='off', top='off', bottom='off')  

scaB2B       = np.array([])
scaDiffCoeff = np.array([])
for b2f in range(len(strB2FList)):
  strB2F = strB2FList[b2f]
  pB2F   = float(strB2F)
  outPrefix = "D_vs_pBB__pBF_"+strB2F +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"

  # D vs pB2B
  #print outfln_DAT
  try:
    pB2B, DiffCoeff = np.loadtxt(outfln_DAT).T
    print pB2B, DiffCoeff
    estiloDiff = lineStyle_axDiff(strB2F)
    axDiff.plot(pB2B, DiffCoeff, **estiloDiff)
    scaXX = pB2B/pB2F
    scaYY = DiffCoeff/pB2F
    axScaDiff.plot(scaXX, scaYY, **estiloDiff)
    scaB2B       = np.append(scaB2B,       scaXX)
    scaDiffCoeff = np.append(scaDiffCoeff, scaYY)
    
  except ValueError as e:
    pass
      
tmpArr = np.vstack((scaB2B, scaDiffCoeff))
XX, YY = np.sort(tmpArr, axis=1)
#axScaDiff.plot( XX, YY, '-or')
i_min, i_max = getLimits(XX, 1e2, 1e7)
popt, pcov = fitPowLaw(XX[i_min:i_max], YY[i_min:i_max])
print "popt ",  popt[0], np.exp(popt[1])
#axScaDiff.plot(XX[i_min:i_max], powLaw(XX[i_min:i_max] , popt[0], popt[1]) , '--r')
#axScaDiff.plot(XX[i_min:i_max], XX[i_min:i_max], '--r')



#handles, labels = axDiff.get_legend_handles_labels()
#axDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.75,0.03), loc="lower left", title="$p_{D}$", prop={'size':10}, frameon=False, numpoints=1)

#figScaDiff.tight_layout()

handles, labels = axDiff.get_legend_handles_labels()
axScaDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.75,0.03), loc="lower left", title="$p_{D}$", prop={'size':10}, frameon=False, numpoints=1)

figScaDiff.tight_layout()

figScaDiff.savefig("Fig3_D"+strHrad+".pdf")
figScaDiff.savefig("DR"+strHrad+".eps")
plt.show()



