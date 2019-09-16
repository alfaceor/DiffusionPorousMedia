#!/usr/bin/python
import numpy as np
from SimulData import *

#flnPrefix="CenterIniCond/rw3D_RegularPorousDepletion"
#flnPrefix="RandomIniCond_02/teste/rw3D_RegularPorousDepletion"
flnPrefix="Over_100000_Steps/rw3D_RegularPorousDepletion"
strNtrials = str(10**5)


def ScaleFunction(sd):
  scaXX = sd.ttime/sd.Hrad**2.0
  scaYY = sd.DelR2/(sd.ttime)
  return scaXX, scaYY


strF2B  = "1.0" 
strB2F  = "0.00001" #"0.001"
strB2B  = "0.001" #"0.000001"
strHrad = "10" 
strBrad = "11" 

lblF2B = latex_float02( float(strF2B) )
lblB2F = latex_float02( float(strB2F) )
lblB2B = latex_float02( float(strB2B) )

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
###### BEGIN PLOT DECORATION VARIABLES
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES

import subprocess
import json

figR2, ax = plt.subplots(1,1, figsize=(8,6))

axR2    = ax
#axR2.set_title("$p_{FB} = "+lblF2B+", p_{BF} ="+lblB2F+", p_{BB} ="+lblB2B+", R_H = "+strHrad+" $")
axR2.text(10, 100, "$ p_{D} ="+lblB2F+", p_{S} ="+lblB2B+", R = "+strHrad+" $")
#axR2Sca.set_title("$p_{FB} = "+lblF2B+", p_{BB} \\gg p_{BF} ="+lblB2F+" $")
axR2.set_xscale("log")
axR2.set_yscale("log")

#axR2.set_xlim(1e0, 2e3)
#axR2.set_ylim(1e0, 2e3)


axR2.set_xlabel("$t$", fontsize='30', rotation=0)
axR2.set_ylabel("$ r_2 $", fontsize='30')

#axR2Sca.tick_params(axis='y',which='minor',left='off')
axR2.tick_params(axis='x',which='minor',bottom='off')
axR2.tick_params(axis='y',which='minor',right='off')
axR2.tick_params(axis='y',which='minor',left='off')
axR2.tick_params(axis='y',which='major',right='on')

#left, bottom, width, height = [0.65, 0.56, 0.3, 0.3]
#left, bottom, width, height = [0.58, 0.527, 0.35, 0.35]
#axR2 = figR2.add_axes([left, bottom, width, height])
#axR2.patch.set_alpha(0.0)
#axR2.set_xlabel("$t$", fontsize='18', rotation=0)
#axR2.set_ylabel("$ \\langle \\Delta r^2 \\rangle $", fontsize='18')
#
#axR2.set_xscale("log")
#axR2.set_yscale("log")
#axR2.xaxis.get_offset_text().set_fontsize(10)
#axR2.yaxis.get_offset_text().set_fontsize(10)
#
##axR2.set_xlim(1e0,1e8)
##axR2.set_ylim(1e0,9e6)
#
#axR2.tick_params(axis='x', labelsize=10)
#axR2.tick_params(axis='y', labelsize=10)  
#axR2.tick_params(axis='x',which='minor',bottom='off')
#axR2.tick_params(axis='y',which='minor',left='off')


#############################################
def lineStyDAT(sd):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  etiqueta   = "$"+sd.strHrad+" $"
  pB2B       = float(sd.strB2B)
  lblB2B     = latex_float( pB2B )
  expBB      = np.abs(np.log10( pB2B ))

  #colval = (expBF + 10.*expBB) / 36.
  if sd.strHrad == "10":
    colmap = plt.cm.Greens
  elif sd.strHrad == "20":
    colmap = plt.cm.Oranges
  elif sd.strHrad == "40":
    colmap = plt.cm.Blues
  elif sd.strHrad == "80":
    colmap = plt.cm.RdPu
  elif sd.strHrad == "160":
    colmap = plt.cm.Wistia # Yellows?
  else:
    colmap = plt.cm.Greens
  #colval  = 1.0 - expBF/6.0
  colval = (expBB+4.0)/10.0 #+3.0
  #print "colval = ",colval
  colorcito  = colmap(colval) #0.2+0.8*expBB)

#  if expBB == 0 : 
#    tipoLinea = '-'
#  elif expBB == 1 : 
#    tipoLinea = '--'
#  elif expBB == 2 : 
#    tipoLinea = '-.' 
#  elif expBB == 3 : 
#    tipoLinea = ':' 
#  elif expBB == 4 : 
#    tipoLinea = '-.' 
#  else:
#    tipoLinea = ':' 
  
  
  
  estilo = {'linestyle'   : tipoLinea, # (0, (5, 10)), 
        'color'           : colorcito,
        'linewidth'       : anchoLinea,
        'label'           : etiqueta,
        'alpha'           : alfa
#        'marker'          : marcador,
#        'markersize'      : 9,
#        'markeredgewidth' : 2,
#        'markeredgecolor' : colorcito,
        } 
  return estilo #'solid'


def plotSimul(axR2,flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials ):
  lblB2B = latex_float( float(strB2B) )
  sd = SimulData()
  sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=True )
  #sd.getData()
  estilo = lineStyDAT(sd)
  linea, = axR2.plot(sd.ttime, sd.DelR2, **estilo)
  popt, pcov, i_min, i_max = fitPowLawFromTo(sd.ttime, sd.DelR2, xmin=1e4, xmax=1e6, returnIndx=True)
  print popt[0], popt[1]
  #axR2.plot(sd.ttime, sd.DelR2[0]*sd.ttime, "--", color='r')
  XXX = sd.ttime[i_min:i_max]
  YYY = powLaw(XXX, popt[0], popt[1]) 
  #axR2.plot(XXX, YYY, ":", color='k')
  xpos = 10**(0.5*(np.log10(XXX[0]) + np.log10(XXX[-1])))
  ypos = 10**(0.4*(np.log10(YYY[0]) + np.log10(YYY[-1])))
  xpos = 1e4
  ypos = 2
  lblExponente = str(popt[0])
  #axR2.text(xpos, ypos, "$ \\sim t^{"+lblExponente+"} "+"$")
  

  #scaXX, scaYY  = ScaleFunction(sd)
  #linea, = axR2Sca.plot(scaXX, scaYY, **estilo)
  return True, linea, ""

###########################################################################

flag, tmpLinea, tmpEtiqueta = plotSimul(axR2, flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials )


if flag:
  print "ok"
  ### lineas   .append(tmpLinea)
  ### etiquetas.append(tmpEtiqueta)
else:
  print strF2B, strB2F, strB2B, strHrad, strBrad, "# MISSING SIMULATION"
    
#legs= axR2.legend(lineas, etiquetas, bbox_to_anchor=(0.8, 0.4), loc="upper left", title="$R_H $", prop={'size':12}, frameon=False, handlelength=2)

figR2.tight_layout()
#figR2.savefig(outfln_PDF)
### print legs
### if legs:
###   axR2.add_artist(legs)
###   legs.get_title().set_fontsize('20')
###   legs.draggable()


# Save Figure
outPrefix    = "MSD_vs_time"+"__B2F_"+strB2F+"__B2B_"+strB2B    #"D_vs_pBB__pBF_"+strB2F #+"__Hrad_"+strHradius
outfln_PDF   = outPrefix+".eps"

figR2.savefig(outfln_PDF)
plt.show()




