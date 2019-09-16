#!/usr/bin/python
import numpy as np
from SimulData import *

#flnPrefix="CenterIniCond/rw3D_RegularPorousDepletion"
flnPrefix="RandomIniCond_02/rw3D_RegularPorousDepletion"
strNtrials = str(10**5)


def ScaleFunction(sd):
  scaXX = sd.ttime/sd.Hrad**2.0
  scaYY = sd.DelR2/(sd.ttime)
  return scaXX, scaYY


strF2B  = "1.0" 
strB2FList  = ["1.0", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001"] #"0.001"
strB2B  = "0.001" #"0.000001"
strHrad = "10" 
strBrad = "11" 

lblF2B = latex_float( float(strF2B) )
#lblB2F = latex_float( float(strB2F) )
lblB2B = latex_float( float(strB2B) )

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
axR2.set_title("$p_{FB} = "+lblF2B+
               #", p_{BF} ="+lblB2F+
               ", p_{BB} ="+lblB2B+
               ", R_H = "+strHrad+
               " $")
#axR2Sca.set_title("$p_{FB} = "+lblF2B+", p_{BB} \\gg p_{BF} ="+lblB2F+" $")
axR2.set_xscale("log")
axR2.set_yscale("log")

#axR2.set_xlim(1e0, 2e3)
#axR2.set_ylim(1e0, 2e3)


axR2.set_xlabel("$t$", fontsize='30', rotation=0)
axR2.set_ylabel("$ \\langle \\Delta r^2 \\rangle $", fontsize='30')

#axR2Sca.tick_params(axis='y',which='minor',left='off')
axR2.tick_params(axis='x',which='minor',bottom='off')
axR2.tick_params(axis='y',which='minor',right='on')
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
  pB2F       = float(sd.strB2F)
  lblB2F     = latex_float( pB2F )
  etiqueta   = "$"+lblB2F+" $"
  expBF      = np.abs(np.log10( pB2F ))

  #colval = (expBF + 10.*expBB) / 36.
  if sd.strHrad == "10":
    colmap = plt.cm.Greens_r
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
  colval = (expBF+4.0)/12.0 #+3.0
  #print "colval = ",colval
  colorcito  = colmap(colval) #0.2+0.8*expBB)

#  if expBF == 0 : 
#    tipoLinea = '-'
#  elif expBF == 1 : 
#    tipoLinea = '--'
#  elif expBF == 2 : 
#    tipoLinea = '-.' 
#  elif expBF == 3 : 
#    tipoLinea = ':' 
#  elif expBF == 4 : 
#    tipoLinea = '-.' 
#  else:
#    tipoLinea = ':' 
  tipoLinea = '-'
  
  
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
  axR2.plot(sd.ttime, sd.ttime, "--")
  axR2.plot(sd.ttime, 1.5*sd.ttime, "--")
  #scaXX, scaYY  = ScaleFunction(sd)
  #linea, = axR2Sca.plot(scaXX, scaYY, **estilo)
  return True, linea, ""

###########################################################################

#flag, tmpLinea, tmpEtiqueta = plotSimul(axR2, flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials )
#if flag:
#  print "ok"
#  ### lineas   .append(tmpLinea)
#  ### etiquetas.append(tmpEtiqueta)
#else:
#  print strF2B, strB2F, strB2B, strHrad, strBrad, "# MISSING SIMULATION"

for b2f in range(len(strB2FList)):
  strB2F = strB2FList[b2f]
  lblB2F = latex_float( float(strB2F) )
  sd = SimulData()
  sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=True )
  estilo = lineStyDAT(sd)
  axR2.plot(sd.ttime, sd.DelR2, **estilo)

linea, = axR2.plot(sd.ttime, sd.ttime, "--", color='r')
drD = DraggableLogLogLine(linea)
drD.connect()

legs = axR2.legend( bbox_to_anchor=(0.8, 0.4), loc="upper left", title="$p_{BF} $", prop={'size':12}, frameon=False, handlelength=2)
#legs= axR2.legend(lineas, etiquetas, bbox_to_anchor=(0.8, 0.4), loc="upper left", title="$R_H $", prop={'size':12}, frameon=False, handlelength=2)

figR2.tight_layout()
#figR2.savefig(outfln_PDF)
### print legs
if legs:
  axR2.add_artist(legs)
  legs.get_title().set_fontsize('20')
  legs.draggable()


# Save Figure
outPrefix    = "MSD_vs_time"+"__B2F_"+strB2F+"__B2B_"+strB2B    #"D_vs_pBB__pBF_"+strB2F #+"__Hrad_"+strHradius
outfln_PDF   = outPrefix+".pdf"

figR2.savefig(outfln_PDF)
plt.show()




