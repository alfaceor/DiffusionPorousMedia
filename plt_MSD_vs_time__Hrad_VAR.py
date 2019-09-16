#!/usr/bin/python
import numpy as np
from SimulData import *

flnPrefix="rw3D_RegularPorousDepletion"
strNtrials = str(10**5)
dirList=["longTime/longTime2/", "longTime/", "./", "../", "../rw3D_RegularPorousMODIFIED/"]


def ScaleFunction(sd):
  scaXX = sd.ttime/sd.Hrad**2.0
  scaYY = sd.DelR2/(sd.ttime)
  return scaXX, scaYY


strF2B      = "1.0" 
#strB2F      = "0.000001"
strB2F      = "1.0"
#strB2BList  = ["0.0001", "0.01"] #, "0.1", "1.0"]
strB2B      = "1.0"

strHradList = ["10", "20", "40" ] #, "80", "160" ] #, "500" ]
strBradList = ["11", "21", "41" ] #, "81", "161" ] #, "501" ]

DiffCoeff  = np.zeros(len(strHradList))

lblF2B = latex_float( float(strF2B) )
lblB2F = latex_float( float(strB2F) )
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

outPrefix    = "MSD_vs_time"+"__B2F_"+strB2F+"__B2B_"+strB2B    #"D_vs_pBB__pBF_"+strB2F #+"__Hrad_"+strHradius
#fitParamsFln = outPrefix+".json"
outfln_PDF   = outPrefix+".pdf"
#outfln_DAT   = outPrefix+".DAT"

#fig, ax = plt.subplots(2,2, sharex=True)
figR2, ax = plt.subplots(1,1, figsize=(8,6))

axR2    = ax
axR2.set_title("$p_{FB} = "+lblF2B+", p_{BF} ="+lblB2F+", p_{BB} ="+lblB2B+" $")
#axR2Sca.set_title("$p_{FB} = "+lblF2B+", p_{BB} \\gg p_{BF} ="+lblB2F+" $")
axR2.set_xscale("log")
axR2.set_yscale("log")

axR2.set_xlim(1e0, 2e3)
axR2.set_ylim(1e0, 2e3)


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


#pB2B = map(float,strB2B)
DiffCoeff = -1.0*np.ones(len(strHradList))
D_err_mim = -1.0*np.ones(len(strHradList))
D_err_max = -1.0*np.ones(len(strHradList))
Figuras   =  [None for i in range(len(strHradList))]
Conexions =  [None for i in range(len(strHradList))]

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

  if expBB == 0 : 
    tipoLinea = '-'
  elif expBB == 1 : 
    tipoLinea = '--'
  elif expBB == 2 : 
    tipoLinea = '-.' 
  elif expBB == 3 : 
    tipoLinea = ':' 
  elif expBB == 4 : 
    tipoLinea = '-.' 
  else:
    tipoLinea = ':' 
  
  
  
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


def plotSimul(axR2, axR2Sca, dirList, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials ):
  lblB2B = latex_float( float(strB2B) )
  flag = False
  idir = 0
  while idir < len(dirList):
    try:
      sd = SimulData()
      #print dirList[idir]
      flnPrefix02 = dirList[idir]+flnPrefix 
      sd.metadata(flnPrefix02, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials )
      sd.getData()
      estilo = lineStyDAT(sd)
#      titulo =  "$p_{FB} = "+lblF2B+", p_{BF} = "+lblB2F+", p_{BB} = "+lblB2B+", R_H = "+strHrad+"$"
#      if fit_params[strHrad]["lambFlag"]:
#        Figuras[hr], Conexions[hr] = sd.calcDwithLamb(fit_params[strHrad], colorcito, titulo)
#      else:
#        sd.calcDiffCoeff(fit_params[strHrad])
#      DiffCoeff[hr] = sd.DiffCoeff
#      D_err_mim[hr] = sd.D_err_min
#      D_err_max[hr] = sd.D_err_max
      axR2.plot(sd.ttime, sd.DelR2, **estilo)
      #axR2.plot(sd.ttime[sd.imin:sd.imax], line(sd.ttime[sd.imin:sd.imax], sd.popt[0], sd.popt[1]), '--', color='k' )
      
      scaXX, scaYY  = ScaleFunction(sd)
      #linea, = axR2Sca.plot(sd.pB2B*sd.ttime/sd.Hrad**2.0, sd.DelR2/(sd.pB2B*sd.ttime), **estilo)
      linea, = axR2Sca.plot(scaXX, scaYY, **estilo)
      flag = True
      return flag, linea, estilo['label']
      break;
    except IOError as e:
      flag = False
      #print e
      pass
    idir = idir + 1
  
  return False, None, None


###########################################################################

lineas    = [] #[None for i in range(len(strB2FList))]
etiquetas = [] #[None for i in range(len(strB2FList))]
for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  strBrad = strBradList[hr]
  flag, tmpLinea, tmpEtiqueta = plotSimul(axR2, axR2Sca, dirList, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials )
  if flag:
    lineas   .append(tmpLinea)
    etiquetas.append(tmpEtiqueta)
  else:
    print strF2B, strB2F, strB2B, strHrad, strBrad, "# MISSING SIMULATION"
    
legs= axR2.legend(lineas, etiquetas, bbox_to_anchor=(0.8, 0.4), loc="upper left", title="$R_H $", prop={'size':12}, frameon=False, handlelength=2)

figR2.tight_layout()
#figR2.savefig(outfln_PDF)
print legs
if legs:
  axR2.add_artist(legs)
  legs.get_title().set_fontsize('20')
  legs.draggable()


figR2.savefig(outfln_PDF)
plt.show()




