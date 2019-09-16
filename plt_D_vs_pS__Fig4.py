#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 

strB2FList  = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
#strB2F  = "0.01"

HradList = [10] #[ 10, 20, 40, 80, 160 ]
BradList = [i+1 for i in HradList]
strHradList = map(str, HradList)
strBradList = map(str, BradList)


def ScaleFunction(XX, YY, Hrad):
  scaFactor = (1.0/Hrad)
  scaXX     =  XX*scaFactor
  scaYY     =  YY*scaFactor
  return scaXX, scaYY

lblScaleFunc = "(1/R_H)"


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

def lineStyle_axDiff(strB2F, strHrad):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 0.8 #0.5
  pB2F       = float(strB2F)
  lblB2F     = latex_float02( pB2F )
  etiqueta   = "$"+lblB2F+" $"
  expBF      = np.abs(np.log10( pB2F ))
  colval     = (expBF+4.0)/12.0
  
  if strHrad == "10":
    tipoLinea = '-'
    colmap     = plt.cm.Greens
  elif strHrad == "20":
    colmap     = plt.cm.Oranges
  elif strHrad == "40":
    colmap     = plt.cm.Blues
  elif strHrad == "80":
    colmap     = plt.cm.Purples
  elif strHrad == "160":
    colmap     = plt.cm.YlOrBr
  else:
    colmap     = plt.cm.Reds
  
  colmap = plt.cm.gnuplot2
  colorcito  = colmap(colval)
  
  case = expBF
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
        'markersize'      : 12, #9,
        'markeredgewidth' : 2,
        'markeredgecolor' : colorcito
        } 
  return estilo #'solid'

#######################  Diffusion vs pBF  ############################
# TODO: Plot a figure and save the Difusion coefficient in organized files.
#### D vs p_S FIGURE
figDiff, axDiff = plt.subplots(figsize = (8, 6))
#lblB2F = latex_float(float(strB2F))
#axDiff.set_title("$ p_{D} = "+lblB2F+" $ ")#+", $ p_{BF} = "+strB2F+" $")
#axDiff.set_title("$ R_H = "+strHrad+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiff.set_xscale("log")
axDiff.set_yscale("log")
axDiff.set_xlabel("$p_{S}$", fontsize='30')
axDiff.set_ylabel("$D$", fontsize='30', rotation=0)

axDiff.set_ylim((1e-7, 4e-1))

axDiff.tick_params(axis='x', which='minor', bottom='off')
axDiff.tick_params(axis='x', which='major', top='on')

#axDiff.tick_params(axis='y', which='minor', right='off')
axDiff.tick_params(axis='y', which='major', right='on')


#### SCALED FIGURE
figScaDiff, axScaDiff = plt.subplots(figsize = (8, 6))
#lblB2F = latex_float(float(strB2F))
#axScaDiff.set_title("$ p_{D} = "+lblB2F+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axScaDiff.set_xscale("log")
axScaDiff.set_yscale("log")
axScaDiff.set_xlabel("$p_{S}/(p_{D} R_H)$", fontsize='30')
axScaDiff.set_ylabel("$D/(p_{D} R_H)$",     fontsize='30') #, rotation=0)

axScaDiff.set_xlim((2e-9, 2e5))
axScaDiff.set_ylim((1e-3, 1e4))

axScaDiff.tick_params(axis='x', which='minor', bottom='off')
axScaDiff.tick_params(axis='x', which='major', top='on')

#axScaDiff.tick_params(axis='y', which='minor', right='off')
axScaDiff.tick_params(axis='y', which='major', right='on')

# legends
legendDict = {} #[None for i in range(len(strHradList))]

inpPrefix  = "D_vs_pBB" #+"__pBB_"+strB2B +"__Hrad_"+strHrad
outPrefix  = "D_vs_pS"+"__pD_"+"VAR" #+"__Hrad_"+strHrad
for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  strBrad = strBradList[hr]
  Hrad    = HradList   [hr]
  lineas     = [None for i in range(len(strB2FList))]
  etiquetas  = [None for i in range(len(strB2FList))]
  #legendDict[strHrad] = {'lineas' : [None for i in range(len(strB2FList))], 'etiquetas' : [None for i in range(len(strB2FList))]}
  for b2f in range(len(strB2FList)):
    strB2F = strB2FList[b2f]
    pB2F   = float(strB2F)
    inpPrefix  = "D_vs_pBB"
    inpfln_DAT = inpPrefix+"__pBF_"+strB2F +"__Hrad_"+strHrad+".DAT"  
    # D vs pB2B
    pB2B, DiffCoeff = np.loadtxt(inpfln_DAT).T
    estiloDiff      = lineStyle_axDiff(strB2F, strHrad)
    #legendDict[strHrad]['lineas']   [b2f], = axDiff.plot(pB2B, DiffCoeff, **estiloDiff)
    #legendDict[strHrad]['etiquetas'][b2f]  = estiloDiff['label']
    lineas   [b2f], = axDiff.plot(pB2B, DiffCoeff, **estiloDiff)
    etiquetas[b2f]  = estiloDiff['label']
    sca_pB2B, sca_DiffCoeff = ScaleFunction(pB2B/pB2F, DiffCoeff/pB2F, Hrad)
    axScaDiff.plot(sca_pB2B, sca_DiffCoeff, **estiloDiff)
    
  # legendDict[strHrad] = {'lineas' : lineas, 'etiquetas' : etiquetas}
  # handles, labels = axDiff.get_legend_handles_labels()
  #legendDict[strHrad]['legenda'] = axScaDiff.legend(legendDict[strHrad]['lineas'], legendDict[strHrad]['etiquetas'], bbox_to_anchor=(0.1+float(hr)/5.0,0.03), loc="lower left", title="$ R_H = "+strHrad+", p_{D}$", prop={'size':12}, frameon=False, numpoints=1)
  legendDict[strHrad] = axScaDiff.legend(lineas[::-1], etiquetas[::-1], bbox_to_anchor=(0.6*hr/float(len(strHradList)), 1.0), loc="upper left", title="$ R_H = "+strHrad+", p_{D}$", prop={'size':12}, frameon=False, numpoints=1)
#  print legendDict
#  for a in legendDict[strHrad].get_texts():
#    a.set_ha("right")
#    a.set_position((20,0))


for shr in legendDict:
  print shr
  print 20*'+'
  axScaDiff.add_artist( legendDict[shr] )
  #legendDict[shr].set_bbox_to_anchor((0.01,1.0) , transform=None)
  legendDict[shr].draggable()

# outfln_PDF = outfln_PDF +"__Hrad_"+strHrad
# figDiff.tight_layout()
# outfln_PDF = outPrefix+".pdf"
# figDiff.savefig(outfln_PDF)

#handles, labels = axScaDiff.get_legend_handles_labels()
#legenda = axScaDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.8,0.03), loc="lower left", title="$ R_H$", prop={'size':12}, frameon=False, numpoints=1)
#for a in legenda.get_texts():
#  a.set_ha("right")
#  a.set_position((20,0))

figScaDiff.tight_layout()
outfln_PDF = outPrefix+"_sca.pdf"
print outfln_PDF
figScaDiff.savefig(outfln_PDF)

plt.show()
