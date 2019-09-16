#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 

#strB2BList  = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2B  = "0.000001"

HradList = [ 10, 20, 40, 80, 160 ]
BradList = [i+1 for i in HradList]
strHradList = map(str, HradList)
strBradList = map(str, BradList)


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

def lineStyle_axDiff(strB2B, strHrad):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2B       = float(strB2B)
  lblB2B     = latex_float( pB2B )
  etiqueta   = "$"+strHrad+" $"
  expBB      = np.abs(np.log10( pB2B ))
  colval     = (expBB+4.0)/10.0
  
  if strHrad == "10":
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

figDiff, axDiff = plt.subplots(figsize = (8, 6))
lblB2B = latex_float(float(strB2B))
axDiff.set_title("$ p_{S} = "+lblB2B+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiff.set_xscale("log")
axDiff.set_yscale("log")
axDiff.set_xlabel("$p_{D}$", fontsize='30')
axDiff.set_ylabel("$D$", fontsize='30', rotation=0)
#axDiff.set_ylim((1e-7, 4e-1))

axDiff.tick_params(axis='x', which='minor', bottom='off')
axDiff.tick_params(axis='x', which='major', top='on')

#axDiff.tick_params(axis='y', which='minor', right='off')
axDiff.tick_params(axis='y', which='major', right='on')

inpPrefix  = "D_vs_pBF" #+"__pBB_"+strB2B +"__Hrad_"+strHrad
outPrefix  = "D_vs_pD"+"__pS_"+strB2B+"__Hrad_"+"VAR"

for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  strBrad = strBradList[hr]
  #strB2B = strB2BList[b2b]
  pB2B   = float(strB2B)
  inpPrefix  = "D_vs_pBF"
  
  inpfln_DAT = inpPrefix+"__pBB_"+strB2B +"__Hrad_"+strHrad+".DAT"  
  # D vs pB2B
  pB2F, DiffCoeff = np.loadtxt(inpfln_DAT).T
  estiloDiff = lineStyle_axDiff(strB2B, strHrad)
  #axDiff.plot(pB2F, DiffCoeff, marker="o", markersize=9, linewidth=2, markeredgewidth=2, markeredgecolor=colorcito, color=colorcito, label="$ "+strB2B +" $")
  axDiff.plot(pB2F, DiffCoeff, **estiloDiff)

handles, labels = axDiff.get_legend_handles_labels()
legenda = axDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.7,0.03), loc="lower left", title="$ R_H$", prop={'size':12}, frameon=False, numpoints=1)
for a in legenda.get_texts():
  a.set_ha("right")
  a.set_position((20,0))

figDiff.tight_layout()
outfln_PDF = outPrefix+".pdf"
figDiff.savefig(outfln_PDF)

plt.show()
