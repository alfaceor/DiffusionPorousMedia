#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 

strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01"] #, "0.1" ] #, "0.1", "1.0" ]
strB2FList = strB2FList[::-1]
#strB2BList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01" ] #, "0.1", "1.0" ]
#strB2BList = strB2BList[::-1]

D0 = 0.25

#strB2B  = "0.000001" #["0.2", "0.4", "0.6", "0.8", "1.0"] 

HradList = [10, 20, 40, 80, 160]
strHradList = map(str, HradList)
strBradList = [str(i+1) for i in HradList]
#strHrad = "10"
#strBrad = "11"
Directory = "RandomIniCond_02/"
flnPrefix = Directory + "rw3D_RegularPorousDepletion"
strNtrials = str(10**5)

################## SQL FUNCTIONS ##################
import sqlite3

def getFitParamsFromDB(Dir, sd):
  db = sqlite3.connect("SimulData.db")
  cursor = db.cursor()
  sqldata = cursor.execute('''SELECT FitMin, FitMax FROM Simulation WHERE Directory = ? AND strF2B = ? AND strB2F = ? AND strB2B = ? AND strHrad = ? AND strBrad =?''', (Dir, sd.strF2B, sd.strB2F, sd.strB2B, sd.strHrad, sd.strBrad))
  row = sqldata.fetchone()
  flag = False
  if (row != None ):
    FitMin, FitMax = row
    #print row
    sd.imin, sd.imax = getLimits(sd.ttime, FitMin, FitMax)
    flag = True
  db.close()
  return flag
  
def set_D_in_DB(Dir, sd):
  db = sqlite3.connect("SimulData.db")
  cursor = db.cursor()
  sqldata = cursor.execute('''UPDATE Simulation SET D = ? WHERE Directory = ? AND strF2B = ? AND strB2F = ? AND strB2B = ? AND strHrad = ? AND strBrad =?''', (sd.DiffCoeff, Dir, sd.strF2B, sd.strB2F, sd.strB2B, sd.strHrad, sd.strBrad))
  db.commit()
  db.close()
  
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
  tipoLinea = "none" 
  anchoLinea = 2 
  alfa = 1.0 
  
  pB2F   = float(strB2F)
  lblB2F = latex_float02( pB2F )
  etiqueta = "$"+lblB2F+" $"
  marcFaceCol = "none"
  
  expBF  = np.abs(np.log10( pB2F ))
  if expBF == 0 : 
    marcador = "o"
  elif expBF == 1 : 
    marcador = "X"
  elif expBF == 2 : 
    marcador = "s" 
  else:
    marcador = (int(expBF), 0, 0)
  print marcador

  #colval = (expBF + 10.*expBB) / 36.
  if strHrad == "10":
    colmap = plt.cm.Greens
    #tipoLinea = '-'
  elif strHrad == "20":
    colmap = plt.cm.Oranges
  elif strHrad == "40":
    colmap = plt.cm.Blues
  elif strHrad == "80":
    colmap = plt.cm.RdPu
  elif strHrad == "160": 
    colmap = plt.cm.Wistia # Yellows?
  else:
    colmap = plt.cm.Greens
  #colval  = 1.0 - expBF/6.0
  colval = (expBF+4.0)/10.0 #+3.0
  #print "colval = ",colval
  colorcito  = colmap(colval) #0.2+0.8*expBB)
  marcEdgeCol  = colorcito
  
  if strHrad == "80":
    marcFaceCol = colorcito
    marcEdgeCol = 'k'
    alfa =0.5
  elif strHrad == "160":
    marcFaceCol = colorcito
    marcEdgeCol = 'k'
    alfa =0.75
  
  
  estilo = {'linestyle'   : tipoLinea, # (0, (5, 10)), 
        'color'           : colorcito,
        'marker'          : marcador,
        'markersize'      : 12,
        'markeredgewidth' : 2,
        'markerfacecolor' : marcFaceCol,
        'markeredgecolor' : marcEdgeCol,
        'linewidth'       : anchoLinea,
        'label'           : etiqueta,
        'alpha'           : alfa} 
  return estilo #'solid'





def lineStyDAT(strB2B, strHrad):
  tipoLinea  = 'None' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2B       = float(strB2B)
  lblB2B     = latex_float02( pB2B )
  etiqueta   = "$"+lblB2B+" $"
  expBB      = np.abs(np.log10( pB2B ))
  colval     = (expBB+4.0)/10.0
  
  if   strHrad ==  "10":
    #tipoLinea  = '-' 
    colmap   = plt.cm.Blues
  elif strHrad ==  "20":
    colmap   = plt.cm.Greens
  elif strHrad ==  "40":
    colmap   = plt.cm.Oranges
  elif strHrad ==  "80":
    colmap   = plt.cm.Reds
  elif strHrad == "160":
    colmap   = plt.cm.Greys
  else:
    colmap   = plt.cm.Purples
  
  colorcito  = colmap(colval)
  
  #case = expBB
  case = np.abs(np.log2( float(strHrad)/10 ))
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
        'markersize'      : 10,
        'markeredgewidth' : 2,
        'markeredgecolor' : colorcito,
        'markerfacecolor' : 'None'
        } 
  return estilo #'solid'


#######################  Diffusion vs pBF  ############################
# TODO: Plot a figure and save the Difusion coefficient in organized files.
fig, axDiff = plt.subplots(1,1, figsize=(8,6))

xlbl = "$p_{S}/(R p_{D})$"
ylbl = "$ D/(R p_{D} D_{0})$"
axDiff.set_xlabel(xlbl, fontsize='30', rotation=0)
axDiff.set_ylabel(ylbl, fontsize='30')

axDiff.set_xscale("log")
axDiff.set_yscale("log")

#axDiff.set_xlim((1e-4,1e6))
#axDiff.set_ylim((8e-1,2e6))


lblF2B = latex_float02(float(strF2B))
#lblB2F = latex_float02(float(strB2F))
#axR2.set_title("$ p_{D} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")
#axDiff.text(2e-4, 5e-1, "$ p_{D} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")


axDiff.tick_params(axis='x',which='minor', top='off')
axDiff.tick_params(axis='x',which='major', top='on')
axDiff.tick_params(axis='x',which='minor', bottom='off')
axDiff.tick_params(axis='x',which='major', bottom='on')

axDiff.tick_params(axis='y',which='minor', right='off')
axDiff.tick_params(axis='y',which='major', right='on')
axDiff.tick_params(axis='y',which='minor', left='off')
axDiff.tick_params(axis='y',which='major', left='on')

#### Inset Figure
left, bottom, width, height = [0.26, 0.43, 0.38, 0.33]
axDiffInset = fig.add_axes([left, bottom, width, height])
#axDiffInset.set_title("$ R_H = "+strHrad+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiffInset.set_xscale("log")
#axDiffInset.set_yscale("log")
axDiffInset.set_xlabel(xlbl) #, fontsize='30')
axDiffInset.set_ylabel(ylbl) #, rotation=0) #, fontsize='30')

axDiffInset.set_xlim((1e-7,1e0))
axDiffInset.set_ylim((0.1, 0.3 )) #1e-1,1e0))
axDiffInset.set_yticks([0.1, 0.2, 0.3]) #1e-1,1e0))

axDiffInset.tick_params(axis='x', labelsize=10)
axDiffInset.tick_params(axis='y', labelsize=10)  

axDiffInset.tick_params(axis='both', which='major', left='on', right='on', top='on', bottom='on')
axDiffInset.tick_params(axis='both', which='minor', left='off', right='off', top='off', bottom='off')


# legends
legendDict = {}

for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  strBrad = strBradList[hr]
  print strHrad, strBrad
  legendDict[strHrad] = {}
  lineas    = [None for i in range(len(strB2FList))]
  etiquetas = [None for i in range(len(strB2FList))]
  for b2f in range(len(strB2FList)):
    strB2F  = strB2FList[b2f]
    pB2F    = float(strB2F)
    Hrad    = float(strHrad)
    outPrefix = "D_vs_pBB__pBF_"+strB2F +"__Hrad_"+strHrad
    outfln_PDF = outPrefix+".pdf"
    outfln_DAT = outPrefix+".DAT"  
    # D vs pB2B
    #print outfln_DAT
    try:
      pB2B, DiffCoeff = np.loadtxt(outfln_DAT).T
      #print pB2B, DiffCoeff
      estiloDiff = lineStyle_axDiff(strB2F, strHrad)
      scaXX = pB2B / (Hrad*pB2F)
      scaYY = DiffCoeff / (Hrad*pB2F*D0)
      lineas[b2f],   = axDiff     .plot(scaXX, scaYY, **estiloDiff)
      etiquetas[b2f] = estiloDiff['label']
      if axDiffInset:
        axDiffInset.plot(scaXX, scaYY, **estiloDiff)
      
    except ValueError as e:
      pass
        
  
  legendDict[strHrad] = axDiff.legend(lineas, etiquetas, bbox_to_anchor=(0.75,0.03), loc="upper left", title="$R = "+strHrad+", p_{D}$", prop={'size':10}, frameon=False, numpoints=1, ncol=2)
#  handles, labels     = axDiff.get_legend_handles_labels()
#  legendDict[strHrad] = axDiff.legend(handles, labels, bbox_to_anchor=(0.75,0.03), loc="lower left", title="$R = "+strHrad+", p_{D}$", prop={'size':10}, frameon=False, numpoints=1)

#ii = 0
#for shr in legendDict:
  print strHrad
  print 20*'+'
  axDiff.add_artist( legendDict[strHrad] )
  #legendDict[shr].set_bbox_to_anchor=((0.1,0.5), transform=None)
  legendDict[strHrad].set_bbox_to_anchor((0.01+0.3*hr,1.0) , transform=None)
  legendDict[strHrad].draggable()
  

legendDict[ "80"].set_bbox_to_anchor((0.72, 0.5) , transform=None)
legendDict["160"].set_bbox_to_anchor((0.72, 0.25) , transform=None)

fig.tight_layout()

#figDiff.savefig("D_vs_pBB__pBF_"+ "VAR" +"__Hrad_"+strHrad+".pdf")
fig.savefig("Fig4_paper.pdf")
fig.savefig("Fig4_paper.eps")
fig.savefig("Dgeral.eps")

plt.show()

