#!/usr/bin/python
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 
#strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2FList = ["0.00001", "0.0001", "0.001", "0.01"] 
strB2FList = strB2FList [::-1]
#strB2BList = ["0.000001"]  #["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2B  = "0.000001" #["0.2", "0.4", "0.6", "0.8", "1.0"] 
strHradList = ["10", "20", "40", "80", "160"] 

Directory = "RandomIniCond_02/"
flnPrefix = Directory + "rw3D_RegularPorousDepletion"
strNtrials = str(10**5)

  
################## FIGURES CONFIGURATION ##################
###### BEGIN PLOT DECORATION VARIABLES
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 24} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES


figR2, axR2 = plt.subplots(2, 1, figsize=(8,12))
axR2_a = axR2[0]
axR2_b = axR2[1]


def lineStyle_axR2_a(sd):
  lineTypes = ['-', '--', ':', '-.']
  #tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2F       = float(sd.strB2F)
  lblB2F     = latex_float02( pB2F )
  etiqueta   = "$"+lblB2F+" $"
  expBF      = np.abs(np.log10( pB2F ))
  colval     = (expBF+4.0)/10.0
  colmap     = plt.cm.Blues
  
  case = expBF
  if case ==    2:
    tipoLinea  = ':' 
  elif case ==  3:
    tipoLinea  = '-.' 
  elif case ==  4:
    tipoLinea  = '--' 
  elif case ==  5:
    tipoLinea  = '-' 
  else:
    tipoLinea  = '-'
  
  case = strHrad
  if case ==    "10":
    #tipoLinea  = '-' 
    colmap     = plt.cm.Greens
  elif case ==  "20":
    #tipoLinea  = '--' 
    colmap     = plt.cm.Oranges
  elif case ==  "40":
    #tipoLinea  = ':' 
    colmap     = plt.cm.Blues
  elif case ==  "80":
    #tipoLinea  = '-.' 
    colmap     = plt.cm.RdPu
  elif case == "160":
    #tipoLinea  = 'steps--' 
    colmap     = plt.cm.Wistia #gnuplot2
  colorcito  = colmap(colval)
  
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


#######################  Diffusion vs pBF  ############################
# TODO: Plot a figure and save the Difusion coefficient in organized files.

#axR2_a = ax
axR2_a.set_xlabel("$p_{D}t/R$", fontsize='30', rotation=0)
axR2_a.set_ylabel("$ r_2 /(R p_{D} t) $", fontsize='30')

axR2_a.set_xscale("log")
axR2_a.set_yscale("log")

#  axR2_a.set_xlim((8e-1,2e7))
#  axR2_a.set_ylim((8e-1,2e6))
#axR2_b.set_xlim(1e-3,2e4)
axR2_a.set_xlim(3e-5,2e4)
#axR2_a.set_ylim(1e1,1e6)

axR2_a.tick_params(axis='x',which='minor', top='off')
axR2_a.tick_params(axis='x',which='major', top='on')
axR2_a.tick_params(axis='x',which='minor', bottom='off')
axR2_a.tick_params(axis='x',which='major', bottom='on')

axR2_a.tick_params(axis='y',which='minor', right='off')
axR2_a.tick_params(axis='y',which='major', right='on')
axR2_a.tick_params(axis='y',which='minor', left='off')
axR2_a.tick_params(axis='y',which='major', left='on')

#  axR2_a.text(1e1, 1e6,"a)")

lblF2B = latex_float02(float(strF2B))
lblB2B = latex_float02(float(strB2B))
#axR2_a.set_title("(a) $ p_{S} = "+lblB2B+" $") #+", $R_H = "+strHrad+"$")
titulo_a = "(a) $ p_{S} = "+lblB2B+" $"
axR2_a.text(0.08, 0.15, titulo_a , transform=axR2_a.transAxes, verticalalignment='top', fontsize=26)
  
#legends
legendDict = {}

for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  Hrad    = int(strHrad)
  strBrad = str(Hrad+1)
  outPrefix = "D_vs_pBF__pBB_"+strB2B +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"
  legendDict[strHrad] = {}
  
  if strHrad in ["20", "40", "80"] :
    strB2FList = ["0.00001", "0.001"] 
    strB2FList = strB2FList [::-1]
  else:
    strB2FList = ["0.00001", "0.0001", "0.001", "0.01"] 
    strB2FList = strB2FList [::-1]
    
  
  
  lineas    = [None for i in range(len(strB2FList))]
  etiquetas = [None for i in range(len(strB2FList))]
  
  for b2f in range(len(strB2FList)):
    strB2F = strB2FList[b2f]
    flagPlot = True
    if (strHrad == "20") or (strHrad == "40") or (strHrad == "80"):
      if strB2F == "0.01":
        flagPlot = False
    elif (strHrad == "160") :
      if (strB2F == "0.01") or (strB2F == "0.001"):
        flagPlot = False
    pB2F   = float(strB2F)
    sd = SimulData()
    try:
      sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=False)
      if sd.isDatafile():
        sd.getData()
      else:
        sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, "1000", loadData=True)
        
      estilo = lineStyle_axR2_a(sd)
      #axR2_a. plot(sd.ttime, sd.DelR2, linewidth=2, color=plt.cm.Blues(colvals[b2f]), label="$ "+strB2F+" $")
      if flagPlot :
        ttime = sd.ttime[sd.DelR2 > Hrad**2/10.0]
        DelR2 = sd.DelR2[sd.DelR2 > Hrad**2/10.0]
        #lineas[b2f], = axR2_a. plot(sd.ttime*pB2F/Hrad, sd.DelR2/(Hrad*pB2F*sd.ttime), **estilo)
        lineas[b2f], = axR2_a. plot(ttime*pB2F/Hrad, DelR2/(Hrad*pB2F*ttime), **estilo)
        etiquetas[b2f] = estilo['label']
      else:
        estilo['alpha'] = 0.0
        lineas[b2f], = axR2_a. plot(sd.ttime*pB2F/Hrad, sd.DelR2/(Hrad*pB2F*sd.ttime), **estilo)
        #lineas[b2f] = None
        etiquetas[b2f] = ''
      
    except IOError as e:
      pass
    
  #axR2_a. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )
  print lineas, etiquetas
  handles, labels = axR2_a.get_legend_handles_labels()
  
  if strHrad in ["20", "40", "80"] : 
    legendDict[strHrad] = axR2_a.legend(lineas, etiquetas, bbox_to_anchor=(0.75,0.03), loc="upper left", title="$R = "+strHrad+", p_{D}$", prop={'size':20}, frameon=False, numpoints=1, handlelength=0.5)
  else:
    legendDict[strHrad] = axR2_a.legend(lineas, etiquetas, bbox_to_anchor=(0.75,0.03), loc="upper left", title="$R = "+strHrad+", p_{D}$", prop={'size':20}, frameon=False, numpoints=1, ncol=2, handlelength=0.5, columnspacing=1)
  #leg = axR2_a.legend(handles, labels, bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{D}$", prop={'size':12}, frameon=False)
  #figs[b2b].savefig("MSD_vs_t__pBB_"+strB2B +"__Hrad_"+strHrad+".pdf")
  axR2_a.add_artist( legendDict[strHrad] )
  legendDict[strHrad].set_bbox_to_anchor((0.01+0.3*hr,1.0) , transform=None)
  legendDict[strHrad].draggable()

#legendDict[ "10"].set_bbox_to_anchor((0.4,1.0) , transform=None)
#legendDict[ "20"].set_bbox_to_anchor((0.7,1.0) , transform=None)
#legendDict[ "40"].set_bbox_to_anchor((0.4,0.8) , transform=None)
#legendDict[ "80"].set_bbox_to_anchor((0.7,0.8) , transform=None)
#legendDict["160"].set_bbox_to_anchor((0.7,0.6) , transform=None)

#legendDict[ "10"].set_bbox_to_anchor((0.35,1.0) , transform=None)
#legendDict[ "20"].set_bbox_to_anchor((0.67,1.0) , transform=None)
#legendDict[ "40"].set_bbox_to_anchor((0.35,0.8) , transform=None)
#legendDict[ "80"].set_bbox_to_anchor((0.67,0.8) , transform=None)
#legendDict["160"].set_bbox_to_anchor((0.67,0.6) , transform=None)

legendDict[ "10"].set_bbox_to_anchor((0.25,1.0) , transform=None)
legendDict[ "20"].set_bbox_to_anchor((0.67,1.0) , transform=None)
legendDict[ "40"].set_bbox_to_anchor((0.28,0.7) , transform=None)
legendDict[ "80"].set_bbox_to_anchor((0.67,0.7) , transform=None)
legendDict["160"].set_bbox_to_anchor((0.67,0.4) , transform=None)




################################################################
################################################################
################################################################
################################################################






################## INPUT PARAMETERS ##################
strF2B = "1.0" 
#strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2BList = ["0.001", "0.01", "0.1"] #, "1.0"] 
strB2BList = strB2BList [::-1]
#strB2FList = ["0.000001"]  #["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2F     = "0.000001"
#strB2B  = "0.000001" #["0.2", "0.4", "0.6", "0.8", "1.0"] 
strHradList = ["10", "20", "40", "80", "160"] 

Directory = "RandomIniCond_02/"
flnPrefix = Directory + "rw3D_RegularPorousDepletion"
strNtrials = str(10**5)

  
################## FIGURES CONFIGURATION ##################
###### BEGIN PLOT DECORATION VARIABLES
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 24} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES

def lineStyle_axR2_b(sd):
  lineTypes = ['-', '--', ':', '-.']
  #tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2B       = float(sd.strB2B)
  lblB2B     = latex_float02( pB2B )
  etiqueta   = "$"+lblB2B+" $"
  expBB      = np.abs(np.log10( pB2B ))
  colval     = (expBB+4.0)/10.0
  colmap     = plt.cm.Blues
  
  case = expBB
  if case ==    0:
    tipoLinea  = '-'
    anchoLinea = 4
    alfa       = 0.5
  elif case ==    2:
    tipoLinea  = ':' 
  elif case ==  3:
    tipoLinea  = '-.' 
  elif case ==  4:
    tipoLinea  = '--' 
  elif case ==  5:
    tipoLinea  = '-' 
  else:
    tipoLinea  = '-'
    
  
  case = strHrad
  if case ==    "10":
    #tipoLinea  = '-' 
    colmap     = plt.cm.Greens
  elif case ==  "20":
    #tipoLinea  = '--' 
    colmap     = plt.cm.Oranges
  elif case ==  "40":
    #tipoLinea  = ':' 
    colmap     = plt.cm.Blues
  elif case ==  "80":
    #tipoLinea  = '-.' 
    colmap     = plt.cm.RdPu
  elif case == "160":
    #tipoLinea  = 'steps--' 
    colmap     = plt.cm.Wistia #gnuplot2
  colorcito  = colmap(colval)
  
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


#######################  Diffusion vs pBF  ############################
# TODO: Plot a figure and save the Difusion coefficient in organized files.

#figR2, axR2_b = plt.subplots(1,1, figsize=(8,6))
#axR2_b = ax
axR2_b.set_xlabel("$p_{S} t/R^{2} $", fontsize='30', rotation=0)
axR2_b.set_ylabel("$ r_2 /(p_{S} t) $", fontsize='30')

axR2_b.set_xscale("log")
axR2_b.set_yscale("log")

axR2_b.set_xlim(1e-3,2e4)
axR2_b.set_ylim(2e-1,2e2)

axR2_b.tick_params(axis='x',which='minor', top='off')
axR2_b.tick_params(axis='x',which='major', top='on')
axR2_b.tick_params(axis='x',which='minor', bottom='off')
axR2_b.tick_params(axis='x',which='major', bottom='on')

axR2_b.tick_params(axis='y',which='minor', right='off')
axR2_b.tick_params(axis='y',which='major', right='on')
axR2_b.tick_params(axis='y',which='minor', left='off')
axR2_b.tick_params(axis='y',which='major', left='on')

#  axR2_b.text(1e1, 1e6,"a)")

lblF2B = latex_float02(float(strF2B))
lblB2F = latex_float02(float(strB2F))
#axR2_b.set_title("(b) $ p_{D} = "+lblB2F+" $") #+", $R_H = "+strHrad+"$")
titulo_b = "(b) $ p_{D} = "+lblB2F+" $"
axR2_b.text(0.08, 0.15, titulo_b , transform=axR2_b.transAxes, verticalalignment='top', fontsize=26)


#legends
legendDict = {}

for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  Hrad    = int(strHrad)
  strBrad = str(Hrad+1)
  outPrefix = "D_vs_pBB__pBF_"+strB2F +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"
  legendDict[strHrad] = {}
  lineas    = [None for i in range(len(strB2BList))]
  etiquetas = [None for i in range(len(strB2BList))]
  
  for b2b in range(len(strB2BList)):
    strB2B = strB2BList[b2b]
    flagPlot = True
    if (strHrad == "80") or (strHrad == "160"):
      if strB2B == "0.001":
        flagPlot = False
    pB2B   = float(strB2B)
    sd = SimulData()
    try:
      sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=False)
      if sd.isDatafile():
        sd.getData()
      else:
        sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, "1000", loadData=True)
        
      estilo = lineStyle_axR2_b(sd)
      #axR2_b. plot(sd.ttime, sd.DelR2, linewidth=2, color=plt.cm.Blues(colvals[b2f]), label="$ "+strB2F+" $")
      if flagPlot :
        ttime = sd.ttime[sd.DelR2 > Hrad**2/10.0]
        DelR2 = sd.DelR2[sd.DelR2 > Hrad**2/10.0]
        lineas[b2b], = axR2_b. plot(ttime*pB2B/Hrad**2., DelR2/(pB2B*ttime), **estilo)
        etiquetas[b2b] = estilo['label']
      else:
        estilo['alpha'] = 0.0
        lineas[b2b], = axR2_b. plot(sd.ttime*pB2B, sd.DelR2/(pB2B*sd.ttime), **estilo)
        #lineas[b2f] = None
        etiquetas[b2b] = ''
      
    except IOError as e:
      pass
    
  #axR2_b. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )
  print lineas, etiquetas
  handles, labels = axR2_b.get_legend_handles_labels()
  legendDict[strHrad] = axR2_b.legend(lineas, etiquetas, bbox_to_anchor=(0.75,0.03), loc="upper left", title="$R = "+strHrad+", p_{S}$", prop={'size':20}, frameon=False, numpoints=1, ncol=2, handlelength=0.5, columnspacing=1)
  #leg = axR2_b.legend(handles, labels, bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{D}$", prop={'size':12}, frameon=False)
  #figs[b2b].savefig("MSD_vs_t__pBB_"+strB2B +"__Hrad_"+strHrad+".pdf")
  axR2_b.add_artist( legendDict[strHrad] )
  legendDict[strHrad].set_bbox_to_anchor((0.01+0.3*hr,1.0) , transform=None)
  legendDict[strHrad].draggable()

legendDict[ "10"].set_bbox_to_anchor((0.20,1.0) , transform=None)
legendDict[ "20"].set_bbox_to_anchor((0.62,1.0) , transform=None)
legendDict[ "40"].set_bbox_to_anchor((0.28,0.7) , transform=None)
legendDict[ "80"].set_bbox_to_anchor((0.66,0.7) , transform=None)
legendDict["160"].set_bbox_to_anchor((0.66,0.4) , transform=None)

figR2.tight_layout()

figR2.savefig("crossoverplot.eps")
    
#plt.show()



