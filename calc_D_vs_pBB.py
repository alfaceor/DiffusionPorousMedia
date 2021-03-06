#!/usr/bin/python
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 
strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2BList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
#strB2B  = "0.000001" #["0.2", "0.4", "0.6", "0.8", "1.0"] 
strHrad = "10"
strBrad = "11"
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

def lineStyle_axR2(sd):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2B       = float(sd.strB2B)
  lblB2B     = latex_float02( pB2B )
  etiqueta   = "$"+lblB2B+" $"
  expBB      = np.abs(np.log10( pB2B ))
  colval     = (expBB+4.0)/10.0
  colmap     = plt.cm.Reds
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

figDiff, axDiff = plt.subplots(figsize=(8,6))
axDiff.set_title("$ R_H = "+strHrad+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiff.set_xscale("log")
axDiff.set_yscale("log")
axDiff.set_xlabel("$p_{BB}$", fontsize='30')
axDiff.set_ylabel("$D$", fontsize='30', rotation=0)
#axDiff.set_ylim((0,0.2))

figs = [None for b2f in range(len(strB2FList))]
figsNames = [None for b2f in range(len(strB2FList))]
for b2f in range(len(strB2FList)):
  strB2F = strB2FList[b2f]
  outPrefix = "D_vs_pBB__pBF_"+strB2F +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"

  figs[b2f], axR2 = plt.subplots(1,1, figsize=(8,6))
  #axR2 = ax
  axR2.set_xlabel("$t$", fontsize='30', rotation=0)
  axR2.set_ylabel("$ \\langle \\Delta r^2 \\rangle $", fontsize='30')
  
  axR2.set_xscale("log")
  axR2.set_yscale("log")
  
  #axR2.set_xlim(1e2,1e6)
  #axR2.set_ylim(1e1,1e6)
  
  lblF2B = latex_float02(float(strF2B))
  lblB2F = latex_float02(float(strB2F))
  axR2.set_title("$ p_{FB} = "+lblF2B+" $, $ p_{BF} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")
  colvals = np.abs(np.log( map(float, strB2BList)))+3.0
  colvals = (colvals/colvals.max() )[::-1]
  #print colvals
  
  ofile = open(outfln_DAT,"w")
  linea="#prob DiffCoeff\n"
  ofile.write(linea)
  for b2b in range(len(strB2BList)):
    strB2B = strB2BList[b2b]
    sd = SimulData()
    try:
      sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=False)
      if sd.isDatafile():
        sd.getData()
      else:
        sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, "1000", loadData=True)
        
      estilo = lineStyle_axR2(sd)
      #axR2. plot(sd.ttime, sd.DelR2, linewidth=2, color=plt.cm.Blues(colvals[b2f]), label="$ "+strB2F+" $")
      axR2. plot(sd.ttime, sd.DelR2, **estilo)
      
      if ( getFitParamsFromDB(Directory, sd) ):
        sd.popt, sd.pcov = fitLine(sd.ttime[sd.imin:sd.imax], sd.DelR2[sd.imin:sd.imax])
        sd.DiffCoeff = sd.popt[0]/6.0
        set_D_in_DB(Directory, sd)
        axR2. plot(sd.ttime[sd.imin:sd.imax], line(sd.ttime[sd.imin:sd.imax], sd.popt[0], sd.popt[1]), ':', linewidth=2.5, color='k')
        linea = ""+sd.strB2B +" "+str(sd.DiffCoeff) +"\n"
        ofile.write(linea)
        
    except IOError as e:
      pass
    
  ofile.close()
  
  axR2. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )
  
  handles, labels = axR2.get_legend_handles_labels()
  leg = axR2.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{BB}$", prop={'size':12}, frameon=False)
  filenameMSD = "MSD_vs_t__pBF_"+strB2F +"__Hrad_"+strHrad+".pdf"
  figs[b2f].savefig(filenameMSD)
  figsNames[b2f] = filenameMSD
  
  # D vs pB2B
  print outfln_DAT
  try:
    pB2B, DiffCoeff = np.loadtxt(outfln_DAT).T
    print pB2B, DiffCoeff
    estiloDiff = lineStyle_axDiff(strB2F)
    axDiff.plot(pB2B, DiffCoeff, **estiloDiff)
  except ValueError as e:
    pass
    

handles, labels = axDiff.get_legend_handles_labels()
axDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.7,0.03), loc="lower left", title="$p_{BF}$", prop={'size':12}, frameon=False, numpoints=1)

figDiff.tight_layout()

filenameDiff="D_vs_pBB__pBF_"+ "VAR" +"__Hrad_"+strHrad+".pdf"
figDiff.savefig(filenameDiff)

#plt.show()

#### Gerar Reporte FINAL ####
import subprocess

outputReport = "D_vs_pBB__Hrad_"+strHrad+".pdf"
figsNames[b2f] = filenameMSD
auxFigs = ""
for i in range(len(figsNames)):
  auxFigs = auxFigs + " " + figsNames[i]
cmd = "pdftk "+filenameDiff+" "+auxFigs+" cat output "+outputReport
subprocess.call(cmd, shell=True)





