#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 
strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2BList = ["0.001"]  #["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
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

def lineStyle_axR2(sd):
  tipoLinea  = '-' 
  anchoLinea = 2 
  alfa       = 1.0 
  pB2F       = float(sd.strB2F)
  lblB2F     = latex_float02( pB2F )
  etiqueta   = "$"+lblB2F+" $"
  expBF      = np.abs(np.log10( pB2F ))
  colval     = (expBF+4.0)/10.0
  colmap     = plt.cm.Blues
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
axDiff.set_xlabel("$p_{D}$", fontsize='30')
axDiff.set_ylabel("$D$", fontsize='30', rotation=0)



figs = [None for b2b in range(len(strB2BList))]
for b2b in range(len(strB2BList)):
  strB2B = strB2BList[b2b]
  outPrefix = "D_vs_pBF__pBB_"+strB2B +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"

  figs[b2b], axR2 = plt.subplots(1,1, figsize=(8,6))
  #axR2 = ax
  axR2.set_xlabel("$t$", fontsize='30', rotation=0)
  axR2.set_ylabel("$ r_2  $", fontsize='30')
  
  axR2.set_xscale("log")
  axR2.set_yscale("log")
  
  axR2.set_xlim((8e-1,2e7))
  axR2.set_ylim((8e-1,2e6))
  #axR2.set_xlim(1e2,1e6)
  #axR2.set_ylim(1e1,1e6)
  
  axR2.tick_params(axis='x',which='minor', top='off')
  axR2.tick_params(axis='x',which='major', top='on')
  axR2.tick_params(axis='x',which='minor', bottom='off')
  axR2.tick_params(axis='x',which='major', bottom='on')
  
  axR2.tick_params(axis='y',which='minor', right='off')
  axR2.tick_params(axis='y',which='major', right='on')
  axR2.tick_params(axis='y',which='minor', left='off')
  axR2.tick_params(axis='y',which='major', left='on')
  
#  axR2.text(1e1, 1e6,"a)")
  
  lblF2B = latex_float02(float(strF2B))
  lblB2B = latex_float02(float(strB2B))
  axR2.set_title("(a) $ p_{S} = "+lblB2B+" $"+", $R_H = "+strHrad+"$")
  colvals = np.abs(np.log( map(float, strB2FList)))+3.0
  colvals = (colvals/colvals.max() )[::-1]
  #print colvals
  
  ofile = open(outfln_DAT,"w")
  linea="#prob DiffCoeff\n"
  ofile.write(linea)
  for b2f in range(len(strB2FList)):
    strB2F = strB2FList[b2f]
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
        linea = ""+sd.strB2F +" "+str(sd.DiffCoeff) +"\n"
        ofile.write(linea)
        
    except IOError as e:
      pass
    
  ofile.close()
  
  axR2. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )
  
  handles, labels = axR2.get_legend_handles_labels()
  leg = axR2.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{D}$", prop={'size':12}, frameon=False)
  #figs[b2b].savefig("MSD_vs_t__pBB_"+strB2B +"__Hrad_"+strHrad+".pdf")
  figs[b2b].savefig("Fig2a.pdf")
  figs[b2b].savefig("Fig2a.eps")
  
  # D vs pB2B
  print outfln_DAT
  try:
    pB2F, DiffCoeff = np.loadtxt(outfln_DAT).T
    print pB2F, DiffCoeff
    estiloDiff = lineStyle_axDiff(strB2B)
    axDiff.plot(pB2F, DiffCoeff, **estiloDiff)
  except ValueError as e:
    pass
    

handles, labels = axDiff.get_legend_handles_labels()
axDiff.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.7,0.03), loc="lower left", title="$p_{BB}$", prop={'size':12}, frameon=False, numpoints=1)

figDiff.tight_layout()

#figDiff.savefig("D_vs_pBF__pBB_"+ "VAR" +"__Hrad_"+strHrad+".pdf")

plt.show()




#left, bottom, width, height = [0.6, 0.25, 0.3, 0.25]
#axIns = fig.add_axes([left, bottom, width, height])
#axIns.set_xlim(0, 40)  
#axIns.set_ylim(0, 40)  
#axIns.tick_params(axis='x', labelsize=10)
#axIns.tick_params(axis='y', labelsize=10)  
#axIns.set_xlabel("$t$")
##axIns.set_xscale("log")

#axIns.set_ylabel("$ \\langle \\Delta r^2 \\rangle $")
##axIns.set_yscale("log")

#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#axR2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))



