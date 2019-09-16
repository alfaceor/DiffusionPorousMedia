#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 
strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
strB2B  = "0.001"
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

def lineStyle_axR2_b(sd):
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


fig, axR2 = plt.subplots(1,2, figsize=(14,6), sharey=True)
fig.subplots_adjust(wspace=0)
#axR2 = ax
axR2_a = axR2[0]
axR2_b = axR2[1]

axR2_a.set_xlabel("$t$", fontsize='30', rotation=0)
axR2_a.set_ylabel("$ r_2  $", fontsize='30')

axR2_a.set_xscale("log")
axR2_a.set_yscale("log")

axR2_a.set_xlim((8e-1,2e7))
axR2_a.set_ylim((8e-1,2e6))
#axR2.set_xlim(1e2,1e6)
#axR2.set_ylim(1e1,1e6)

axR2_a.tick_params(axis='x',which='minor', top='off')
axR2_a.tick_params(axis='x',which='major', top='on')
axR2_a.tick_params(axis='x',which='minor', bottom='off')
axR2_a.tick_params(axis='x',which='major', bottom='on')

axR2_a.tick_params(axis='y',which='minor', right='off')
axR2_a.tick_params(axis='y',which='major', right='on')
axR2_a.tick_params(axis='y',which='minor', left='off')
axR2_a.tick_params(axis='y',which='major', left='on')


lblF2B = latex_float02(float(strF2B))
lblB2B = latex_float02(float(strB2B))
axR2_a.set_title("(a) $ p_{S} = "+lblB2B+" $"+", $R_H = "+strHrad+"$")

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
    axR2_a. plot(sd.ttime, sd.DelR2, **estilo)
    
    if ( getFitParamsFromDB(Directory, sd) ):
      sd.popt, sd.pcov = fitLine(sd.ttime[sd.imin:sd.imax], sd.DelR2[sd.imin:sd.imax])
      sd.DiffCoeff = sd.popt[0]/6.0
      set_D_in_DB(Directory, sd)
      axR2_a. plot(sd.ttime[sd.imin:sd.imax], line(sd.ttime[sd.imin:sd.imax], sd.popt[0], sd.popt[1]), ':', linewidth=2.5, color='k')
      
  except IOError as e:
    pass

axR2_a. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )

handles, labels = axR2_a.get_legend_handles_labels()
leg_a = axR2_a.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{D}$", prop={'size':12}, frameon=False)



################## INPUT PARAMETERS ##################
strF2B = "1.0" 
strB2F = "0.001"
strB2BList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ] 
#strB2B  = "0.000001" #["0.2", "0.4", "0.6", "0.8", "1.0"] 
strHrad = "40"
strBrad = "41"

axR2_b.set_xlabel("$t$", fontsize='30', rotation=0)
#axR2_b.set_ylabel("$ r_2 $", fontsize='30')

axR2_b.set_xscale("log")
axR2_b.set_yscale("log")

axR2_b.set_xlim((8e-1,2e7))
axR2_b.set_ylim((8e-1,2e6))
#axR2_b.set_xlim(1e2,1e6)
#axR2_b.set_ylim(1e1,1e6)

axR2_b.tick_params(axis='x',which='minor', top='off')
axR2_b.tick_params(axis='x',which='major', top='on')
axR2_b.tick_params(axis='x',which='minor', bottom='off')
axR2_b.tick_params(axis='x',which='major', bottom='on')

axR2_b.tick_params(axis='y',which='minor', right='off')
axR2_b.tick_params(axis='y',which='major', right='on')
axR2_b.tick_params(axis='y',which='minor', left='off')
axR2_b.tick_params(axis='y',which='major', left='on')

lblF2B = latex_float02(float(strF2B))
lblB2F = latex_float02(float(strB2F))
axR2_b.set_title("(b) $ p_{D} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")

for b2b in range(len(strB2BList)):
  strB2B = strB2BList[b2b]
  sd = SimulData()
  try:
    sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=False)
    if sd.isDatafile():
      sd.getData()
    else:
      sd.metadata(flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, "1000", loadData=True)
      
    estilo = lineStyle_axR2_b(sd)
    #axR2. plot(sd.ttime, sd.DelR2, linewidth=2, color=plt.cm.Blues(colvals[b2f]), label="$ "+strB2F+" $")
    axR2_b. plot(sd.ttime, sd.DelR2, **estilo)
    
    if ( getFitParamsFromDB(Directory, sd) ):
      sd.popt, sd.pcov = fitLine(sd.ttime[sd.imin:sd.imax], sd.DelR2[sd.imin:sd.imax])
      sd.DiffCoeff = sd.popt[0]/6.0
      set_D_in_DB(Directory, sd)
      axR2_b. plot(sd.ttime[sd.imin:sd.imax], line(sd.ttime[sd.imin:sd.imax], sd.popt[0], sd.popt[1]), ':', linewidth=2.5, color='k')
      
  except IOError as e:
    pass

axR2_b. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )

handles, labels = axR2_b.get_legend_handles_labels()
leg_b = axR2_b.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.02,0.35), loc="lower left", title="$p_{S}$", prop={'size':12}, frameon=False)




fig.savefig("Fig2.pdf")
fig.savefig("msdgeral.eps")

plt.show()
