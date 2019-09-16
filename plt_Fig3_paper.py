#!/usr/bin/python
#import SimulResult as sr
import numpy as np
from SimulData import *

################## INPUT PARAMETERS ##################
strF2B = "1.0" 
strB2F2Plot = "0.000001"

strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ]
strB2FList = strB2FList[::-1]
strB2BList = ["0.0001", "0.001", "0.01", "0.1", "1.0" ] # ["0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0" ]
strB2BList = strB2BList[::-1]
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
strB2F = strB2F2Plot
fig, axR2 = plt.subplots(1,1, figsize=(8,6))
#axR2 = ax

axR2.set_xlabel("$p_S t$", fontsize='30', rotation=0)
axR2.set_ylabel("$ r_2 / (p_S \\times t)$", fontsize='30')

axR2.set_xscale("log")
axR2.set_yscale("log")

axR2.set_xlim((1e-4,1e6))
#axR2.set_ylim((8e-1,2e6))


lblF2B = latex_float02(float(strF2B))
lblB2F = latex_float02(float(strB2F))
#axR2.set_title("$ p_{D} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")
axR2.text(2e-4, 5e-1, "$ p_{D} = "+lblB2F+" $"+", $R_H = "+strHrad+"$")

axR2.tick_params(axis='x',which='minor', top='off')
axR2.tick_params(axis='x',which='major', top='on')
axR2.tick_params(axis='x',which='minor', bottom='off')
axR2.tick_params(axis='x',which='major', bottom='on')

axR2.tick_params(axis='y',which='minor', right='off')
axR2.tick_params(axis='y',which='major', right='on')
axR2.tick_params(axis='y',which='minor', left='off')
axR2.tick_params(axis='y',which='major', left='on')



colvals = np.abs(np.log( map(float, strB2BList)))+3.0
colvals = (colvals/colvals.max() )[::-1]
#print colvals

#ofile = open(outfln_DAT,"w")
#linea="#prob DiffCoeff\n"
#ofile.write(linea)

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
    axR2. plot(sd.pB2B * sd.ttime, sd.DelR2/(sd.pB2B * sd.ttime), **estilo)
    
#        if ( getFitParamsFromDB(Directory, sd) ):
#          sd.popt, sd.pcov = fitLine(sd.ttime[sd.imin:sd.imax], sd.DelR2[sd.imin:sd.imax])
#          sd.DiffCoeff = sd.popt[0]/6.0
#          set_D_in_DB(Directory, sd)
#          axR2. plot(sd.ttime[sd.imin:sd.imax], line(sd.ttime[sd.imin:sd.imax], sd.popt[0], sd.popt[1]), ':', linewidth=2.5, color='k')
#          #linea = ""+sd.strB2B +" "+str(sd.DiffCoeff) +"\n"
#          #ofile.write(linea)
 
  except IOError as e:
    passfig, axR2 = plt.subplots(1,1, figsize=(8,6))

  
#ofile.close()

#axR2. plot(sd.ttime, 1.5*sd.ttime, '--', color='r' )

handles, labels = axR2.get_legend_handles_labels()
leg = axR2.legend(handles, labels, bbox_to_anchor=(0.01,0.08), loc="lower left", title="$p_{S}$", prop={'size':12}, frameon=False, handlelength=0.8)
#figs[b2f].savefig("MSD_vs_t__pBF_"+strB2F +"__Hrad_"+strHrad+".pdf")
fig.savefig("Fig3_Main.pdf")




#### Inset Figure

left, bottom, width, height = [0.45, 0.48, 0.48, 0.44]
axDiff = fig.add_axes([left, bottom, width, height])
#axDiff.set_title("$ R_H = "+strHrad+" $ ")#+", $ p_{BF} = "+strB2F+" $")

axDiff.set_xscale("log")
axDiff.set_yscale("log")
axDiff.set_xlabel("$p_{S}$") #, fontsize='30')
axDiff.set_ylabel("$D$", rotation=0) #, fontsize='30')
#axDiff.set_ylim((0,0.2))

axDiff.tick_params(axis='x', labelsize=10)
axDiff.tick_params(axis='y', labelsize=10)  


axDiff.tick_params(axis='both', which='major', left='on', right='on', top='on', bottom='on')  
axDiff.tick_params(axis='both', which='minor', left='off', right='off', top='off', bottom='off')  


figScaDiff, axScaDiff = plt.subplots(1,1, figsize=(8,6))
axScaDiff.set_xlabel("$p_S t$", fontsize='30', rotation=0)
axScaDiff.set_ylabel("$ r_2 / (p_S \\times t)$", fontsize='30')

axScaDiff.set_xscale("log")
axScaDiff.set_yscale("log")

axScaDiff.set_xlabel("$p_{S}/p_{D}$", fontsize='30')
axScaDiff.set_ylabel("$D/p_{D}$", fontsize='30')

#strB2FList = ["0.000001", "0.00001", "0.0001", "0.001", "0.01" ]

for b2f in range(len(strB2FList)):
  strB2F = strB2FList[b2f]
  pB2F   = float(strB2F)
  outPrefix = "D_vs_pBB__pBF_"+strB2F +"__Hrad_"+strHrad
  outfln_PDF = outPrefix+".pdf"
  outfln_DAT = outPrefix+".DAT"

  # D vs pB2B
  #print outfln_DAT
  try:
    pB2B, DiffCoeff = np.loadtxt(outfln_DAT).T
    print pB2B, DiffCoeff
    estiloDiff = lineStyle_axDiff(strB2F)
    axDiff.plot(pB2B, DiffCoeff, **estiloDiff)
    axScaDiff.plot(pB2B/pB2F, DiffCoeff/pB2F, **estiloDiff)
  except ValueError as e:
    pass
      

handles, labels = axDiff.get_legend_handles_labels()
axDiff.legend(handles, labels, bbox_to_anchor=(0.75,0.03), loc="lower left", title="$p_{D}$", prop={'size':10}, frameon=False, numpoints=1)

fig.tight_layout()
fig.savefig("Fig3_paper.pdf")

axScaDiff.legend(handles, labels, bbox_to_anchor=(0.75,0.03), loc="lower left", title="$p_{D}$", prop={'size':10}, frameon=False, numpoints=1)

figScaDiff.tight_layout()


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



