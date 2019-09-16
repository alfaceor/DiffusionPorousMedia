#!/usr/bin/python
# Gaussian
def normalDistri(x, D, t):
#  return (1.0/(np.sqrt(2.0*np.pi*KT)))*np.exp(-0.5*(x**2.0)/KT)
  d    =  1.0;
  Norm =  (np.pi*4.0*D*t)**(-d/2.0)
  return Norm*np.exp((-x**2.0)/(4.0*D*t))

################## SQL FUNCTIONS ##################
import sqlite3
def get_D_from_DB(Dir, sd):
  db = sqlite3.connect("SimulData.db")
  cursor = db.cursor()
  sqldata = cursor.execute('''SELECT D FROM Simulation WHERE Directory = ? AND strF2B = ? AND strB2F = ? AND strB2B = ? AND strHrad = ? AND strBrad =?''', (Dir, sd.strF2B, sd.strB2F, sd.strB2B, sd.strHrad, sd.strBrad))
  row = sqldata.fetchone()
#  flag = False
#  if (row != None ):
#    FitMin, FitMax = row
#    #print row
#    sd.imin, sd.imax = getLimits(sd.ttime, FitMin, FitMax)
#    flag = True
#  db.close()
  return row[0]

from SimulData import *
import numpy as np 

################## INPUT PARAMETERS ##################
strHistoTimes =  ["1", "10", "100", "1000", "10000", "100000"] #, "100000000"];
flnPrefix = "RandomIniCond_02/rw3D_RegularPorousDepletion"  #"TotalFree/rw3D_RegularPorousDepletion"
strF2B = "1.0"
strB2F = "0.001" #"0.001"
strB2B = "0.00001"
strHrad = "20"
strBrad = "21"
strNtrials = "100000"


################## FIGURES CONFIGURATION ##################
import matplotlib.pyplot as plt
###### BEGIN PLOT DECORATION VARIABLES
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES
from scipy.stats import norm


lblF2B  = latex_float(float(strF2B))
lblB2F  = latex_float(float(strB2F))
lblB2B  = latex_float(float(strB2B))
titulo       = "$p_{FB} = "+lblF2B+", p_{BF} = "+lblB2F+", p_{BB} = "+lblB2B+", R_H = "+strHrad+"$"
#titulo       = "$p_{FB} = "+lblF2B+", p_{BB} = "+lblB2B+", R_H = "+strHrad+"$"
# -------- FIG HIS -------- #
figHis, axHis = plt.subplots(figsize=(8, 6))
axHis.set_title(titulo)
axHis.set_yscale("log")

axHis.set_xlim(-100, 150)
#axHis.set_ylim(1e-6, 1e0)

axHis.set_xlabel("$ x $"   , size=30)
axHis.set_ylabel("$ P(x) $", size=30)


## -------- FIG INSET -------- #
left, bottom, width, height = [0.64, 0.6, 0.3, 0.3]
axR2 = figHis.add_axes([left, bottom, width, height])
axR2.patch.set_alpha(0.0)
axR2.set_xlabel("$t$", fontsize='18', rotation=0)
axR2.set_ylabel("$ \\langle \\Delta r^2 \\rangle $", fontsize='18')

axR2.set_xscale("log")
axR2.set_yscale("log")
axR2.xaxis.get_offset_text().set_fontsize(10)
axR2.yaxis.get_offset_text().set_fontsize(10)

#axR2.set_xlim(1e0,1e8)
#axR2.set_ylim(1e0,9e6)

axR2.tick_params(axis='x', labelsize=10)
axR2.tick_params(axis='y', labelsize=10)  
axR2.tick_params(axis='x',which='minor',bottom='off')
axR2.tick_params(axis='y',which='minor',left='off')


# -------- FIG HIS SCA -------- #
figHisSca, axHisSca = plt.subplots(figsize=(8,6))
axHisSca.set_title(titulo)
axHisSca.set_yscale("log")

axHisSca.set_xlim(-3, 4)
axHisSca.set_ylim(1e-3, 1e1)

axHisSca.set_xlabel("$ x / \\sqrt{t} $"  , size=30)
axHisSca.set_ylabel("$ P(x)  \\sqrt{t} $", size=30)


colvals = np.abs(np.log10(np.array(map(float,strHistoTimes))))
colvals = colvals/colvals.max()
print colvals


def lineStySimul(strHistoTime, strB2F):
  tipoLinea  = "-"
  anchoLinea = 2 
  histoTime  = float(strHistoTime)
  pB2F = float(strB2F)
  lblB2F  = latex_float(pB2F)
  etiqueta   = "$"+latex_float(histoTime)+ "$" # ", "+lblB2F+"$"
  expo       = int(np.abs(np.log10(histoTime)))
  colval     = expo/8.0
  colorcito  = plt.cm.Blues_r(colval)
  #colorcito  = 'k'
  alfa = 0.5
  print expo
  if expo == 0:
    marcador = "o"
  elif expo == 1:
    marcador = "x"
  elif expo == 2:
    marcador = "s"
  else:
    marcador = (int(expo), 0, 0)
    
#  colval     = np.abs(np.log10(histoTime)) / np.abs(np.log10(1000)) #histoTime/1000.
#  if strB2F == "1.0":
#    colorcito  = plt.cm.Reds(colval)
#  if strB2F == "0.001":
#    colorcito  = plt.cm.Blues(colval)
  
  estilo = {'linestyle'   : tipoLinea, # (0, (5, 10)), 
        'marker'          : marcador,
        'markersize'      : 8,
        'markeredgewidth' : 2,
        'markeredgecolor' : colorcito,
        'color'           : colorcito,
        'linewidth'       : anchoLinea,
        'label'           : etiqueta,
        'alpha'           : alfa} 
  return estilo 


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

for i in range(len(strHistoTimes)):
  strHistoTime = strHistoTimes[i]
  histoTime    = float(strHistoTimes[i])
  estilo = lineStySimul(strHistoTime, strB2F)
  
  dataFilePrefix = flnPrefix+"_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHrad+"__Brad_"+strBrad+"__ts_"+strHistoTime
  dataFilex       = dataFilePrefix +".hisx"
  dataFiley       = dataFilePrefix +".hisy"
  dataFilez       = dataFilePrefix +".hisz"
  datax = np.loadtxt(dataFilex)
  datay = np.loadtxt(dataFiley)
  dataz = np.loadtxt(dataFilez)
  
  freq = (datax[:,2]+datay[:,2]+dataz[:,2])/3.0
  posx = 0.5*(datax[:,0]+datax[:,1])
  axHis   .plot(posx, freq, **estilo ) # "-o", color=colorcito, label="$"+strHistoTime+"$")
  
  posxSca = posx/np.sqrt(histoTime)
  freqSca = freq*np.sqrt(histoTime)
  axHisSca.plot(posxSca, freqSca, **estilo ) #"-o", color=colorcito, label="$"+strHistoTime+"$")
  
########### INSET PLOT ############3
sd = SimulData()
#print dirList[idir]
flnPrefix02 = "./"+flnPrefix 
sd.metadata(flnPrefix02, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials )
sd.getData()
estilo = lineStyDAT(sd)
axR2.plot(sd.ttime, sd.DelR2, linestyle='-', color='k') #**estilo)  
axR2.plot(sd.ttime, sd.ttime, linestyle='--', color='r') #estilo['color'])  
DifCoeff = get_D_from_DB("RandomIniCond_02/", sd)

  
#popt, pcov = fitGaussian(posxSca, freqSca)
#print popt
#XX = np.linspace(-0.25,0.25, num=100)
#YY = gaussian(XX, popt[0], popt[1])


#axHisSca    .plot(XX, YY, '--', color='r', linewidth=2, alpha=0.5)
#axHisScaZoom.plot(XX, YY, '--', color='r', linewidth=2, alpha=0.5)

XX = np.linspace(-4,4, num=500)
YY = norm.pdf(XX)
print DifCoeff
YY = normalDistri(XX, DifCoeff, 1)
axHisSca    .plot(XX, YY, ':', color='r', linewidth=2) #, alpha=0.5)
#axHisScaZoom.plot(XX, YY, '--', color='b', linewidth=2, alpha=0.5)
handles, labels = axHisSca.get_legend_handles_labels()
legSca = axHisSca.legend(handles, labels, bbox_to_anchor=(0.85,0.05), loc="lower left", title="$t$", prop={'size':11}, frameon=False, handlelength=2, numpoints=1)

handles, labels = axHis.get_legend_handles_labels()
leg = axHis.legend(handles, labels, bbox_to_anchor=(0.85,0.05), loc="lower left", title="$t$", prop={'size':11}, frameon=False, handlelength=2, numpoints=1)

if leg :
  leg.draggable()

figHis       .tight_layout()
figHisSca    .tight_layout()
#figHisScaZoom.tight_layout()

figHis       .savefig("Histogram_time"+".pdf")
figHisSca    .savefig("Histogram_time"+"_Sca"+".pdf")
#figHisScaZoom.savefig("Histogram_time"+"_ScaZoom"+".pdf")

plt.show()



