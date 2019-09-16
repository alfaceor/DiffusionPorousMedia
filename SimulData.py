# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 17:29:54 2018

@author: alfaceor
"""
import numpy as np
import matplotlib.pyplot as plt
import os

class SimulData:
  def __init__(self):
    self.strF2B     = "" #strF2B
    self.strB2F     = "" #strB2F
    self.strB2B     = "" #strB2B
    self.strHrad    = "" #strHrad
    self.strBrad    = "" #strBrad
    self.strNtrials = "" #strNtrials
    self.flnPrefix  = "" #"rw3D_RegularPorousDepletion"
    self.arguments  = "" #" --pF2B "+self.strF2B+" --pB2F "+self.strB2F+" --pB2B "+self.strB2B+" --Hradius "+self.strHrad+" --Bradius "+self.strBrad
    self.qsubArgs   = ""
    self.netFile    = "" #self.flnPrefix+"_F2B_"+self.strF2B+"__B2F_"+self.strB2F+"__B2B_"+self.strB2B+"__Hrad_"+self.strHrad+"__Brad_"+self.strBrad+".net"
    self.dataFile   = "" #flnPrefix+"_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHrad+"__Brad_"+strBrad+"__nt_"+strNtrials+".dat"
    self.path       = ""
    
    # Converted values
    self.pF2B       = -1 #float(strF2B)
    self.pB2F       = -1 #float(strB2F)
    self.pB2B       = -1 #float(strB2B)
    self.Hrad       = -1 #float(strHrad)
    self.Brad       = -1 #float(strBrad)
    self.Ntrials    = -1 #int  (strNtrials)
    
    
    #calculated values
    #data             = [] #np.loadtxt(self.dataFile)
    self.ttime       = [] #data[:,0]
    self.DelR2       = [] #data[:,7]
    self.ttimeLambda = []
    self.EffecDiff   = []
    
    # FIT DATA
    self.fitparams = None
    self.tmin      = -1 #0
    self.tmax      = -1 #len(self.ttime)
    self.imin      = -1 #0
    self.imax      = -1 #len(self.ttime)
    self.DiffCoeff = -1 #0.0
    self.D_err_min = -1
    self.D_err_max = -1
    self.popt      = []
    self.pcov      = []

    # Fabio Method
    self.lambFlag  = False
    self.strLamb1  = ""
    self.strLamb2  = ""
    self.strLamb3  = ""
    self.lamb1     = -1
    self.lamb2     = -1
    self.lamb3     = -1
    self.DiffLamb1 = -1
    self.DiffLamb2 = -1
    self.DiffLamb3 = -1
    
    
  def metadata(self, flnPrefix, strF2B, strB2F, strB2B, strHrad, strBrad, strNtrials, loadData=False):
    self.flnPrefix  = flnPrefix
    self.strF2B     = strF2B
    self.strB2F     = strB2F
    self.strB2B     = strB2B
    self.strHrad    = strHrad
    self.strBrad    = strBrad
    self.strNtrials = strNtrials
    
    self.convertStrParams()
    self.defineFiles()
    if (loadData == True):
      self.getData()
    
    
  def metadataFromCSVentry(self, flnPrefix, csv_entry):
    self.flnPrefix  = flnPrefix
    self.strF2B     = csv_entry["strF2B"    ]
    self.strB2F     = csv_entry["strB2F"    ]
    self.strB2B     = csv_entry["strB2B"    ]
    self.strHrad    = csv_entry["strHrad"   ]
    self.strBrad    = csv_entry["strBrad"   ]
    self.strNtrials = csv_entry["strNtrials"]
    
    self.convertStrParams()
    self.defineFiles()

  def convertStrParams(self):
    self.pF2B       = float(self.strF2B)
    self.pB2F       = float(self.strB2F)
    self.pB2B       = float(self.strB2B)
    self.Hrad       = int  (self.strHrad)
    self.Brad       = int  (self.strBrad)
    self.Ntrials    = int  (self.strNtrials)
    
  def defineFiles(self):
    self.arguments  = " --pF2B "+self.strF2B+" --pB2F "+self.strB2F+" --pB2B "+self.strB2B+" --Hradius "+self.strHrad+" --Bradius "+self.strBrad
    self.netFile    = self.flnPrefix+"_F2B_"+self.strF2B+"__B2F_"+self.strB2F+"__B2B_"+self.strB2B+"__Hrad_"+self.strHrad+"__Brad_"+self.strBrad+".net"
    self.dataFile   = self.flnPrefix+"_F2B_"+self.strF2B+"__B2F_"+self.strB2F+"__B2B_"+self.strB2B+"__Hrad_"+self.strHrad+"__Brad_"+self.strBrad+"__nt_"+self.strNtrials+".dat"
    self.qsubArgs  =" -v F2B="+self.strF2B+",B2F="+self.strB2F+",B2B="+self.strB2B+",Hrad="+self.strHrad+",Brad="+self.strBrad

  
  def getData(self):
    #self.defineFiles()
    #calculated values
    data       = np.loadtxt(self.dataFile)
    self.ttime = data[:,0]
    self.DelR2 = data[:,7]
    
    # FIT DATA
    self.tmin  = 0
    self.tmax  = len(self.ttime)
    self.DiffCoeff  = 0.0
    self.popt  = []
    self.pcov  = []
    
  def isDone(self):
    self.defineFiles()
    if os.path.isfile(self.dataFile) :
      return True  #program+" "+self.arguments # "JOB TO RUN
    else:
      if os.path.isfile(self.netFile) :
        return  True #"echo JOB is RUNNING OR IMCOMPLETE!!" 
      else:
        return False # "echo JOB is DONE"
  
  def isDatafile(self):
    self.defineFiles()
    if os.path.isfile(self.dataFile) :
      return True  #program+" "+self.arguments # "JOB TO RUN
    else:
      return False # "echo JOB is DONE"  
  
  def getCMD(self, program):
    self.defineFiles()
    if os.path.isfile(self.netFile) :
      if os.path.isfile(self.dataFile) :
        return  "echo JOB is DONE"
      else:
        return "echo JOB is probabibly RUNNING OR IMCOMPLETE!!"
    else:
      return program+" "+self.arguments # "echo JOB TO SUBMIT!!!"

  def calcDiffCoeff(self, fitparams):
    self.getData()
    self.fitparams = fitparams
    tmin, tmax     = fitparams["fitrange"]
    self.imin, self.imax     = getLimits(self.ttime, tmin, tmax)
    self.popt, self.pcov = curve_fit(line, self.ttime[self.imin:self.imax], self.DelR2[self.imin:self.imax])
    self.DiffCoeff = self.popt[0]/6.0
    self.D_err_min = 0.5*self.pcov[0][0]/6.0#D_aux.mean() - D_aux.min()
    self.D_err_max = 0.5*self.pcov[0][0]/6.0#D_aux.mean() - D_aux.min()
    

  def calcDwithLamb(self, fitparams, colorcito, titulo):
    self.calcDiffCoeff(fitparams)
    #fitparams
    tmin, tmax     =  fitparams["fitrange"]
    imin, imax     =  getLimits(self.ttime, tmin, tmax)
    self.lamb1     =  fitparams["lamb1"]
    self.lamb2     =  fitparams["lamb2"]
    self.lamb3     =  fitparams["lamb3"]
    self.strLamb1  =  str(self.lamb1)
    self.strLamb2  =  str(self.lamb2)
    self.strLamb3  =  str(self.lamb3)
    EffecDiff      =  self.DelR2[imin:imax]/(6.*self.ttime[imin:imax])
    # lambda 1
    ttimeLambda1 = 1./self.ttime[imin:imax]**self.lamb1
    ttimeLambda2 = 1./self.ttime[imin:imax]**self.lamb2
    ttimeLambda3 = 1./self.ttime[imin:imax]**self.lamb3
    
    popt, pcov = fitLine(ttimeLambda1,EffecDiff)
    slope1 = popt[0]
    self.DiffLamb1 = popt[1]
    popt, pcov = fitLine(ttimeLambda2,EffecDiff)
    slope2 = popt[0]
    self.DiffLamb2 = popt[1]
    popt, pcov = fitLine(ttimeLambda3,EffecDiff)
    slope3 = popt[0]
    self.DiffLamb3 = popt[1]

    D_aux = np.array([self.DiffCoeff, self.DiffLamb1, self.DiffLamb2, self.DiffLamb3])
    self.D_err_min = D_aux.mean() - D_aux.min()
    self.D_err_max = D_aux.max()  - D_aux.mean()
    
    # plot the lambdas
    fig, ax = plt.subplots(2,2, figsize=(8,6))
    fig.suptitle(titulo)
    axR2    = ax[0][0]
    axDiff1 = ax[0][1]
    axDiff2 = ax[1][0]
    axDiff3 = ax[1][1]
    
    axR2   .tick_params(axis='both', labelsize=10)
    axDiff1.tick_params(axis='both', labelsize=10)
    axDiff2.tick_params(axis='both', labelsize=10)
    axDiff3.tick_params(axis='both', labelsize=10)
    
    
    axR2.set_xscale("log")
    axR2.set_yscale("log")
    axR2.tick_params(axis='x', labelsize=10)
    axR2.tick_params(axis='y', labelsize=10)  
    
    axR2.plot(self.ttime, self.DelR2, linewidth=0,  color='b', marker="o", markersize=1)
    axR2.plot(self.ttime[imin:imax], line(self.ttime[imin:imax], self.popt[0], self.popt[1]), '--', color='k' )
    
    linea, = axR2. plot(self.ttime, self.ttime , '--', linewidth=2, color='r', alpha=0.5 )
    dr = DraggableLogLogLine(linea)
    dr.connect()
    
    
    Y1 = D_aux.min()   *np.ones(len(ttimeLambda1))
    Y2 = D_aux.max()   *np.ones(len(ttimeLambda1))
    XX = np.linspace(0, ttimeLambda1.max(), num=len(ttimeLambda1))
    axDiff1.fill_between(XX, Y1, Y2, facecolor="white", alpha=0.5, hatch="x", edgecolor='g', linestyle='dashed')
    XX = np.linspace(0, ttimeLambda2.max(), num=len(ttimeLambda2))
    axDiff2.fill_between(XX, Y1, Y2, facecolor="white", alpha=0.5, hatch="x", edgecolor='g', linestyle='dashed')
    XX = np.linspace(0, ttimeLambda3.max(), num=len(ttimeLambda3))
    axDiff3.fill_between(XX, Y1, Y2, facecolor="white", alpha=0.5, hatch="x", edgecolor='g', linestyle='dashed')
    
    axDiff1.set_xlabel("$ t^{-"+self.strLamb1+"}$")
    axDiff1.ticklabel_format( useMathText=True, style='sci', axis='both', scilimits=(0,0) )
    axDiff1.plot(ttimeLambda1, EffecDiff, color=colorcito, marker='o' )
    XX = np.linspace(0,ttimeLambda1.max(), num=100 )
    axDiff1.plot(XX, line(XX, slope1, self.DiffLamb1), '--', color='k' )
    axDiff1.axhline(y=self.DiffCoeff, color='k')
    
    
    axDiff2.set_xlabel("$ t^{-"+self.strLamb2+"}$")
    axDiff2.ticklabel_format( useMathText=True, style='sci', axis='both', scilimits=(0,0) )
    axDiff2.plot(ttimeLambda2, EffecDiff, color=colorcito, marker='o' )
    XX = np.linspace(0,ttimeLambda2.max(), num=100 )
    axDiff2.plot(XX, line(XX, slope2, self.DiffLamb2), '--', color='k' )
    axDiff2.axhline(y=self.DiffCoeff, color='k')
    
    axDiff3.set_xlabel("$ t^{-"+self.strLamb3+"}$")
    axDiff3.ticklabel_format( useMathText=True, style='sci', axis='both', scilimits=(0,0) )
    axDiff3.plot(ttimeLambda3, EffecDiff, color=colorcito, marker='o' )
    XX = np.linspace(0,ttimeLambda3.max(), num=100 )
    axDiff3.plot(XX, line(XX, slope3, self.DiffLamb3), '--', color='k' )
    axDiff3.axhline(y=self.DiffCoeff, color='k')
    
    
    
    axDiff1.xaxis.get_offset_text().set_fontsize(10)
    axDiff2.xaxis.get_offset_text().set_fontsize(10)
    axDiff3.xaxis.get_offset_text().set_fontsize(10)
    
    axDiff1.yaxis.get_offset_text().set_fontsize(10)
    axDiff2.yaxis.get_offset_text().set_fontsize(10)
    axDiff3.yaxis.get_offset_text().set_fontsize(10)
    
    fig.tight_layout()
    fig.savefig("FIT_F2B_"+self.strF2B+"__B2F_"+self.strB2F+"__B2B_"+self.strB2B+"__Hrad_"+self.strHrad+"__Brad_"+self.strBrad+"__nt_"+self.strNtrials+".png")
    
    return fig, dr
  
  
  def getCSVstrLine(self):
    return  \
    "\"" +self.strF2B    +"\","+ \
    "\"" +self.strB2F    +"\","+ \
    "\"" +self.strB2B    +"\","+ \
    "\"" +self.strHrad   +"\","+ \
    "\"" +self.strBrad   +"\","+ \
    "\"" +self.strNtrials+"\""  


from scipy.optimize import curve_fit
from scipy.stats import norm

def line(x, m, b): 
  return m*x+b

def fitLine(xdata, ydata):
  popt, pcov = curve_fit(line, xdata, ydata)
  return popt, pcov

def powLaw(x, m, b): 
  return np.exp(b)*x**m

def fitPowLaw(xdata,ydata):
  popt, pcov = curve_fit(line, np.log(xdata), np.log(ydata))
  return popt, pcov

def fitPowLawFromTo(xdata,ydata, xmin=None, xmax=None, returnIndx=False):
  if (xmin == None) and (xmax == None):
    popt, pcov = curve_fit(line, np.log(xdata), np.log(ydata))
    return popt, pcov
  else:
    if xmin == None:
      xmin = xdata[0]
    if xmax == None:
      xmax = xdata[-1]
    i_min, i_max = getLimits(xdata, xmin, xmax)
    popt, pcov = curve_fit(line, np.log(xdata[i_min:i_max]), np.log(ydata[i_min:i_max]))
    if returnIndx:
      return popt, pcov, i_min, i_max
    else:
      return popt, pcov
    


def gaussian(x, loc, scale):
  return norm.pdf(x, loc, scale)

def fitGaussian(xdata, ydata):
  popt, pcov = curve_fit(gaussian, xdata, ydata)
  return popt, pcov

def getLimits(arr, amin, amax):
  # Array must be sorted from min to max
  i_min=0
  i_max=len(arr)-1
  #if arr[i_min] > amin or arr[i_max] < amax:
  #  print "Out of bounds!!"
  #  exit()
  ii    = i_min
  while (arr[ii] <= amin ):
    ii = ii+1
    i_min = ii

  ii    = i_max
  while (arr[ii] >= amax ):
    ii = ii-1
    i_max = ii

  return i_min, i_max


def latex_float(f):
  float_str = "{0:.2g}".format(f)
  if "e" in float_str:
    base, exponent = float_str.split("e")
    #print float_str, base, exponent
    if base == "1" :
      #print base
      return r"10^{{{0}}}".format(int(exponent))
    else:
      return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
  else:
    return float_str


def latex_float02(f):
  float_str = "{0:.2e}".format(f)
  base, exponent = float_str.split("e")
  if (float(base) == 1):
    if (float(exponent) == 0):
      return r"1"
    return r"10^{{{0}}}".format(int(exponent))
  return r"10^X"


class DraggableLogLogLine:
  def __init__(self, linea):
    self.linea = linea
    self.press = None

  def connect(self):
    'connect to all the events we need'
    self.cidpress = self.linea.figure.canvas.mpl_connect(
      'button_press_event', self.on_press)
    self.cidrelease = self.linea.figure.canvas.mpl_connect(
      'button_release_event', self.on_release)
    self.cidmotion = self.linea.figure.canvas.mpl_connect(
      'motion_notify_event', self.on_motion)

  def on_press(self, event):
    'on button press we will see if the mouse is over us and store some data'
    if event.inaxes != self.linea.axes: return

    contains, attrd = self.linea.contains(event)
    if not contains: return
    #print('event contains', self.linea.xy)
    #x0, y0 = self.linea.xy
    x0 = self.linea.get_xdata()
    y0 = self.linea.get_ydata()
    self.press = x0, y0, event.xdata, event.ydata

  def on_motion(self, event):
    'on motion we will move the linea if the mouse is over us'
    if self.press is None: return
    if event.inaxes != self.linea.axes: return
    x0, y0, xpress, ypress = self.press
    dx = event.xdata / xpress
    dy = event.ydata / ypress
    #print('x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f' %
    #   (x0, xpress, event.xdata, dx, x0+dx))
    self.linea.set_xdata(x0*dx)
    self.linea.set_ydata(y0*dy)

    self.linea.figure.canvas.draw()


  def on_release(self, event):
    'on release we reset the press data'
    self.press = None
    self.linea.figure.canvas.draw()

  def disconnect(self):
    'disconnect all the stored connection ids'
    self.linea.figure.canvas.mpl_disconnect(self.cidpress)
    self.linea.figure.canvas.mpl_disconnect(self.cidrelease)
    self.linea.figure.canvas.mpl_disconnect(self.cidmotion)


