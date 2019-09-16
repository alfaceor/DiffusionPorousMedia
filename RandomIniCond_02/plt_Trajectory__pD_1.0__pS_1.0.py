#!/usr/bin/python

strF2B="1.0"
strB2F="1.0"
strB2B="1.0"
strHrad="10"
strBrad="11"

Hrad = int(strHrad)

import numpy as np
from mayavi import mlab

mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0))

#dataTraj = np.loadtxt("trajectoryTest.csv", dtype=float)
filenameTraj= "SingleTrajectory/rw3D_RegularPorousDepletion_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHrad+"__Brad_"+strBrad+".trj"
dataTraj = np.loadtxt(filenameTraj, dtype=float)
dataTraj = dataTraj[:1500]


#print x
#s = 2 + np.sin(t)

## TO AVOID PROBLEMS WITH THE PLOT REMOVE REPETIVE SEQUENTIAL VALUES
dataNew  = []
dataPrev = [None, None, None]
count =0
for i in range(len(dataTraj)):
  dataCurr = [ dataTraj[i, 1], dataTraj[i, 2], dataTraj[i, 3] ] 
  if (dataPrev != dataCurr):
    dataNew.append(dataTraj[i])
    dataPrev = dataCurr
  else:
    #print dataPrev, dataCurr
    count=count+1
    #print "-------", i

dataNew = np.array(dataNew)
#print dataNew
x = dataNew[:, 1] 
y = dataNew[:, 2] 
z = dataNew[:, 3] 

#print dataNew
#mlab.points3d(x, y, z) # s, colormap="copper", scale_factor=.25, mode='cube')
mlab.plot3d(x, y, z, tube_radius=0.25, color=(1,0,0)) #, representation='points')


############## GEOMETRY ################# 

dataGeom = np.loadtxt("rw3D_Geometry__Hrad_"+strHrad+"__Brad_"+strBrad+"__HARDWALL.csv", dtype=int)

XX = dataGeom[:,0]
YY = dataGeom[:,1]
ZZ = dataGeom[:,2]

Delta = (2*Hrad +1)
cubeArgs={"mode"         :  'cube',
          "color"        :  (210./256., 180./256., 140./256.), 
          "opacity"      :  0.2, 
          "scale_factor" :  1.0
         }
for ii in range(-3,0):
  for jj in range(-1,2):
    for kk in range(-2,1):
      mlab.points3d(XX + ii*Delta, YY + jj*Delta, ZZ + kk*Delta, **cubeArgs) 


#mlab.savefig("figura.eps")

mlab.orientation_axes()

mlab.show()

