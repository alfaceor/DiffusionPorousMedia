#!/usr/bin/python

strF2B="1.0"
strB2F="0.1"
strB2B="0.0001"
strHrad="20"
strBrad="21"

Hrad = int(strHrad)

UcellxSize=(2*(Hrad)+1)
UcellySize = UcellxSize
UcellzSize = UcellxSize

import numpy as np
from mayavi import mlab

mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0))

#dataTraj = np.loadtxt("trajectoryTest.csv", dtype=float)
filenameTraj= "SingleTrajectory/rw3D_RegularPorousDepletion_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHrad+"__Brad_"+strBrad+".trj"
dataTraj = np.loadtxt(filenameTraj, dtype=float)
dataTraj = dataTraj[2000:5260]
###   #dataTraj = dataTraj[:1260]

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

i_rnd = int( dataNew[0][1])
j_rnd = int( dataNew[0][2])
k_rnd = int( dataNew[0][3])

indx = i_rnd % UcellxSize;
if ( indx < 0 ): indx = indx + UcellxSize;
indy = j_rnd % UcellySize;
if ( indy < 0 ): indy = indy + UcellySize;
indz = k_rnd % UcellzSize;
if ( indz < 0 ): indz = indz + UcellzSize;

x0 = dataNew[0, 1] 
y0 = dataNew[0, 2] 
z0 = dataNew[0, 3] 

dataNew = dataNew -np.array([0, x0, y0, z0]) +  np.array([0, indx, indy, indz])

xp = dataNew[:, 1] 
yp = dataNew[:, 2] 
zp = dataNew[:, 3] 
#print dataNew
#mlab.points3d(x, y, z) # s, colormap="copper", scale_factor=.25, mode='cube')

print x[0], y[0], z[0]
print xp[0], yp[0], zp[0]

#mlab.plot3d(x, y, z, tube_radius=0.25, color=(0,0,1)) #, representation='points')
mlab.plot3d(xp, yp, zp, tube_radius=0.25, color=(1,0,0)) #, representation='points')
#mlab.points3d(x[-1], y[-1], z[-1] , mode='sphere', color=(1,0,0))
mlab.points3d(xp[-1], yp[-1], zp[-1] , mode='sphere', color=(0,0,1))


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
for ii in range(-1,1):
  for jj in range(0,2):
    for kk in range(0,2):
      mlab.points3d(XX + ii*Delta, YY + jj*Delta, ZZ + kk*Delta, **cubeArgs) 


#mlab.savefig("figura.eps")

#mlab.orientation_axes()

mlab.show()

