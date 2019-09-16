#!/usr/bin/python

strF2B="1.0"
strB2F="0.001"
strB2B="0.0001"
strHrad="10"
strBrad="11"

Hrad = int(strHrad)

import numpy as np
from mayavi import mlab

mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0))

dataGeom = np.loadtxt("rw3D_Geometry__Hrad_"+strHrad+"__Brad_"+strBrad+"__HARDWALL.csv", dtype=int)

XX = dataGeom[:,0]
YY = dataGeom[:,1]
ZZ = dataGeom[:,2]

Delta = (2*Hrad +1)
cubeArgs={"mode"         :  'cube',
          "color"        :  (210./256., 180./256., 140./256.), 
          "opacity"      :  1.0, 
          "scale_factor" :  1.0
         }
for ii in range(-2, 2):
  for jj in range(-2, 2):
    for kk in range(-2,2):
      mlab.points3d(XX + ii*Delta, YY + jj*Delta, ZZ + kk*Delta, **cubeArgs) 


#mlab.savefig("figura.eps")

mlab.orientation_axes()

mlab.show()

