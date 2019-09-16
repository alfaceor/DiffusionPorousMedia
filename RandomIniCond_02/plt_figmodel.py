#!/usr/bin/python

import numpy as np
from mayavi import mlab

mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0))



############## GEOMETRY ################# 

dx_1st = np.array( [  1, -1,  0,  0,  0,  0 ] )
dy_1st = np.array( [  0,  0,  1, -1,  0,  0 ] )
dz_1st = np.array( [  0,  0,  0,  0,  1, -1 ] )

dx_2nd = np.array( [  1, -1, -1,  1,   0,  0,  0,  0,   1, -1, -1,  1    ])
dy_2nd = np.array( [  1,  1, -1, -1,   1, -1, -1,  1,   0,  0,  0,  0    ])
dz_2nd = np.array( [  0,  0,  0,  0,   1,  1, -1, -1,   1,  1, -1, -1    ])


cubeArgs={"mode"         :  'cube',
          #"color"        :  (210./256., 180./256., 140./256.), 
          "opacity"      :  0.5, 
          "scale_factor" :  1.0
          #"representation" : 'wireframe'
          }
mlab.points3d(0, 0, 0 , mode='sphere', color=(0, 0, 1)) 
mlab.points3d(dx_1st, dy_1st, dz_1st , color=(1, 1, 0),  opacity=1.0, scale_factor=1.0, mode='cube', line_width=4.0) 
#mlab.points3d(dx_2nd, dy_2nd, dz_2nd , color=(1, 1, 0),  opacity=1.0, scale_factor=1.0, mode='cube') 

"""
XX = dataGeom[:,0]
YY = dataGeom[:,1]
ZZ = dataGeom[:,2]

Delta = (2*Hrad +1)
         }
for ii in range(-1, 1):
  for jj in range(-2, 0):
    for kk in range(0,2):
      mlab.points3d(XX + ii*Delta, YY + jj*Delta, ZZ + kk*Delta, **cubeArgs) 


#mlab.savefig("figura.eps")

mlab.orientation_axes()


"""
mlab.show()
