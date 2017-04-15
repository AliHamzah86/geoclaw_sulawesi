"""
Shows how to plot the CCtopo grid as contours

"""

from pylab import *

from clawpack.geoclaw.topotools import Topography

# Read topography for contour lines:
CCtopo = Topography()
CCtopo.read('../DataFiles/CCtopo_1_3sec.tt3',3)

contour(CCtopo.X, CCtopo.Y, CCtopo.Z, arange(0,21,2), colors='k')

