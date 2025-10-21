"""
Shows how to plot the CCtopo grid as contours

"""

import matplotlib.pyplot as plt
import numpy as np

from clawpack.geoclaw.topotools import Topography

# Read topography for contour lines:
CCtopo = Topography()
CCtopo.read('../DataFiles/CCtopo_1_3sec.tt3',3)

plt.contour(CCtopo.X, CCtopo.Y, CCtopo.Z, np.arange(0,21,2), colors='k')
plt.savefig('CCtopo_1_3sec.png', dpi=300)
plt.show()

