# Calculating the order of accuracy of an RK scheme from the Butcher Tableau

import numpy as np
import sys
from ButcherTableau import *

RK3 = ButcherTableau('RK3',
                     [[0,0,0,0],
                      [1,0,0,0],
                      [0.25,0.25,0,0],
                      [1/6,1/6,2/3,0]],
                      [1/6,1/6,2/3,0])

RK3i = ButcherTableau('RK3i',
                      [[0,0,0,0], [0,1,0,0],
                      [0.25,0,0.25,0],
                      [1/6,1/6,2/3,0]],
                      [1/6,1/6,2/3,0])
IMEX = doubleButcher(RK3, RK3i)
IMEX.order()
RK3i.describe()
tmp = RK3i.plotStabilityFunction(np.linspace(0,10,101)**2)

RK24i = ButcherTableau('RK24i',
                      [[0,0,0,0], [0,1,0,0],
                      [0.25,0,0.25,0],
                      [0.5,0,0,0.5]],
                      [0.5,0,0,0.5])
IMEX = doubleButcher(RK3, RK24i)
IMEX.order()
RK24i.describe()

SSP2332i = ButcherTableau('SSP2332i',
                            [[0.25,0,0],
                           [0,0.25,0],
                           [1/3,1/3,1/3]],
                           [1/3,1/3,1/3])

ARS3i = ButcherTableau('ARS3i',
                        [[0.5,0,0,0],
                        [1/6,0.5,0,0],
                        [-0.5,0.5,0.5,0],
                        [1.5,-1.5,0.5,0.5]],
                        [1.5,-1.5,0.5,0.5])


g = 1- 2**(-0.5)
SSP3332i = ButcherTableau('SSP3332i',
                        [[g,0,0],
                           [1-2*g,g,0],
                           [0.5-g,0,g]],
                           [1/6,1/6,2/3])

g = (3+np.sqrt(3))/6
ARS3233i = ButcherTableau('ARS3233i',[[g,0],[1-2*g,g]],[0.5,0.5])

# From Kennedy and Carpenter
c2 = 4/7
c3 = 1
BSDIRK = ButcherTableau('BSDIRK',
                        [[1/6,0,0],
                        [c2-1/6,1/6,0],
[(6*c3-1)*(3-24*c2+24*c2**2+2*c3)/(6*(4*c2-3)*(6*c2-1)),
 (6*c3-1)*(c2-c3)/(3*(4*c2-3)*(6*c2-1)),
 1/6]],
 [6*(2-3*c2-3*c3+6*c2*c3)/((6*c2-1)*(6*c3-1)),
  (3-4*c3)/(2*(6*c2-1)*(c2-c3)),
  (4*c2-3)/(2*(c2-c3)*(6*c3-1))])

