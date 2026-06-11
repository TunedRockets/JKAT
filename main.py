


import src.jkat as jkat
import src.jkat.utils as utils
from src.jkat.ephemeris.JPLHorizons import horizons_request
import math as m
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt


# jkat.add_solar_system()
# o = jkat.examples.Omuamua
# jkat.plot(o, t=o.tp, t_bounds=(-m.inf, o.tp), max_distance=6*AU)
# jkat.show()

''' 
need trajectories: 
intercept, 
rendezvous
'''

ap = jkat.AU * 5.45
pe = 10*jkat.SUN_RADIUS
a, e = utils.apse2ae(ap,pe)
p = utils.a2p(a,e)
ob = jkat.Orbit(p,e,0,0,0,0,jkat.SUN_MU)
v1 = ob.v(m.pi)

ap = jkat.AU * 5.45
pe = 5*jkat.SUN_RADIUS
a, e = utils.apse2ae(ap,pe)
p = utils.a2p(a,e)
ob = jkat.Orbit(p,e,0,0,0,0,jkat.SUN_MU)
v2 = ob.v(m.pi)

print(f"difference: {v2-v1} km/s")
