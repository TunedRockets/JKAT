


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

print('sending request')
ob = horizons_request('499', center='0', t_start=0, t_end=5*jkat.YEAR)
ob2 = horizons_request('499', center='10', t_start=0, t_end=5*jkat.YEAR)
print('request done')

tt =np.linspace(0,5*jkat.YEAR, 100)
rr = []
oo = []
for t in tt:
    ob.current_time = t
    ob2.current_time = t
    rr.append(ob.i)
    oo.append(ob2.i)

plt.plot(tt,rr)
plt.plot(tt,oo)
plt.show()
