



import src.jkat as jkat
import src.jkat.utils as utils
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


# ISO = jkat.examples.ATLAS

ISO = jkat.Orbit(
    utils.pe2p(0.7*jkat.AU, 1.4),
    1.4,
    0.5,
    0.4,
    1.3,
    jkat.JULIAN_YEAR,
    jkat.SUN_MU
)

Earth = jkat.Earth
tp = ISO.tp
dt = jkat.JULIAN_YEAR

res = jkat.direct_transfer(Earth,ISO, (tp-dt, tp+dt, tp-dt, tp+dt), dv1_w = 1, dv2_w = 1)

trans = jkat.orbit_from_lambert_transfer(Earth,ISO,res['ts'],res['te'])
pprint(res)

print(trans)
jkat.plot(trans,t_bounds=(res['ts'], res['te']), color="purple")
jkat.plot(ISO, t=res['te'], t_bounds=(-m.inf, res['te']+dt), max_distance=8*jkat.AU, color="deeppink")
jkat.add_solar_system(symbols=True, t=res['te'])
jkat.show()
# timerange = np.linspace(tp-dt, tp+dt)
# dv = jkat.trajectories.porkchop_plot(Earth, ISO, timerange, timerange, prograde=None, rendezvous=False)

# plt.imshow(dv.T, origin="lower", extent=(timerange[0],timerange[-1],timerange[0],timerange[-1]))
# plt.scatter(res['ts'], res['te'])
# plt.show()


# st = np.linspace(0.5,2)*JULIAN_YEAR
# et = np.linspace(1.25,2.5)*JULIAN_YEAR

# dv = porkchop_plot(Earth, Mars, st,et)
# # plt.ion()
# plt.imshow(dv.T, origin="lower", extent=(st[0],st[-1],et[0],et[-1]), vmax=20)
# try: 
#     res = direct_transfer(Earth,Mars, bounds=(st[0],st[-1],et[0],et[-1]))
#     plt.scatter(res['ts'], res['te'])
#     print(res)
# except: print("BAD BAD BAD")

# plt.xlabel("start time [seconds]")
# plt.ylabel("end time [seconds]")
# plt.show()
