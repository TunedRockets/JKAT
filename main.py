


import jkat
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt



tle = 'ISS (ZARYA)\n' \
'1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927\n' \
'2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'

ob = jkat.orbit_from_tle(tle)

jkat.plot(ob)
jkat.plotting.center()
jkat.show()
jkat.add_solar_system()
jkat.show()
