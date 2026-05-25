


from src.jkat.utils import *
from src.jkat.kep import *
from src.jkat.ephemeris import *
from src.jkat.trajectories import porkchop_plot, direct_transfer
from src.jkat.ephemeris.examples import Omuamua

import src.jkat.plotting as oplt
import matplotlib.pyplot as plt
import numpy as np
from tests.test_timings import test_uni_round_trip_hyper
from tests.test_curtis import test_curtis_4_7
from pprint import pprint


import src.jkat as jkat

jkat.add_solar_system()
o = jkat.examples.Omuamua
jkat.plot(o, t=o.tp, t_bounds=(-m.inf, o.tp), max_distance=6*AU)
jkat.show()


# st = np.linspace(0.5,2, 100)*JULIAN_YEAR
# et = np.linspace(1.25,2.5, 100)*JULIAN_YEAR

# dv = porkchop_plot(Earth, Mars, st,et)

# res = direct_transfer(Earth,Mars, bounds=(st[0],st[-1],et[0],et[-1]))


# plt.imshow(dv.T, origin="lower", extent=(st[0]/JULIAN_YEAR,st[-1]/JULIAN_YEAR,et[0]/JULIAN_YEAR,et[-1]/JULIAN_YEAR), vmax=20)
# plt.scatter(res['ts'], res['te'])
# plt.xlabel("start time [Years]")
# plt.ylabel("end time [Years]")
# plt.show()
