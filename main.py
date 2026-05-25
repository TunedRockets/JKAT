


from src.utils import *
from src.kep import *
from src.ephemeris import *
from src.trajectories import porkchop_plot

import src.plotting as oplt
import matplotlib.pyplot as plt
import numpy as np
from tests.test_timings import test_uni_round_trip_hyper
from tests.test_curtis import test_curtis_4_7
from pprint import pprint




st = np.linspace(0.5,2, 200)*JULIAN_YEAR
et = np.linspace(1.25,2.5, 200)*JULIAN_YEAR

dv = porkchop_plot(Earth, Mars, st,et)
pprint(dv)
plt.imshow(dv.T, origin="lower", extent=(st[0]/JULIAN_YEAR,st[-1]/JULIAN_YEAR,et[0]/JULIAN_YEAR,et[-1]/JULIAN_YEAR), vmax=20)
plt.xlabel("start time [Years]")
plt.ylabel("end time [Years]")
plt.show()
