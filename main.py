


from src.utils import *
from src.kep import *
from src.ephemeris import *
import src.plotting as oplt
import matplotlib.pyplot as plt
import numpy as np
from tests.test_timings import test_uni_round_trip_hyper
from tests.test_curtis import test_curtis_4_7
from pprint import pprint



pprint(JPL_ephemeris)

oplt.add_solar_system(planets='11110000', initials=True)
oplt.show()

