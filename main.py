


import src.jkat as jkat
import src.jkat.utils as utils
from src.jkat.ephemeris.JPLHorizons import horizons_request
import math as m
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
import datetime as dt



jkat.add_solar_system()
o = jkat.ephemeris.examples.Omuamua
jkat.plot(o, t=o.tp, t_bounds=(-m.inf, o.tp+jkat.YEAR), max_distance=6*jkat.AU, stilt_spacing='range', stilt_number=30)
jkat.show()

