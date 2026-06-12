


import src.jkat as jkat
import src.jkat.utils as utils
from src.jkat.ephemeris.JPLHorizons import horizons_request
import math as m
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
import datetime as dt

# jkat.add_solar_system()
# o = jkat.examples.Omuamua
# jkat.plot(o, t=o.tp, t_bounds=(-m.inf, o.tp), max_distance=6*AU)
# jkat.show()

''' 
need trajectories: 
intercept, 
rendezvous
'''
from jkat.ephemeris.JPLHorizons import horizons_request

horizons_request('Eris', t_start=jkat.ephemeris.to_time(dt.datetime(2035,1,1)), t_end=jkat.ephemeris.to_time(dt.datetime(2037,1,1)))