


import jkat
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from tqdm import tqdm
parameter = 29508 # [km]
eccentricity = 0.0002931 # [-]
inclination = np.radians(56.9925) # [deg]
right_ascention_ascending_node = np.radians(341.627) # [deg]
argument_perigee = np.radians(54.4577) # [deg]
time_periapsis = dt.datetime(2026,6,22, 15,10,13)
time_periapsis = jkat.to_time(time_periapsis) # we convert the date into seconds

ob = jkat.Orbit(
        p = parameter,
        e = eccentricity,
        i = inclination,
        raan = right_ascention_ascending_node,
        argp = argument_perigee,
        tp = time_periapsis,
        mu = jkat.EARTH_MU
)

