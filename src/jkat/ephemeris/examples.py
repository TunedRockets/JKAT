''' 
Examples of actual objects, to use as a reference or for trajectory analysis

'''
from ..kep import Orbit, orbit_from_ephemeris
from ..utils import JULIAN_CENTURY, AU, SUN_MU, pe2p
from .JPLHorizons import horizons_request
from .realtime import to_time
from datetime import datetime

from pathlib import Path
import math as m
__all__ = [
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    "Jupiter",
    'Saturn',
    'Uranus',
    'Neptune',
    'Pluto'
]

t_start = datetime(1950,1,1)
t_end = datetime(2050,1,1)
N = 200

Mercury =  horizons_request(199, t_start=t_start, t_end=t_end, steps=N)
Venus =  horizons_request(299, t_start=t_start, t_end=t_end, steps=N)
Earth =  horizons_request(399, t_start=t_start, t_end=t_end, steps=N)
Mars =  horizons_request(499, t_start=t_start, t_end=t_end, steps=N)
Jupiter =  horizons_request(599, t_start=t_start, t_end=t_end, steps=N)
Saturn =  horizons_request(699, t_start=t_start, t_end=t_end, steps=N)
Uranus =  horizons_request(799, t_start=t_start, t_end=t_end, steps=N)
Neptune =  horizons_request(899, t_start=t_start, t_end=t_end, steps=N)
Pluto =  horizons_request(999, t_start=t_start, t_end=t_end, steps=N)


# ISOs:

Omuamua = Orbit(
    p = pe2p(0.255916*AU, 1.20113),
    e = 1.20113,
    i = m.radians(122.74),
    argp = m.radians(241.811),
    raan = m.radians(24.597),
    tp = to_time(datetime(2017,9,9)),
    mu=SUN_MU
)

Borisov = Orbit(
    p = pe2p(2.00652*AU, 3.3565),
    e=3.3565,
    i=m.radians(44.053),
    raan=m.radians(308.15),
    argp=m.radians(209.12),
    tp=to_time(datetime(2019,12,8)),
    mu=SUN_MU
)

ATLAS = Orbit(
    p = pe2p(1.35645*AU, 6.14135),
    e=6.14135,
    i = m.radians(175.12),
    raan=m.radians(322.17),
    argp = m.radians(128.02),
    tp=to_time(datetime(2025,10,29)),
    mu=SUN_MU
)