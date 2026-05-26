''' 
Examples of actual objects, to use as a reference or for trajectory analysis

'''
from ..kep import Orbit, orbit_from_ephemeris
from ..utils import JULIAN_CENTURY, AU, SUN_MU, pe2p
from .realtime import to_time
from datetime import date

from pathlib import Path
import math as m
__all__ = [
    'JPL_ephemeris',
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    "Jupiter",
    'Saturn',
    'Uranus',
    'Neptune'
]



def _get_jpl():
    '''quick function for getting rough ephemerides''' 
    filepath = Path(__file__).parent/ Path('JPL_ephemeris.txt') # TODO: move to dedicated data folder?
    with open(filepath,'r') as file:
        lines = file.readlines()
        headings = lines[2].split()
        data = lines[5:21]
    
    for main, derivs in zip(data[::2], data[1::2]):
        main = main.split()
        derivs = derivs.split()
        name = main.pop(0)
        dic = {}
        for i, var in enumerate(headings):
            dic[var] = float(main[i])
            if i == 0: dic[var] *= AU
            elif i > 1: dic[var] *= m.pi/180

            dic['d' + var] = float(derivs[i])
            if i == 0: dic['d' + var] *= AU
            elif i > 1: dic['d' + var] *= m.pi/180
            dic['d' + var] /= JULIAN_CENTURY
        JPL_ephemeris[name] = dic
        

JPL_ephemeris = {}
'''dictionary of JPL ephemeris values valid for 1800AD to 2050AD'''
_get_jpl()


def from_JPL(name,t): 
    d = JPL_ephemeris[name]
    return orbit_from_ephemeris(
        d['a'] + d['da']*t,
        d['e'] + d['de']*t,
        d['i'] + d['di']*t,
        d['L'] + d['dL']*t,
        d['longp'] + d['dlongp']*t,
        d['raan'] + d['draan']*t,
        SUN_MU
    )

Mercury =  from_JPL('Mercury',0)
Venus =  from_JPL('Venus',0)
Earth =  from_JPL('Earth',0)
Mars =  from_JPL('Mars',0)
Jupiter =  from_JPL('Jupiter',0)
Saturn =  from_JPL('Saturn',0)
Uranus =  from_JPL('Uranus',0)
Neptune =  from_JPL('Neptune',0)


# ISOs:

Omuamua = Orbit(
    p = pe2p(0.255916*AU, 1.20113),
    e = 1.20113,
    i = m.radians(122.74),
    argp = m.radians(241.811),
    raan = m.radians(24.597),
    tp = to_time(date(2017,9,9)),
    mu=SUN_MU
)

Borisov = Orbit(
    p = pe2p(2.00652*AU, 3.3565),
    e=3.3565,
    i=m.radians(44.053),
    raan=m.radians(308.15),
    argp=m.radians(209.12),
    tp=to_time(date(2019,12,8)),
    mu=SUN_MU
)

ATLAS = Orbit(
    p = pe2p(1.35645*AU, 6.14135),
    e=6.14135,
    i = m.radians(175.12),
    raan=m.radians(322.17),
    argp = m.radians(128.02),
    tp=to_time(date(2025,10,29)),
    mu=SUN_MU
)