''' 
Examples of actual objects, to use as a reference or for trajectory analysis

'''
from ..kep import Orbit, orbit_from_ephemeris
from ..utils import JULIAN_CENTURY, AU, SUN_MU
from .realtime import to_time

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
    filepath = Path(__file__).parent.parent.parent / Path('data/timekeeping/JPL_ephemeris.txt')
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