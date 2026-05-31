'''
direct single burn changes in orbit
'''

import math as m
from ..utils import *
from ..kep.orbits import Orbit
from ..kep import determination as det
import numpy as np
from typing import overload

__all__ = [
    'single_burn',
    'orbit_rotation',
]



@overload
def single_burn(ob:Orbit, dVvec:np.ndarray, t:float, f:None=None)->Orbit:...
@overload
def single_burn(ob:Orbit, dVvec:np.ndarray, t:None=None, f:float=...)->Orbit:...
@overload
def single_burn(ob:Orbit, dVvec:np.ndarray, t:None=None, f:None=None)->Orbit:...


def single_burn(ob:Orbit, dVvec:np.ndarray, t:float|None=None,f:float|None=None)->Orbit:
    ''' make a single burn with the given dv in cartesian coordinates at the given time or angle, and returns new orbit
    '''

    if (not t is None) and (not f is None): raise ValueError("Overconstrained burn, both t and f given")

    if not t is None: f = ob.f(t)
    if f is None: f = ob.f(0)

    if t is None: t = ob.t(f) # since we need t later as well

    rvec, vvec = ob.vectors(f)

    vvec += dVvec

    ob2 = det.orbit_from_rv(rvec,vvec,ob.mu,t)
    return ob2

def orbit_rotation(ob:Orbit, angle:float, t:float|None=None,f:float|None=None, conservative:bool=True)->tuple[np.ndarray, Orbit]:
    '''rotate the orbit around it's r vector at a given point. a conservative rotation keeps the magnitude of the 
    velocity vector constant, keeping the orbit the same, a non-conservative projects the velocity vector in the new direction,
    minimizing expended dV,
    returns the delta v vector and the resulting orbit'''

    # get time and true anomaly:
    if (not t is None) and (not f is None): raise ValueError("Overconstrained burn, both t and f given")
    if not t is None: f = ob.f(t) # set f from t
    if f is None: f = ob.f(0) # no input given
    if t is None: t = ob.t(f) # since we need t later as well

    # get r,v vectors:
    rvec,vvec = ob.vectors(f)

    vvec2 = rodrigues_rot(vvec,rvec,angle)

    if not conservative:
        # project vvec onto vvec2 for new vvec2
        vproj = vvec.dot(vvec2)/(vvec2.dot(vvec2)) * vvec2
        vvec2 = vproj

    dvvec = vvec2 - vvec
    ob2 = det.orbit_from_rv(rvec,vvec2,ob.mu,t)
    return dvvec, ob2
