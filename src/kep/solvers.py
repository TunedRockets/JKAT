''' 
Solvers for some simple astrodynamic problems,
such as lamert or gauss' problem
'''

import numpy as np
import math as m
from ..utils import stumpff_c,stumpff_s, root_finder_bisection

__all__ = [
    'lambert',
]

def lambert(r1vec:np.ndarray, r2vec:np.ndarray, dt:float, mu:float,
            prograde:bool|None=None):
    '''solve lamberts problem between two input vectors an a time,
    will choose the short angle by default, however if prograde is true
    it will enforce a prograde orbit (ccw around z-axis).
    likewise if prograde is false, enforce a retrograde orbit.
    currently only does one revolution solutions'''


    r1 = np.linalg.norm(r1vec)
    r2 = np.linalg.norm(r2vec)
    angle = m.acos((r1vec.dot(r2vec))/(r1*r2))

    if not prograde is None: # figure out angle
        wvec = np.cross(r1vec,r2vec)
        if wvec[2] < 0: angle *= -1 # long way is prograde
        if not prograde: angle *= -1 # want to go retrograde
    f = angle % (2*np.pi) # invert if long way

    # canonise for ease of convergence:
    DU = r1/2 + r2/2 # Distance unit is average of magnitudes
    TU = m.sqrt(DU**3/mu) # Time unit from Kepler's 3rd law
    r1 /= DU; r2 /= DU; dt /= TU; mu = 1 # apply canonization

    # precalculation:
    A = m.sin(f) * m.sqrt(r1*r2 / (1-m.cos(f))) 
    
    # function definitions:
    S = stumpff_s
    C = stumpff_c
    def y(z,S,C):
        try:
            return 2 + A * (z*S-1) / m.sqrt(C) # since r1+r2 = 2 [DU]
        except ZeroDivisionError, ValueError:
            return 2 # when S & C hits 0, best to tread 2nd term as 0

    def F(z):
        Sz = S(z)
        Cz = C(z)
        yz = y(z,Sz,Cz) # calculating before saves a lot of calls
        return (yz/Cz)**(3/2) * Sz + A * m.sqrt(yz) - dt # since dt*sqrt(mu) = dt
    
    a = _yzval((r1+r2)/A) # lower bound
    b = 39.477 # (4pi^2) upper bound
    try:
        z = root_finder_bisection(F,a,b)
    except ValueError: raise ArithmeticError("Lambert failed to find root, possibly too short time")
    
    # use result to find lagrange coefficients:
    Sz = S(z)
    Cz = C(z)
    f = 1 - y(z,Sz,Cz)/r1
    g_dot = 1 - y(z,Sz,Cz)/r2
    g = A*m.sqrt(y(z,Sz,Cz)) * TU # undo canonisation (rest are unitless so dont need)

    v1vec = (r2vec - f*r1vec)/g
    v2vec = (g_dot*r2vec - r1vec)/g
    return v1vec, v2vec


def _yzval(d):
    '''calculate approximate root of lambert y function
    from values of d = r1+r2 / A'''
    # uses a linear interpolation, 
    # so overestimates it slightly,
    # meaning that there is a slight window between (z_min,z_min_approx)
    # that is not covered. if the value is there it will miss it.
    # in future make this function closer to get better results for that
    yarr = [39.478, 9.869, 0,    -27.965, -44.651, -56.166, -64.706]
    darr = [-1.414, 0,     1.414, 10,     20,       30,      39.478]
    y = np.interp(d,darr,yarr)
    return y