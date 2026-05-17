''' 
equations for going from and to vectors from elements
that is, keplerian elements <-> cartesian vectors
see `elements.py` for list of relevant variables


'''

import math as m
import numpy as np

from .elements import finf, f2r

__all__ = [
    'Q_basis',
    'kep2rvec',
    'kep2vvec',
    'kep2vectors',
    'vectors2kep',
]


def Q_basis(raan:float, i:float, argp:float)->np.ndarray:
    '''get perifocal frame basis matrix,
    i.e. [p,q,w] basis vectors'''

    # normalize angles:
    if i < 0:
        i *= -1
        raan += m.pi
        argp += m.pi

    co = m.cos(raan)
    so = m.sin(raan)
    ci = m.cos(i)
    si = m.sin(i)
    cp = m.cos(argp)
    sp = m.sin(argp)

    Q = np.array([
        [co*cp - so*sp*ci, -co*sp -so*cp*ci, so*si],
        [so*cp+co*sp*ci, -so*sp + co*cp*ci, -co*si],
        [sp*si, cp*si, ci]
    ])
    return Q

def kep2rvec(f:float,
             p:float,
             e:float,
             i:float,
             raan:float,
             argp:float)->np.ndarray:
    '''position vector in inertial frame'''


    # check invalid angles:
    if e > 1 and abs(f) > finf(e): raise ArithmeticError(
        f"Invalid true anomaly for hyperbolic orbit. {abs(f)} > {finf(e)}"
    )

    rvec = f2r(f,p,e) * np.array([m.cos(f), m.sin(f), 0]) # rvec in pqw
    rvec = Q_basis(raan,i,argp)@rvec
    return rvec

def kep2vvec(f:float,
             p:float,
             e:float,
             h:float,
             i:float,
             raan:float,
             argp:float)->np.ndarray:
    '''velocity vector in inertial frame'''

    vvec = h/p * np.array([-m.sin(f),e + m.cos(f), 0]) # vvec in pqw
    vvec = Q_basis(raan,i,argp)@vvec

    return vvec

def kep2vectors(f:float,
                p:float,
                e:float,
                h:float,
                i:float,
                raan:float,
                argp:float)->tuple[np.ndarray,np.ndarray]:
    '''position and velocity vector in inertial frame'''
    
    # check invalid angles:
    if e > 1 and abs(f) > finf(e): raise ArithmeticError(
        f"Invalid true anomaly for hyperbolic orbit. {abs(f)} > {finf(e)}"
    )
    Q = Q_basis(raan,i,argp)
    rvec = f2r(f,p,e) * np.array([m.cos(f), m.sin(f), 0]) # rvec in pqw
    rvec = Q@rvec

    vvec = h/p * np.array([-m.sin(f),e + m.cos(f), 0]) # vvec in pqw
    vvec = Q@vvec

    return rvec, vvec

def vectors2kep(rvec:np.ndarray,vvec:np.ndarray,
           mu:float)->tuple[float,float,float,float,float,float]:
    '''keplerian elements from vectors.

    :return: p, e, i, raan, argp, f
    :rtype: tuple[float,float,float,float,float,float]'''

    hvec = np.cross(rvec,vvec) # angular momentum vector
    r = np.linalg.norm(rvec) # radius
    h = np.linalg.norm(hvec) # angular momentum
    evec = np.cross(vvec,hvec)/mu - rvec/r # eccentricity vector
    e:float = np.linalg.norm(evec) # type:ignore
    p = hvec.dot(hvec)/mu # parameter
    vr = rvec.dot(vvec)/r # radial velocity
    k = np.array([0,0,1]) # "up" vector
    nvec = np.cross(k,hvec) # node vector

    # normal calculations:
    i = m.acos(hvec[2]/h) # inclination

    raan = m.acos(nvec[0]/np.linalg.norm(nvec)) # RAAN
    if nvec[1] < 0: raan = 2*m.pi - raan # angle is negative (but we want to keep it in range 0-2pi)
    
    argp = m.acos(nvec.dot(evec)/(np.linalg.norm(nvec)*e)) # argument of periapsis
    if evec[2] < 0: argp = 2*m.pi - argp

    f = m.acos(evec.dot(rvec)/(e*r)) # true anomaly
    if vr < 0: f = 2*m.pi - f

    # special cases:
    if i==0 and e!= 0: # elliptical equitorial
        raan = 0 # "ascending node" on x axis
        argp = m.acos(evec[0]/e) # arg_p from x axis
        if evec[1] < 0: argp = 2*m.pi - argp

    elif e == 0 and i != 0: # circular inclined
        
        argp = 0 # "periapsis" on node
        f = m.acos(nvec.dot(rvec)/(r*np.linalg.norm(nvec)))
        if rvec[2] < 0: f = 2*m.pi - f
    elif e==0 and i== 0: # circular equatorial
        # node and periapsis on x axis
        raan = 0
        argp = 0
        f = m.acos(rvec[0]/r)
        if rvec[1]<0: f = 2*m.pi - f

    return p,e,i,raan,argp,f