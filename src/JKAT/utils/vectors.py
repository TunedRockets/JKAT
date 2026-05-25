''' 
equations for going from and to vectors from elements
that is, keplerian elements <-> cartesian vectors
see `elements.py` for list of relevant variables


'''

import math as m
import numpy as np
from numpy.linalg import norm
from typing import overload

from .elements import finf, f2r
from .misc import stumpff_c, stumpff_s
from .optimizers import root_finder_newton

__all__ = [
    'Q_basis',
    'kep2rvec',
    'kep2vvec',
    'kep2vectors',
    'vectors2kep',
    'propagate_vectors'
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


def kep2rvec(f:float|np.ndarray,
             p:float,
             e:float,
             i:float,
             raan:float,
             argp:float)->np.ndarray:
    '''position vector in inertial frame'''

    if isinstance(f,np.ndarray):
        # truncate invalid angles:
        if e > 1: f = np.maximum(np.minimum(f, finf(e)),-finf(e))
        rr = f2r(f,p,e)
        rr = np.vstack([rr*np.cos(f), rr*np.sin(f), np.zeros_like(f)])
        rr = Q_basis(raan,i,argp)@rr
        return rr.T
    else:
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
    r = norm(rvec) # radius
    h = norm(hvec) # angular momentum
    evec = np.cross(vvec,hvec)/mu - rvec/r # eccentricity vector
    e:float = norm(evec) # type:ignore
    p = hvec.dot(hvec)/mu # parameter
    vr = rvec.dot(vvec)/r # radial velocity
    k = np.array([0,0,1]) # "up" vector
    nvec = np.cross(k,hvec) # node vector

    # normal calculations:
    i = m.acos(hvec[2]/h) # inclination

    raan = m.acos(nvec[0]/norm(nvec)) # RAAN
    if nvec[1] < 0: raan = 2*m.pi - raan # angle is negative (but we want to keep it in range 0-2pi)
    
    argp = m.acos(nvec.dot(evec)/(norm(nvec)*e)) # argument of periapsis
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
        f = m.acos(nvec.dot(rvec)/(r*norm(nvec)))
        if rvec[2] < 0: f = 2*m.pi - f
    elif e==0 and i== 0: # circular equatorial
        # node and periapsis on x axis
        raan = 0
        argp = 0
        f = m.acos(rvec[0]/r)
        if rvec[1]<0: f = 2*m.pi - f

    return p,e,i,raan,argp,f


def propagate_vectors(rvec:np.ndarray, vvec:np.ndarray, dt:float, mu:float)->tuple[np.ndarray,np.ndarray]:
    '''
    Propagate a r and v vector forward in time by dt\n
    Uses the universal variable method, so does not care about the type of orbit
    (barring maybe degenerate ones)\n
    this is essentially equivalent to from_rv().time_to_rv().
    '''

    # using Curtis' numerical method
    r = norm(rvec)
    vr = np.dot(rvec, vvec)/r # initial radial velocity

    alpha = 2/r - vvec.dot(vvec)/mu # inverse semi-major axis
    # quick almost zero check
    root_mu = m.sqrt(mu)
    chi = root_mu*dt*abs(alpha) # good first guess
    def F(chi): 
        z = alpha * chi**2
        return (r*vr/root_mu) * chi**2 * stumpff_c(z) + (1-alpha*r) * chi**3 * stumpff_s(z) + r*chi - root_mu*dt
    
    def F_prime(chi):
        z = alpha * chi**2
        return (r*vr/root_mu) * chi * (1-alpha*chi**2 * stumpff_s(z)) + (1-alpha*r) * chi**2 * stumpff_c(z) + r
    
    chi = root_finder_newton(F, F_prime, chi)

    # lagrange coefficents using chi
    f = 1 - (chi**2/r)*stumpff_c(chi**2 * alpha)
    g = dt - 1/root_mu * chi**3 * stumpff_s(chi**2 * alpha)
    r_1 = f*rvec + g*vvec
    
    f_dot = root_mu/(r * norm(r_1)) * (chi**3 * alpha * stumpff_s(chi**2 * alpha) - chi)
    g_dot = 1 - chi**2/norm(r_1) * stumpff_c(chi**2 * alpha)
    v_1 = f_dot*rvec + g_dot*vvec

    return r_1, v_1