''' 
simple equations for converting between different elements,
and getting simple resulting variables.
this way there is no trivial mistype when using these simple
equations. if this were C++ then they would be inline or smth.

relevant variables are:
p - parameter
a - semi-major axis
h - angular momentum
e - eccentricity
i -inclination
raan - right ascention of ascending node
argp - argument of periapsis
f - true anomaly
longp - longitude of periapsis
arglat - argument of latitude
l - true longitude
r - radius
v - velocity
vr - radial velocity
vt - tangential velocity
fpa - flight path angle
vesc - escape velocity
vcirc - circular velocity
vinf - excess velocity
c3 - characteristic energy
finf - hyperbolic asymptote angle
ap - apoapsis
pe - periapsis


for now, only supports single inputs, but changing to
numpy is trivial
'''
from .anomalies import *
import math as m
import numpy as np
from typing import overload

__all__ = [
    'h2p',
    'p2h',
    'a2p',
    'p2a',
    'h2a',
    'a2ap',
    'a2pe',
    'h2ap',
    'h2pe',
    'apse2ae',
    'pe2p',
    'pev2e',
    'r2f',
    'longp',
    'arglat',
    'l',
    'f2r',
    'f2v',
    'f2vt',
    'f2vr',
    'f2v',
    'f2fpa',
    'vcirc',
    'vesc',
    'vinf',
    'c3',
    'finf',
    'turn_angle',
    'aiming_radius',
    'hohmann_transfer',
    'hohmann_angle',
    'synodic_period',
    'circular_phasing'
]



# === ellipses ===

def h2p(h:float, mu:float)->float:
    '''anglular momentum to parameter'''
    return h**2/mu

def p2h(p:float, mu:float)->float:
    '''parameter to angular momentum'''
    return m.sqrt(p*mu)

def a2p(a:float,e:float)->float:
    '''semi major axis to parameter'''
    return a*(1-e**2)

def p2a(p:float, e:float)->float:
    '''parameter to semi major axis'''
    if e == 1: return m.inf
    else: return p / (1-e**2)

def h2a(h:float, e:float, mu:float)->float:
    '''angular momentum to semi major axis'''
    if e == 1: return m.inf
    else: return h**2 / (mu * (1-e**2))

def a2ap(a:float,e:float)->float:
    '''apoapsis'''
    if e < 1: return a*(1+e)
    else: return m.inf
    

def a2pe(a:float,e:float)->float:
    '''periapsis'''
    return a*(1-e)

def h2ap(h:float,e:float, mu:float):
    '''angular momentum to apoapsis'''
    return (h**2/mu) / (1-e)

def h2pe(h:float,e:float, mu:float):
    '''angular momentum to periapsis'''
    return (h**2/mu) / (1+e)

def apse2ae(ap:float, pe:float)->tuple[float,float]:
    '''semi-major axis and eccentricity from apsides'''
    return 0.5*(ap+pe), abs(ap-pe)/(ap+pe)

def pe2p(pe:float, e:float)->float:
    '''periapsis and eccentricity to parameter
    (supports hyperbolic orbits)'''
    return pe*(1+e)

# === polar equation derivatives ===

@overload
def f2r(f:np.ndarray, p:float, e:float)->np.ndarray:...
@overload
def f2r(f:float, p:float, e:float)->float:...

def f2r(f:float|np.ndarray, p:float, e:float)->float|np.ndarray:
    '''polar equation radius'''
    return p/(1+e*np.cos(f))

@overload
def r2f(r:float, p:float, e:float)->float:...
@overload
def r2f(r:np.ndarray, p:float, e:float)->np.ndarray:...

def r2f(r:float|np.ndarray, p:float, e:float)->float|np.ndarray:
    '''return true anomaly that results in given radius.
    returns nan if radius is impossible'''
    try: return np.arccos((p/r - 1)/e)
    except (ValueError): return np.nan

@overload
def f2vt(f:float, p:float, e:float, h:float)->float:...
@overload
def f2vt(f:np.ndarray, p:float, e:float, h:float)->np.ndarray:...

def f2vt(f:float|np.ndarray, p:float, e:float, h:float)->float|np.ndarray:
    '''tangential velocity'''
    return h/f2r(f,p,e)

@overload
def f2vr(f:np.ndarray, p:float, e:float, h:float)->np.ndarray:...
@overload
def f2vr(f:float, p:float, e:float, h:float)->float:...

def f2vr(f:float|np.ndarray, p:float, e:float, h:float)->float|np.ndarray:
    '''radial velocity'''
    return h/p * e * np.sin(f)

@overload
def f2v(f:float, p:float, e:float, h:float)->float:...
@overload
def f2v(f:np.ndarray, p:float, e:float, h:float)->np.ndarray:...

def f2v(f:float|np.ndarray, p:float, e:float, h:float)->float|np.ndarray:
    '''scalar velocity'''
    return np.sqrt(f2vr(f,p,e,h)**2 + f2vt(f,p,e,h)**2)

def f2fpa(f:float, e:float)->float:
    '''flight path angle'''
    return m.atan(e*m.sin(f)/(1+e*m.cos(f)))

# === alternate elements ====
def longp(raan:float, argp:float)->float:
    '''longitude of periapsis'''
    return raan+argp

def arglat(f:float,argp:float)->float:
    '''argument of latitude'''
    return f + argp

def l(f:float, raan:float, argp:float)->float:
    '''true longitude'''
    return f + raan + argp


# === hyperbolic values ====

def pev2e(pe:float,vp:float, mu:float)->tuple[float,float]:
    '''periapsis altitude and speed to eccentricity and parameter'''

    h = pe*vp
    p = h2p(h,mu)
    e = p/pe - 1
    return e, p

def vcirc(f:float, p:float, e:float, h:float)->float:
    '''circular velocity'''
    return m.sqrt((h**2/p)/f2r(f,p,e))

def vesc(f:float, p:float, e:float, h:float)->float:
    '''escape velocity'''
    return m.sqrt(2*(h**2/p)/f2r(f,p,e))
    
def vinf(e:float, h:float, mu:float)->float:
    '''excess velocity'''
    if e < 1: return m.nan
    else: return (mu/h)*m.sqrt((e**2 - 1))

def c3(e:float, h:float, mu:float)->float:
    '''characteristic energy'''
    return ((mu/h)**2 * (e**2 - 1))

def finf(e:float)->float:
    '''hyperbolic asymptote angle'''
    if e <= 1: return m.nan
    else: return m.acos(-1/e)

def turn_angle(e:float)->float:
    '''double the hyperbolic asymptote angle'''
    if e <= 1: return m.nan
    else: return 2*m.acos(-1/e)

def aiming_radius(p:float,e:float)->float:
    '''hyperbolic aiming radius'''
    if e <= 1: return m.nan
    else: return (p/(e**2-1))*m.sqrt(e**2 - 1)




# === hohmann ===



def hohmann_transfer(r1:float, r2:float ,mu:float)->tuple[float, float, float]:
    '''calculate delta v and travel time for a hohmann transfer,
    return dv1, dv2, travel time'''
    
    dv1 = m.sqrt(mu/r1)*abs( # first burn
        m.sqrt(2*r2/(r1+r2)) - 1)
    dv2 = m.sqrt(mu/r2)*abs( # second burn
        m.sqrt(2*r1/(r1+r2)) - 1)
    T = a2T(apse2ae(r1,r2)[0],mu) # tranfer orbit period
    return dv1,dv2, T/2

def synodic_period(T1:float, T2:float)->float:
    '''synodic period of two periods'''
    if T1 == T2: return m.inf
    return (T1*T2)/abs(T1-T2)

def hohmann_angle(r1:float, r2:float, mu:float)->float:
    '''ideal angle a planet is ahead to start hohmann transfer.
    tranfer time can be provided to prevent redundant computation'''

    T_tra = a2T(0.5*(r1+r2),mu)/2
    return m.pi - T2n(a2T(r1,mu))*T_tra
    
def circular_phasing(df:float, n:int, a:float, mu:float)->tuple[float, float]:
    '''change in phase (+ve is forward in orbit), number of orbits, semi-major axis, and mu
    to single delta v and time'''

    df = df/n # phase per orbit
    T = a2T(a,mu)
    dt = -(df/(2*m.pi))*T # single orbit phasing extra time
    T_new = T+dt
    a_new = T2a(T_new, mu)
    v1 = vcirc(0,a,0,p2h(a,mu))
    a_new,e = apse2ae(a, 2*a_new-a)
    v2 = f2v((0 if a_new > a else m.pi),a2p(a_new,e),e,p2h(a_new,mu))
    dv = abs(v1-v2)
    return dv, T_new*n






