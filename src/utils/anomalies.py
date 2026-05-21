''' 
Converts between different anomalies, uses the following convention
similar to pykep, with a2b.
anomalies are:
M - mean anomaly
E - eccentric anomaly
L - mean longitude
X - universal anomaly
f - true anomaly
t - time
T - period

also relevant:
a - semi-major axis
e - eccentricity
n - mean motion

distinction between hyperbolic,parabolic, and eccentric anomalies
are all handled by the value of e given.
(all are given as E)
"parabolic anomaly" treated as the same as Mean anomaly

for now, only supports single inputs, but changing to
numpy is trivial
'''
import math as m
from .optim import root_finder_newton, root_finder_fallback
from . import elements as elem
from .misc import stumpff_c, stumpff_s

# for star importing:
__all__ = [
    'a2T',
    'T2a',
    'p2char',
    "h2n",
    "M2E",
    "E2M",
    "f2E",
    "E2f",
    "M2t",
    "t2M",
    "t2f",
    "f2t",
    "t2f_kep",
    "f2t_kep",
    't2L',
    'L2t',
    'L2M',
    'M2L',
    't2X',
    'X2t',
    'X2f',
    'f2X'
]

# clamp f for hyperbolic checks
clamp_f = lambda f: ((f+m.pi) % (2*m.pi)) - m.pi



# ===== standard kepler =====

def a2T(a:float, mu:float)->float:
    '''semi major axis to period
    (hyperbolic orbits have infinite period)'''
    if 0 < a < m.inf: return 2*m.pi*m.sqrt(a**3/mu)
    else: return m.inf

def T2a(T:float, mu:float)->float:
    '''period to semi major axis'''
    return m.cbrt(mu * (T/(2*m.pi))**2)

def p2char(p:float,e:float, mu:float)->float:
    '''representative time span of an orbit, for elliptical orbits it's the period,
    for hyperbolic orbits, it's the time between +-pi/2 true anomaly'''
    if e < 1: return a2T(elem.p2a(p,e),mu)
    else: return 2 * f2t_kep(m.pi/2,e,elem.p2a(p,e),mu)
        



def h2n(h:float, e:float, mu:float)->float:
    '''angular momentum to mean motion'''
    if e == 1: return mu**2/(h**3)
    else: return mu**2/(h**3) * abs((1-e**2)**(3/2))

def M2E(M:float, e:float)->float:
    '''Mean anomaly to eccentric anomaly'''
    if e == 0:
        return M
    elif e == 1:
        return M
    elif e < 1:
        f = lambda E: E-e*m.sin(E) - M
        df = lambda E: 1-e*m.cos(E)
        return root_finder_newton(f,df,M)
    else:
        f = lambda F: e*m.sinh(F)-F-M
        df = lambda F: e*m.cosh(F)-1
        return root_finder_newton(f,df,M)

def E2M(E:float, e:float)->float:
    '''Eccentric anomaly to mean anomaly'''
    if e == 0:
        return E
    elif e == 1:
        return E
    elif e < 1:
        return E - e*m.sin(E)
    else:
        return e*m.sinh(E) - E
    
def f2E(f:float, e:float)->float:
    '''true anomaly to eccentric anomaly'''
    if e == 0:
        return f
    elif e == 1:
        return 0.5*m.tan(f/2)+(1/6)*m.tan(f/2)**3
    elif e < 1:
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(f/2)) % (2*m.pi)
    else:
        return 2*m.atanh(m.sqrt((e-1)/(e+1))*m.tan(f/2))

def E2f(E:float, e:float)->float:
    '''eccentric anomaly to true anomaly'''
    if e == 0:
        return E
    elif e == 1:
        z = m.cbrt(3*E + m.sqrt(1+(3*E)**2))
        return 2*m.atan(z-1/z)
    elif e < 1:
        return 2*m.atan(m.sqrt((1+e)/(1-e))*m.tan(E/2)) % (2*m.pi)
    else:
        return 2*m.atan(m.sqrt((e+1)/(e-1))*m.tanh(E/2))
    
def M2t(M:float, e:float, h:float, mu:float)->float:
    ''' Mean anomaly to time'''
    return M / h2n(h,e,mu)

def t2M(t:float, e:float, h:float, mu:float)->float:
    '''time to mean anomaly'''
    return t * h2n(h,e,mu)

def t2f_kep(t:float, e:float, h:float, mu:float)->float:
    '''time to true anomaly, using kepler's method'''
    M = t2M(t,e,h,mu)
    E = M2E(M,e)
    f = E2f(E,e)
    return f

def f2t_kep(f:float, e:float, h:float, mu:float)->float:
    '''true anomaly to time, using kepler's method'''
    if e >= 1: f = clamp_f(f)
    if e > 1 and abs(f) > elem.finf(e): raise ArithmeticError(
        f"Invalid true anomaly for hyperbolic orbit. {abs(f)} > {elem.finf(e)}"
    )
    E = f2E(f,e)
    M = E2M(E,e)
    t = M2t(M,e,h,mu)
    return t


# === alternate anomalies ====

def t2L(t:float, e:float, h:float, mu:float, raan:float, argp:float)->float:
    '''time to mean longitude'''
    M = t2M(t,e,h,mu)
    return M2L(M,raan,argp)

def L2t(L:float, e:float, h:float, mu:float, raan:float, argp:float)->float:
    '''mean longitude to time'''
    M = L2M(L,raan,argp)
    return M2t(M,e,h,mu)

def M2L(M:float,raan:float,argp:float)->float:
    '''mean anomaly to mean longitude'''
    return elem.longp(raan, argp) + M

def L2M(L:float,raan:float,argp:float)->float:
    '''mean longitude to mean anomaly'''
    return L-elem.longp(raan,argp) 

# === universal variables ===

def t2X(t:float, e:float, h:float, mu:float)->float:
    '''time to universal anomaly, assumes X = 0 at t = 0 at periapsis'''

    if t==0: return 0.0 # by definition
    # e==0 might have issues as well. test and check

    a = elem.h2a(h,e,mu)
    rp = elem.h2pe(h,e,mu)
    S = stumpff_s
    C = stumpff_c
    root_mu = m.sqrt(mu)
    F = lambda chi: (1-rp/a)*chi**3 * S(chi**2/a) + rp*chi - root_mu*t
    dF = lambda chi: (1-rp/a)*chi**2 * C(chi**2/a) + rp

    # root finding can have issues, look into what below is needed:
    # (before the min max needed to be expanded for it
    # to converge properly)
    X_max = root_mu*t/rp # min max based on prussing and conway
    X_min = root_mu*t/(h*h/(mu*(1-e))) if e < 1 else 0.0
    if X_max < X_min: X_min,X_max = X_max,X_min # swap if wrong signs

    if X_max - X_min < 1e-10: X = X_max # fix edge case where region is tiny (e~=0)
    else: X = root_finder_fallback(F,dF,X_min,X_max)
    return X
    


def X2t(X:float, e:float, h:float, mu:float)->float:
    '''universal anomaly to time, assumes X = 0 at t = 0 at periapsis'''

    if X==0: return 0.0 # by definition
    # e==0 might have issues as well. test and check

    a = elem.h2a(h,e,mu)
    rp = elem.h2pe(h,e,mu)
    S = stumpff_s
    return (
        (1 - rp/a)*X**3 * S(X**2/a) + rp*X
    )/m.sqrt(mu)


def X2f(X:float,e:float, h:float, mu:float)->float:
    '''universal anomaly to true anomaly'''

    # first find eccentric anomaly, then E2f
    if e == 1: return 2 * m.atan(m.sqrt(mu)*X/h)

    # hyperbolic and eccentric have same eq
    E = m.sqrt(mu*abs(1-e**2))/h * X
    return E2f(E, e)

def f2X(f:float, e:float, h:float, mu:float)->float:
    '''true anomaly to universal anomaly'''

    if e == 1: return h/m.sqrt(mu) * m.tan(f/2)

    E = f2E(f,e)
    return h/m.sqrt(mu*abs(1-e**2)) * E

def f2t(f:float, e:float, h:float, mu:float)->float:
    '''true anomaly to time using universal variable method'''
    if e >= 1: f = clamp_f(f)
    if e > 1 and abs(f) > elem.finf(e): raise ArithmeticError(
        f"Invalid true anomaly for hyperbolic orbit. {abs(f)} > {elem.finf(e)}"
    )
    X = f2X(f,e,h,mu)
    return X2t(X,e,h,mu)

def t2f(t:float, e:float, h:float, mu:float)->float:
    '''time to true anomaly using universal variable method'''
    X = t2X(t,e,h,mu)
    return X2f(X,e,h,mu)