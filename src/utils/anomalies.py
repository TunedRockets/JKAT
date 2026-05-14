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
from .optim import root_finder_newton

# for star importing:
__all__ = [
    'a2T',
    "h2n",
    "M2E",
    "E2M",
    "f2E",
    "E2f",
    "M2t",
    "t2M",
    "t2f",
    "f2t"
]

# ===== standard kepler =====

def a2T(a:float, mu:float)->float:
    '''semi major axis to period
    (hyperbolic orbits have infinite period)'''
    if 0 < a < m.inf: return 2*m.pi*m.sqrt(a**3/mu)
    else: return m.inf

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
        return E - e*m.sinh(E) - E
    
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
    
def M2t(M:float, h:float,e:float, mu:float)->float:
    ''' Mean anomaly to time'''
    return M / h2n(h,e,mu)

def t2M(t:float, h:float,e:float, mu:float)->float:
    '''time to mean anomaly'''
    return t * h2n(h,e,mu)

def t2f(t:float, h:float, e:float, mu:float)->float:
    '''time to true anomaly, using kepler's method'''
    M = t2M(t,h,e,mu)
    E = M2E(M,e)
    f = E2f(E,e)
    return f

def f2t(f:float, h:float, e:float, mu:float)->float:
    '''true anomaly to time, using kepler's method'''
    E = f2E(f,e)
    M = E2M(E,e)
    t = M2t(M,h,e,mu)
    return t



