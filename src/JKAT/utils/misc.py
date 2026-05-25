''' 
Non-astrodynamical functions, such as mathematical ones
'''
import math as m
import numpy as np

# for star importing:
__all__ = [
    "stumpff_s",
    "stumpff_c",
    "unit",
    'elazr2vec',
    'vec2elazr',
    'rodrigues_rot',
    'mod_bounds',
]


# === Stumpff ====
def stumpff_s(z:float)->float:
    '''Stumpff sine function, also known as c_3,
    equivalent to the infinite series:
    ∑(-z)^n / (2n + 3)! for n->infinity

    :param z: input value
    :type z: float
    :return: S(z)
    :rtype: float
    '''
    if z==0:
        return 1/6
    elif z > 0:
        z_sqrt = m.sqrt(z)
        return (z_sqrt - m.sin(z_sqrt))/(z_sqrt**3)
    else:
        z_sqrt = m.sqrt(-z)
        return (-z_sqrt + m.sinh(z_sqrt))/(z_sqrt**3)
    
def stumpff_c(z:float)->float:
    '''Stumpff sine function, also known as c_2,
    equivalent to the infinite series:
    ∑(-z)^n / (2n + 2)! for n->infinity

    :param z: input value
    :type z: float
    :return: C(z)
    :rtype: float
    '''

    if z==0:
        return 1/2
    elif z > 0:
        return (1-m.cos(m.sqrt(z)))/z
    else:
        return (-1 + m.cosh(m.sqrt(-z)))/(-z)
    

# === linear algebra ===

def unit(x:np.ndarray)->np.ndarray:
    '''creates a normal vector'''
    return x/np.linalg.norm(x)

def elazr2vec(el:float, az:float, r:float)->np.ndarray:
    '''elevation, azimuth, and range to a vector'''
    return np.array([np.cos(el) * np.cos(az),
        np.cos(el) * np.sin(az),
        np.sin(el)])*r

def vec2elazr(rvec:np.ndarray)->tuple[float,float,float]:
    '''vector to elevation, azimuth, and range'''
    Az = np.arctan2(rvec[1], rvec[0])
    El = np.pi/2 - np.arctan2(np.sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1]),rvec[2])
    z = float(np.linalg.norm(rvec))
    return El, Az, z

def rodrigues_rot(v:np.ndarray, pivot:np.ndarray, angle:float)->np.ndarray:
    '''rotates one vector around another using Rodrigues' formula, follows 
    right hand rule/CCW convention
    see also: rot_unit'''

    # normalise pivot:
    pivot = unit(pivot)
    return v*np.cos(angle) + (np.cross(pivot,v))*np.sin(angle) + pivot*(pivot.dot(v))*(1-np.cos(angle))

# === other math ===

def mod_bounds(lower:float, value:float, upper:float, modulo:float)->list[float]:
    '''
    Will return all values that are a modulo distance from value
    which are inside the bounds
    '''
    
    # move value below lower:
    while value-modulo > lower: value -= modulo
    # go through to add to list
    l = []
    while value <= upper:
        if lower <= value < upper: l.append(value)
        value += modulo
    return l