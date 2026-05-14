''' 
Non-astrodynamical functions, such as mathematical ones
'''
import math as m
import numpy as np

# for star importing:
__all__ = [
    "stumpff_s",
    "stumpff_c",
    "unit"
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
