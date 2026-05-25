'''
basic hohmann functions and the like for simple transfers
(those often form the basis of mroe complicated analysis)

'''

import math as m
from ..utils import *



__all__ = [
    'hohmann_transfer',
    'synodic_period',
    'hohmann_angle',
    
]



def hohmann_transfer(r1:float, r2:float ,mu:float)->tuple[float, float, float]:
    '''calculate delta v and travel time for a hohmann transfer,
    return dv1, dv2, travel time'''
    
    dv1 = m.sqrt(mu/r1)*abs( # first burn
        m.sqrt(2*r2/(r1+r2)) - 1)
    dv2 = m.sqrt(mu/r1)*abs( # second burn
        m.sqrt(2*r1/(r1+r2)) - 1)
    T = a2T(apse2ae(r1,r2)[0],mu) # tranfer orbit period
    return dv1,dv2, T/2


def synodic_period(T1:float, T2:float)->float:
    '''synodic period of two periods'''
    return (T1*T2)/abs(T1-T2)

def hohmann_angle(r1:float, r2:float, mu:float)->float:
    '''ideal angle a planet is ahead to start hohmann transfer.
    tranfer time can be provided to prevent redundant computation'''

    T_tra = a2T(0.5*(r1+r2),mu)/2
    return m.pi - T2n(a2T(r1,mu))*T_tra
    


