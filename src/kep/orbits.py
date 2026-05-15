

import math as m
import numpy as np
from ..utils import * # change to non-star
__all__ = [
    'Orbit'
]


class Orbit():
    '''
    Main class to represent a keplerian orbit
    values of the orbit are stored as a six-tuple of keplerian elements,\n
    $(p,e,i,\\Omega,\\omega,t_p)$
    as such, a parent body must be represented in the form of a given sgp
    as a 7th element.
    Other variables are accessed through property accessors.

    Special cases, I.e. circular, planar, parabolic, and degenerate, 
    are dealt with specially:
    - Degenerate orbits are currently not implemented
    - Circular orbits assume the "periapsis" is the ascending node
    - planar orbits assume the "ascending node" is on the X axis
    - Circular planar orbits apply both of the above
    - orbits with negative inclination are possible, 
      and are interpreted with the "ascending node"
      being the descending one
    '''

    def __init__(self, p:float, e:float, i:float, raan:float,
                  argp:float, tp:float, mu:float ) -> None:
        
        self.p = p
        '''parameter (semi-latus rectum)'''
        self.e = e
        '''eccentricity'''
        self.i = i
        '''inclination'''
        self.raan = raan
        '''Right ascension of the ascending node'''
        self.argp = argp
        '''argument of periapsis'''
        self.tp = tp
        '''time of periapsis passage'''
        self.mu = mu
        '''Parent body standard gravitational parameter'''


# ======================
# properties
# ======================
# implement relevant setters when possible


    # === aliases ===
    @property
    def parameter(self)->float:
        '''Alias of self.p'''
        return self.p
    
    @property
    def eccentricity(self)->float:
        '''Alias of self.e'''
        return self.e
    
    @property
    def inclination(self)->float:
        '''Alias of self.i'''
        return self.i

    @property
    def right_ascention_ascending_node(self)->float:
        '''Alias of self.raan'''
        return self.raan
    
    @property
    def argument_periapsis(self)->float:
        '''Alias of self.argp'''
        return self.argp
    
    @property
    def time_periapsis(self)->float:
        '''Alias of self.tp'''
        return self.tp

    # === derivative properties ===
    @property
    def a(self)->float:
        '''Semi-major axis\n
        hyperbolic orbits are considered to have negative semi-major axes\n
        
        Changing this variable scales p by a corresponding amount\n
        only works if new a is same sign as old a'''
        return p2a(self.p,self.e)

    @a.setter
    def a(self, a:float)->None:
        if self.a * a > 0 and self.e != 1: 
            # don't change sign accidentally
            self.p = a2p(a,self.e)
        else: raise ValueError("Cannot flip sign of semi-major axis")
    
    @property
    def semi_major_axis(self)->float:
        '''alias of self.a'''
        return self.a

    @property
    def h(self)->float:
        '''angular momentum\n
        changing this changes the parameter of the orbit'''
        return p2h(self.p,self.mu)
    
    @h.setter
    def h(self, h:float):
        # p = h^2/mu
        self.p = h2p(h,self.mu)

    def f(self, t:float = 0)->float:
        '''true anomaly at a given time, (defaults to epoch)'''
        return t2f(t,self.e,self.h,self.mu)
    
    @property
    def true_anomaly(self)->float:
        '''true anomaly at epoch, alias of self.f(0)'''
        return self.f(0)

    @property
    def longp(self)->float:
        '''longitude of periapsis'''
        return longp(self.raan,self.argp)
    
    @property
    def longitude_periapsis(self)->float:
        '''longitude of periapsis, alias of self.longp'''
        return self.longp
    
    def l(self,t:float = 0)->float:
        '''true longitude at a given time (defaults to epoch)'''
        return l(self.f(t),self.raan,self.argp)

    @property
    def true_longitude(self)->float:
        '''true longitude at epoch, alias of self.l(0)'''
        return self.l(0)
    
    # === polar equation derivatives ===

    def r(self,f:float)->float:
        '''polar equaiton radius'''
        return f2r(f,self.p,self.e)
    
    def v(self,f:float)->float:
        '''polar equation velocity'''
        return f2v(f,self.p,self.e,self.h)
    
    def vr(self,f:float)->float:
        '''polar equation radial velocity'''
        return f2vr(f,self.p,self.e,self.h)
    
    def vt(self,f:float)->float:
        '''polar equation tangential velocity'''
        return f2vt(f,self.p,self.e,self.h)
    
    def fpa(self,f:float)->float:
        '''polar equation tangential velocity'''
        return f2fpa(f,self.e)
    
    @property
    def ap(self)->float:
        '''apoapsis'''
        return ap(self.a, self.e)
    
    @property
    def pe(self)->float:
        '''periapsis'''
        return pe(self.a, self.e)
    
    # === hyperbolic values ==== TODO: add aliases

    def vesc(self, f)->float:
        '''escape velocity at given true anomaly'''
        return vesc(f,self.p,self.e,self.h)
    
    def vcirc(self, f)->float:
        '''circular velocity at given true anomaly'''
        return vcirc(f,self.p,self.e,self.h)
    
    @property
    def vinf(self)->float:
        '''excess velocity'''
        return vinf(self.e,self.h)
    
    @property
    def c3(self)->float:
        '''characteristic energy'''
        return c3(self.e,self.h)
    
    @property
    def finf(self)->float:
        return finf(self.e)
    
    @property
    def turn_angle(self)->float:
        return turn_angle(self.e)
    
    @property
    def aiming_radius(self)->float:
        return aiming_radius(self.p,self.e)





    # === timing ===

    def T():...

    def M():...

    def L():...

    def X():...

    def n():...

    def link_tf(self,t:float,f:float)->float:...

    # === vectors ===


    

