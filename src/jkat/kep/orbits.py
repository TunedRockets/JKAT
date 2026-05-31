

import math as m
import numpy as np
from typing import Never,overload
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

    def __repr__(self) -> str: # always need a repr (for debugging at least)
        return f"Orbit:\n {self.p=}\n {self.e=}\n {self.i=}\n {self.raan=}\n {self.argp=}\n {self.tp=}\n {self.mu=}"
    
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

    

    @property
    def longp(self)->float:
        '''longitude of periapsis'''
        return longp(self.raan,self.argp)
    
    @property
    def longitude_periapsis(self)->float:
        '''longitude of periapsis, alias of self.longp'''
        return self.longp
    
    
    
# === polar equation derivatives ===

    @overload
    def r(self, f:np.ndarray)->np.ndarray:...
    @overload
    def r(self, f:float)->float:...

    def r(self,f:float|np.ndarray)->float|np.ndarray:
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
        return a2ap(self.a, self.e)
    
    @property
    def apoapsis(self)->float:
        '''apoapsis, alias of self.ap'''
        return self.ap
    
    @property
    def pe(self)->float:
        '''periapsis'''
        return a2pe(self.a, self.e)
    
    @property
    def periapsis(self)->float:
        '''periapsis, alias of self.pe'''
        return self.pe
    
    @ap.setter
    def ap(self,ap:float):
        self.set_apses(ap=ap)

    @pe.setter
    def pe(self,pe:float):
        self.set_apses(pe=pe)
    
    def set_apses(self, ap:float|None=None, pe:float|None=None):
        '''set new apses, one apsis may be omitted to keep the current value,
        if ap is smaller than pe, the argument of periapsis is flipped'''
        if ap is None: ap = self.ap
        if pe is None: pe = self.pe

        if ap < pe:
            ap,pe = pe,ap
            self.argp -=m.pi

        if m.isfinite(ap):
            a,e = apse2ae(ap,pe)
            self.e = e; self.a = a 
            if e == 0: self.argp=0
        else: # keep e if beyond 1, otherwise set e==1:
            if self.e < 1: self.e = 1
            p = pe2p(pe,self.e)
            self.p = p
    
    def cross_radius(self,r:float)->float:
        '''return the positive true anomaly where orbit crosses the altitude.
        returns NaN if orbit does not cross altitude'''
        return r2f(r, self.p, self.e)




    
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
        return vinf(self.e,self.h, self.mu)
    
    @property
    def c3(self)->float:
        '''characteristic energy'''
        return c3(self.e,self.h, self.mu)
    
    @property
    def finf(self)->float:
        return finf(self.e)
    
    @property
    def turn_angle(self)->float:
        return turn_angle(self.e)
    
    @property
    def aiming_radius(self)->float:
        return aiming_radius(self.p,self.e)


    @property
    def canonical_time_period(self)->float:
        '''representative time span of an orbit, for elliptical orbits it's the period,
        for hyperbolic orbits, it's the time between +-pi/2 true anomaly'''
        return p2char(self.p,self.e,self.mu)



# === timing ===

    def f(self, t:float = 0, after_periapsis:bool=False)->float:
        '''true anomaly at a given time, either from epoch (default) or periapsis.'''
        if not after_periapsis: t -= self.tp
        return t2f(t,self.e,self.h,self.mu)
    
    @property
    def true_anomaly(self)->float:
        '''true anomaly at epoch, alias of self.f(self.tp)'''
        return self.f(self.tp)
    
    def l(self,t:float = 0, after_periapsis:bool=False)->float:
        '''true longitude at a given time, either from epoch (default) or periapsis.'''
        if not after_periapsis: t -= self.tp
        return l(self.f(t),self.raan,self.argp)

    @property
    def true_longitude(self)->float:
        '''true longitude at epoch, alias of self.l(self.tp)'''
        return self.l(self.tp)
    
    def t(self,f:float=0, after_periapsis:bool=False)->float:
        '''time at given true anomaly, either from epoch (default) or periapsis.'''
        t = f2t(f,self.e,self.h,self.mu)
        if not after_periapsis: t += self.tp
        return t

    @property
    def T(self)->float:
        '''period of the orbit'''
        return a2T(self.a,self.mu)
    
    @T.setter
    def T(self, T)->None:
        self.a = T2a(T,self.mu)
    
    @property
    def period(self)->float:
        '''period of the orbit. alias of self.T'''
        return self.T
    

    def t2M(self,t:float, after_periapsis:bool=False)->float:
        '''time to mean anomaly, either from epoch (default) or periapsis.'''
        if not after_periapsis: t -= self.tp
        return t2M(t,self.e,self.h,self.mu)
    
    def M2t(self,M:float, after_periapsis:bool=False)->float:
        '''mean anomaly to time, either from epoch (default) or periapsis.'''
        t = M2t(M,self.e,self.h,self.mu)
        if after_periapsis: t += self.tp
        return t

    def f2M(self,f:float)->float:
        '''true anomaly to mean anomaly'''
        E= f2E(f,self.e)
        return E2M(E,self.e)

    def t2L(self, t:float, after_periapsis:bool=False)->float:
        '''time to mean longitude, either from epoch (default) or periapsis.'''
        if not after_periapsis: t -= self.tp
        return t2L(t,self.e,self.h,self.mu,self.raan,self.argp)

    def t2X(self,t:float, after_periapsis:bool=False)->float:
        '''time to universal anomaly, either from epoch (default) or periapsis.'''
        if not after_periapsis: t -= self.tp
        return t2X(t,self.e,self.h,self.mu)
    
    def f2X(self, f:float)->float:
        '''true anomaly to unversal anomaly'''
        return f2X(f,self.e,self.h,self.mu)

    @property
    def n(self)->float:
        '''mean motion of orbit'''
        return h2n(self.h,self.e,self.mu)

    def link_tf(self,t:float,f:float)->float:
        '''link a certain time and true anomaly by changing tp'''
        t_after = self.t(f,True)
        self.tp = t - t_after
        return self.tp

# === vectors ===

    def rvec(self, f:float|np.ndarray)->np.ndarray:
        '''position vector at given true anomaly'''
        return kep2rvec(f,self.p,self.e,self.i,self.raan,self.argp)
    
    def t2rvec(self,t:float)->np.ndarray:
        '''position vector at given epoch time'''
        if m.isfinite(self.T): t %= self.T
        f = self.f(t)
        return self.rvec(f)


    def vvec(self, f:float)->np.ndarray:
        '''velocity vector at given true anomaly'''
        return kep2vvec(f,self.p,self.e,self.h,self.i,self.raan,self.argp)
    
    def t2vvec(self, t:float)->np.ndarray:
        '''velocity vector at given epoch time'''
        if m.isfinite(self.T): t %= self.T
        f = self.f(t)
        return self.vvec(f)


    def vectors(self,f:float)->tuple[np.ndarray,np.ndarray]:
        '''position and velocity vectors at given true anomaly'''
        return kep2vectors(f,self.p,self.e,self.h,self.i,self.raan,self.argp)
    
    def t2vectors(self,t:float)->tuple[np.ndarray,np.ndarray]:
        '''position and velocity vectors at given epoch time'''
        if m.isfinite(self.T): t %= self.T
        f = self.f(t)
        return self.vectors(f)

    @property
    def Q_basis(self)->np.ndarray:
        '''Perifocal frame basis, i.e. [p,q,w] basis vectors'''
        return Q_basis(self.raan,self.i,self.argp)

    @property
    def evec(self)->np.ndarray:
        '''eccentricity vector'''
        return (self.Q_basis[:,0]*self.e).flatten()
    
    @property
    def hvec(self)->np.ndarray:
        '''angular momentum vector'''
        return (self.Q_basis[:,2]*self.h).flatten()
    
    @property
    def hyperbolic_origin(self)->np.ndarray:
        '''unit vector in the direction of the origin of a hyperbolic orbit'''
        return unit(rodrigues_rot(self.evec,self.hvec,self.finf))

    @property
    def hyperbolic_destination(self)->np.ndarray:
        '''unit vector in the direction of the destination of a hyperbolic orbit'''
        return unit(rodrigues_rot(self.evec,self.hvec,-self.finf))

# === other orbits ===

    def synodic_period(self,other:"Orbit")->float:
        '''period between when the angle of the two orbiting objects is the same.
        I.e. the period between launch windows (on circular orbits)

        :param other: other orbit
        :type other: Orbit
        :return: synodic period
        :rtype: float
        '''
        if self.e >=1 or other.e >= 1:
            raise ValueError("hyperbolic orbits do not have synodic periods")
        return synodic_period(self.T,other.T)
    
    def hohmann_angle(self,other:"Orbit")->float:
        '''calculates the optimal angle the other planets needs to be ahead
        for optimal hohmann transfer.

        :param other: other orbit
        :type other: Orbit
        :return: optimal angle
        :rtype: float
        '''

        if self.e >=1 or other.e >= 1:
            raise ValueError("hyperbolic orbits do not have hohmann transfers")
        return hohmann_angle(self.a,other.a,self.mu)
    
    def hohmann_time(self,other:"Orbit")->float:
        '''compute the optimal time for a hohmann transfer,
        assumes circular orbits so will be incorrect for ellitical orbits.
        this time repeats every synodic period

        :param other: other orbit
        :type other: Orbit
        :return: time of hohmann transfer
        :rtype: float
        '''

        if self.e >=1 or other.e >= 1:
            raise ValueError("hyperbolic orbits do not have hohmann transfers")
        
        # establish base angle:
        phi_0 = other.f(0) - self.f(0)
        phi_opt = self.hohmann_angle(other)
        rel_n = other.n - self.n
        if rel_n == 0: return m.inf # same time, no relative motion

        # want to wait to when it goes from phi_0 to phi'
        # phi = phi_0 + n*t (mod 2pi)
        t_until = (phi_opt - phi_0)/rel_n
        return t_until

    def hohmann(self,other:"Orbit")->tuple[float,float,float]:
        '''hohmann transfer between orbits, return dv1,dv2,time'''

        return hohmann_transfer(self.a,other.a,self.mu)
        

    def plane_crossing(self,other:"Orbit")->float:
        '''returns the true anomaly where this orbit crosses the other orbit's plane,
        returns the "ascending node", "descending node" is pi radians opposite this one

        :param other: Other orbit
        :type other: Orbit
        :return: true anomaly of relative ascending node
        :rtype: float
        '''
        e = self.evec
        N = np.cross(self.hvec,other.hvec) # relative node vector
        theta = m.acos(N.dot(e)/(np.linalg.norm(N) * np.linalg.norm(e))) # anomaly of node vector (i.e. crossing)
        return theta
    
    def distance_to(self,other:"Orbit", time:float)->np.ndarray:
        '''return vector between two orbits at given time,
        starting at self, and going to origin

        :param other: orther orbit to compare distance to
        :type other: Orbit
        :param time: time at which to measure the range
        :type time: float
        :return: vector to the other orbit
        :rtype: np.ndarray
        '''
        return self.t2rvec(time) - other.t2rvec(time)



# === helpers ====

    def point_locus(self, f_start:float = 0, f_end:float = (2*m.pi), num:int=100):
        '''get a locus of points in 3d space'''

        # limit angle:
        if self.e > 1:
            f_start = max(f_start, -self.finf+1e-6)
            f_end = min(f_end, self.finf-1e-6)
        
        Q = self.Q_basis

        ff = np.linspace(f_start,f_end,num)

        rr = self.r(ff)
        rr = np.vstack([rr*np.cos(ff), rr*np.sin(ff), np.zeros_like(ff)])
        
        rr = Q@rr
        return rr.T
        
