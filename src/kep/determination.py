'''
This is functions that generate orbits from ephimerides, elements, or observations
each function will return an orbit object.
methods for finding optimal transfers are not included, just generation of orbits
'''
from .orbits import Orbit
from ..utils import stumpff_c, stumpff_s, unit, a2p,\
p2h, root_finder_bisection, L2M, M2t, f2t, apse2ae
import math as m
import numpy as np

__all__ = [

]

def orbit_from_ephemeris(a:float, e:float, i:float, L:float, longp:float, raan:float, mu:float)->Orbit:
    '''
    creates an orbit from ephemeris values, only works with elliptical orbits
    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param i: inclination
    :type i: float
    :param L: mean longitude
    :type L: float
    :param longp: longitude of perihelion
    :type longp: float
    :param raan: longitude of ascending node
    :type raan: float
    :param mu: Parent body standard gravitational parameter
    :type mu: float
    '''
    if e >= 1: raise ValueError("Hyperbolic ephemeris not implemented")
    p = a2p(a,e)
    M = L - longp # mean anomaly
    argp = longp - raan
    M = L2M(L,raan,argp)
    tp = -M2t(M,e,p2h(p,mu),mu)
    return Orbit(p,e,i,raan,argp,tp,mu)

def orbit_from_keplerian(a:float, e:float, i:float, raan:float, argp:float, f:float, mu:float)->Orbit:
    '''
    Generate an orbit from the 6 keplerian elements. true anomaly assumed to be
    the anomaly at the epoch, (from which the periapsis passage is calculated)

    :param a: Semi-major axis
    :type a: float
    :param e: Eccentricity
    :type e: float
    :param i: Inclination
    :type i: float
    :param raan: Right ascension of the ascending node
    :type raan: float
    :param arg_p: Argument of periapsis
    :type arg_p: float
    :param theta: True anomaly
    :type theta: float
    :param sgp: Parent body standard gravitational parameter
    :type sgp: float
    '''
    if (e < 1 and a < 0) or (e > 1 and a > 0) or (e == 1 and a != m.inf): raise ValueError("Eccentricity and Semi-major axis sign mismatch") 
    p = a2p(a,e)
    tp = -f2t(f,e,p2h(p,mu),mu)
    ob = Orbit(p,e,i,raan,argp,tp,mu)
    return ob

def from_elements(*,p:float|None=None, e:float|None=None, a:float|None=None,
                  ap:float|None=None, pe:float|None=None, h:float|None=None,
                  mu:float|None=None,...)->Orbit:
    '''create orbit from arbitrary list of elements'''

    # geometry:
    if not ap is None and not pe is None: # apses
        a_, e_ = apse2ae(ap,pe)
        if not a is None: a = a_
        else: raise ValueError("Overconstrained definintion (ap, pe, a)")
        if not e is None: e = e_
        else: raise ValueError("Overconstrained definintion (ap, pe, e)")
        p_ = a2p(a,e)
        if not p is None: p = p_
        else: raise ValueError("Overconstrained definintion (ap, pe, p)")
    elif not ap is None: # only apoapsis
        ...
    #




def orbit_from_rv(r:np.ndarray, v:np.ndarray, sgp:float, time:float=0)->Orbit:
    '''
    Creates an orbit given a position and velocity vector.\n
    the time is given to set the orbit within the epoch, if nothing is provided then
    it's assumed the given data is at the epoch.
    '''

    if np.linalg.norm(np.cross(r,v)) == 0:
        raise NotImplementedError("Degenerate orbits are not implemented yet")

    h_vec = np.cross(r,v) # angular momentum vector
    r_mag = np.linalg.norm(r)
    v_r = r.dot(v)/r_mag
    h = np.linalg.norm(h_vec) # angular momentum
    e_vec = np.cross(v,h_vec)/sgp - r/r_mag # eccentricity vector
    e = np.linalg.norm(e_vec) # eccentricty
    p = h_vec.dot(h_vec)/sgp # parameter
    k_hat = np.array([0,0,1]) # "up" vector
    n_vec = np.cross(k_hat,h_vec) # node vector

    # normal calculations:
    i = m.acos(h_vec[2]/h) # inclination

    RAAN = m.acos(n_vec[0]/np.linalg.norm(n_vec)) # RAAN
    if n_vec[1] < 0: RAAN = 2*m.pi - RAAN # angle is negative (but we want to keep it in range 0-2pi)
    
    arg_p = m.acos(n_vec.dot(e_vec)/(np.linalg.norm(n_vec)*e)) # argument of periapsis
    if e_vec[2] < 0: arg_p = 2*m.pi - arg_p

    theta = m.acos(e_vec.dot(r)/(e*r_mag)) # true anomaly
    if v_r < 0: theta = 2*m.pi - theta

    # special cases:
    if i==0 and e!= 0: # elliptical equitorial
        RAAN = 0 # "ascending node" on x axis
        arg_p = m.acos(e_vec[0]/e) # arg_p from x axis
        if e_vec[1] < 0: arg_p = 2*m.pi - arg_p

    elif e == 0 and i != 0: # circular inclined
        
        arg_p = 0 # "periapsis" on node
        theta = m.acos(n_vec.dot(r)/(r_mag*np.linalg.norm(n_vec)))
        if r[2] < 0: theta = 2*m.pi - theta
    elif e==0 and i== 0: # circular equatorial
        # node and periapsis on x axis
        RAAN = 0
        arg_p = 0
        theta = m.acos(r[0]/r_mag)
        if r[1]<0: theta = 2*m.pi - theta

    ob = Orbit(p,e,i,RAAN,arg_p,0,sgp) # type: ignore (for e, which is "floating", not float)
    # now figure out the passage at periapsis:
    ob.link_time_and_theta(theta, time) # i swear to god you must work!
    return ob

def orbit_from_lambert(r1:np.ndarray, r2:np.ndarray, start_time:float,
                        end_time:float, sgp:float, short_way:bool = True)->Orbit:
    '''creates an orbit by solving lambert's problem\n
    start_time is the time at r1, end_time is the time at r2. can choose between long and short way.\n
    for getting only the vectors, use lambert_vectors()\n
    start_time also used to set orbit in proper epoch'''
    v1, _ = lambert_vectors(r1,r2,(end_time-start_time),sgp,short_way)
    ob = orbit_from_rv(r1,v1, sgp, start_time)
    return ob

def orbit_from_lambert_transfer(origin:Orbit, destination:Orbit, start_time:float,
                        end_time:float, short_way:bool = True)->Orbit:
    r1 = origin.time_to_rv(start_time)[0]
    r2 = destination.time_to_rv(end_time)[0]
    return orbit_from_lambert(r1,r2,start_time,end_time,origin.sgp,short_way)
    

def orbit_from_gauss(observations:list[np.ndarray],
                        times:list[float], 
                        positions:list[np.ndarray],
                        sgp:float)->"Orbit":
    '''Uses Gauss' method to find an orbit from 3 observations, taken at 3 different times.
    you also need to include the positions of the three observations at those 3 given times
    observations assumed to be normal vectors'''
    # TODO UNTESTED!!!!!

    # Basic checks before we start the process:
    if len(observations) != 3 or len(times) != 3 or len(positions) != 3:
        raise ValueError(f"Incorrect number of arguments for Gauss' method"
                            f"({len(observations)}/3) observations"
                            f"({len(times)}/3) times"
                            f"({len(positions)}/3) positions")
    
    
    # unpack the values
    rho1, rho2, rho3 = observations
    t1, t2, t3 = times
    R1, R2, R3 = positions

    #ensure normalised observations
    rho1 = unit(rho1)
    rho2 = unit(rho2)
    rho3 = unit(rho3)

    # calculate intermediates
    tau1 = t1-t2
    tau3 = t3-t2
    tau = t3-t1
    p1 = np.cross(rho2, rho3)
    p2 = np.cross(rho1, rho3)
    p3 = np.cross(rho1, rho2)
    D0 = np.cross(rho1, p1)

    #calculate the D matrix (0 indexed unlike the Curtis variables)
    D = np.outer(np.array([R1,R2,R3]),np.array([p1,p2,p3]))
    
    # if anything is wrong, check D first

    # calculate more intermediates
    A = 1/D0 * (-D[0,1]*tau3/tau + D[1,1] + D[2,1]*tau1/tau)
    B = 1/(6*D0) * (D[0,1]*(tau3**2 - tau**2)*tau3/tau + D[2,1]*(tau**2 - tau1**2)*tau1/tau)
    E = R2.dot(rho2)
    a = -(A**2 + 2*A*E + R2**2)
    b = -2*sgp*B*(A + E)
    c = -sgp**2 * B**2

    # find root
    root = np.roots([1,0,a,0,0,b,0,0,c])[0]
    # in theory there can be multiple roots, but we simply pick the 1st
    r2 = root

    # find slant ranges
    range1 = 1/D0 * (
        (6*(D[2,0]*tau1/tau3 + D[1,0]*tau/tau3)*r2**3 + sgp*D[2,0]*(tau**2 - tau1**2)*tau1/tau3) /
        (6*r2**3 + sgp*(tau**2 - tau3**2))
        - D[0,0]
    )
    range2 = A + sgp*B/(r2**3)
    range3 = 1/D0 * (
        (6*(D[0,2]*tau3/tau1 + D[1,2]*tau/tau1)*r2**3 + sgp*D[0,2]*(tau**2 - tau3**2)*tau3/tau1) /
        (6*r2**3 + sgp*(tau**2 - tau1**2))
        - D[2,2]
    )

    # find positions (overwriting r2)
    r1 = R1 + range1*rho1
    r2 = R2 + range2*rho2
    r3 = R3 + range3*rho3

    # find lagrange coefficients
    f1 = 1 - 1/2 * sgp/r2**3 * tau1**2
    g1 = tau1 - 1/6 * sgp/r2**3 * tau1**3
    f3 = 1 - 1/2 * sgp/r2**3 * tau3**2
    g3 = tau3 - 1/6 * sgp/r2**3 * tau3**3

    # find v2
    v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1)

    # refine solution using universal anomaly
    chi_l_bound = -np.pi
    chi_u_bound = np.pi # ??? TODO
    error = 1
    while (error > 1e-8):

        old_v2 = v2
        old_r2 = r2

        magr2 = np.linalg.norm(r2)
        a = 1/(2/magr2 - (v2.dot(v2)/sgp))
        vr2 = v2.dot(unit(r2))

        # solve for universal anomaly
        # let's try bisection like in lambert vectors
        def fn(x): return (magr2*vr2/(np.sqrt(sgp)) * x**2 * stumpff_c(x**2/a) +
                            (1-magr2/a) * x**3 * stumpff_s(x**2/a) + magr2*x)
        def fn1(x): return fn(x) - np.sqrt(sgp)*tau1
        def fn3(x): return fn(x) - np.sqrt(sgp)*tau3

        chi1 = root_finder_bisection(fn1, chi_l_bound, chi_u_bound)
        chi3 = root_finder_bisection(fn3, chi_l_bound, chi_u_bound)

        # get new lagrange coefficients
        f1 = 1 - chi1**2/magr2 * stumpff_c(chi1**2/a)
        g1 = tau1 - 1/np.sqrt(sgp) * chi1**3 * stumpff_s(chi1**2/a)
        f3 = 1 - chi3**2/magr2 * stumpff_c(chi3**2/a)
        g3 = tau3 - 1/np.sqrt(sgp) * chi3**3 * stumpff_s(chi3**2/a)

        # get new slant ranges with this (using intermediates)
        c1 = g3/(f1*g3 - f3*g1)
        c3 = -g1/(f1*g3 - f3*g1)
        range1 = 1/D0 * (-D[0,0] + D[1,0]/c1 - c3/c1 * D[2,0])
        range2 = 1/D0 * (-c1*D[0,1] + D[1,1] - c3 * D[2,1])
        range3 = 1/D0 * (-c1/c3*D[0,2] + D[1,2]/c3 - D[2,2])

        # find positions again
        r1 = R1 + range1*rho1
        r2 = R2 + range2*rho2
        r3 = R3 + range3*rho3

        # find v2
        v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1)

        # calculate error:
        error = np.linalg.norm(old_v2 - v2) + np.linalg.norm(old_r2 - r2)
        # not the most theoretically sound error but it works (hopefully)
        continue
    # now we return the improved measurements
    return orbit_from_rv(r2, v2, sgp, t2)

def point_to_point(p1:np.ndarray, p2:np.ndarray, radius:float, start_time:float, end_time:float, sgp:float, angular_speed:float, epoch_angle:float)->Orbit:
    '''Create a point to point orbit between two coordinates on a sphere with given radius\n
    points given in elevation/azimuth.\n
    essentially a wrapper of the lambert function'''

    long_shift_1 = epoch_angle + start_time*angular_speed
    long_shift_2 = epoch_angle + end_time*angular_speed
    r1 = elaz_vector(p1[0], p1[1] + long_shift_1, radius)
    r2 = elaz_vector(p2[0], p2[1] + long_shift_2, radius)
    v1,_ = lambert_vectors(r1,r2,end_time-start_time,sgp,True)
    # ensure orbit isn't inside the planet by checking we go upwards:
    if v1.dot(r1) <= 0: raise ValueError("too short time given")

    return orbit_from_rv(r1,v1, sgp, start_time)
