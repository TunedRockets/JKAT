''' 
Direct transfers from one orbit to another, as well as the lambert transfer 
it is based on (move this?)

current transfers:
- direct transfer
    simple orbiting body to other orbiting body optimizer
wanted transfers:
- inclination augmented direct transfer
    allow change inclination before performing transfer
- time agnostic orbit transfer
    transfer from otbit to orbit, not caring about position in final orbit


'''
from ..kep import Orbit
import math as m
import numpy as np
from ..utils import *
from ..utils import stumpff_c,stumpff_s, root_finder_bisection

# import matplotlib.pyplot as plt

__all__ = [
    'direct_transfer',
    'lambert'
]



def lambert(r1vec:np.ndarray, r2vec:np.ndarray, dt:float, mu:float,
            prograde:bool|None=None):
    '''solve lamberts problem between two input vectors an a time,
    will choose the short angle by default, however if prograde is true
    it will enforce a prograde orbit (ccw around z-axis).
    likewise if prograde is false, enforce a retrograde orbit.
    currently only does one revolution solutions'''


    r1 = np.linalg.norm(r1vec)
    r2 = np.linalg.norm(r2vec)
    angle = m.acos((r1vec.dot(r2vec))/(r1*r2))

    if not prograde is None: # figure out angle
        wvec = np.cross(r1vec,r2vec)
        if wvec[2] < 0: angle *= -1 # long way is prograde
        if not prograde: angle *= -1 # want to go retrograde
    f = angle % (2*np.pi) # invert if long way

    # canonise for ease of convergence:
    DU = r1/2 + r2/2 # Distance unit is average of magnitudes
    TU = m.sqrt(DU**3/mu) # Time unit from Kepler's 3rd law
    r1 /= DU; r2 /= DU; dt /= TU; mu = 1 # apply canonization

    # precalculation:
    A = m.sin(f) * m.sqrt(r1*r2 / (1-m.cos(f))) 
    
    # function definitions:
    S = stumpff_s
    C = stumpff_c
    def y(z,S,C):
        try:
            return 2 + A * (z*S-1) / m.sqrt(C) # since r1+r2 = 2 [DU]
        except ZeroDivisionError, ValueError:
            return 2 # when S & C hits 0, best to tread 2nd term as 0

    def F(z):
        Sz = S(z)
        Cz = C(z)
        yz = y(z,Sz,Cz) # calculating before saves a lot of calls
        return (yz/Cz)**(3/2) * Sz + A * m.sqrt(yz) - dt # since dt*sqrt(mu) = dt
    
    try:
        a = _yzval((r1+r2)/A) if A > 0 else -39.477 # lower bound

        b = 39.477 # (4pi^2) upper bound
        z = root_finder_bisection(F,a,b)
    except ValueError: raise ArithmeticError("Lambert failed to find root, possibly too short time")

    # use result to find lagrange coefficients:
    Sz = S(z)
    Cz = C(z)
    f = 1 - y(z,Sz,Cz)/r1
    g_dot = 1 - y(z,Sz,Cz)/r2
    g = A*m.sqrt(y(z,Sz,Cz)) * TU # undo canonisation (rest are unitless so dont need)

    v1vec = (r2vec - f*r1vec)/g
    v2vec = (g_dot*r2vec - r1vec)/g
    return v1vec, v2vec


def _yzval(d):
    '''calculate approximate root of lambert y function
    from values of d = r1+r2 / A'''
    # uses a linear interpolation, 
    # so overestimates it slightly,
    # meaning that there is a slight window between (z_min,z_min_approx)
    # that is not covered. if the value is there it will miss it.
    # in future make this function closer to get better results for that
    aarr = [39.478, 9.869, 0,    -27.965, -44.651, -56.166, -64.706]
    darr = [-1.414, 0,     1.414, 10,     20,       30,      39.478]
    a = np.interp(d,darr,aarr)
    return a

# origin, destination, start bound, end bound, **kwargs (weights and min-max)
# variables to consider (each has _min, _max, _w)
# dv1, dv2, ts, tt, te, r

def direct_transfer(
        origin:Orbit,
        destination:Orbit,
        bounds:tuple[float,float,float,float]|None = None,
        **kwargs,
)->dict:
    '''function for calculating the optimal transfer between two orbiting objects,
    bounds is the search area, **kwargs determine bounds and weights
    returns dict with: ts,te,dv1,dv2,r'''
    
    np.seterr(all='ignore')
    
    # set **kwargs defaults
    kwargs.setdefault('dv1_min', 0)
    kwargs.setdefault('dv1_max', m.inf)
    kwargs.setdefault('dv1_w', 1) # intercept dV
    kwargs.setdefault('dv2_min', 0)
    kwargs.setdefault('dv2_max', m.inf)
    kwargs.setdefault('dv2_w', 1) # rendezvous dV (rendezvous by default)

    default_max__time = max(origin.canonical_time_period, destination.canonical_time_period)
    kwargs.setdefault('ts_min',0)
    kwargs.setdefault('ts_max', default_max__time)
    kwargs.setdefault('ts_w', 0)
    kwargs.setdefault('te_min',0)
    kwargs.setdefault('te_max', default_max__time)
    kwargs.setdefault('te_w', 0)
    kwargs.setdefault('tt_min',0)
    kwargs.setdefault('tt_max', m.inf)
    kwargs.setdefault('tt_w', 0)

    kwargs.setdefault('r_min',0)
    kwargs.setdefault('r_max', m.inf)
    kwargs.setdefault('r_w', 0)

    kwargs.setdefault('prograde', None)

    if not bounds is None:
        kwargs['ts_min'] = bounds[0]
        kwargs['ts_max'] = bounds[1]
        kwargs['te_min'] = bounds[2]
        kwargs['te_max'] = bounds[3]
    
        
    # first define optimizer function:
    def F(st:np.ndarray)->float: # start + travel time
        s = st[0]; t = st[1]
        # time exclusions:
        if not (kwargs['ts_min'] < s < kwargs['ts_max']): return m.inf
        if not (kwargs['tt_min'] < t < kwargs['tt_max']): return m.inf
        if not (kwargs['te_min'] < s + t < kwargs['te_max']): return m.inf # ensure we're not outside bounds
        r1,v1 = origin.t2vectors(s)
        r2,v2 = destination.t2vectors(s+t)
        try:
            vl1,vl2 = lambert(r1,r2,t,origin.mu, kwargs['prograde'])
        except (ArithmeticError, ValueError): return m.inf # trajectories doesn't work

        
        dv1 = np.linalg.norm(vl1-v1)
        dv2 = np.linalg.norm(vl2-v2)
        r = np.linalg.norm(r2)

        # result exclusions:
        if not (kwargs['dv1_min'] < dv1 < kwargs['dv1_max']): return m.inf
        if not (kwargs['dv2_min'] < dv2 < kwargs['dv2_max']): return m.inf
        if not (kwargs['r_min'] < r < kwargs['r_max']): return m.inf


        weight = (
            s*kwargs['ts_w'] +
            t*kwargs['tt_w'] +
            (s+t)*kwargs['te_w'] +
            dv1*kwargs['dv1_w'] +
            dv2*kwargs['dv2_w'] +
            r*kwargs['r_w']
        )
        return weight
    
    
    # get points?
    points = _pois(origin,destination,bounds=(
        kwargs['ts_min'], kwargs['ts_max'], kwargs['te_min'], kwargs['te_max']
    ))

    points[:,1] -= points[:,0] # make travel time

    Fpoints = []
    for p in points: Fpoints.append(F(p))

    points = points[np.argsort(Fpoints)]

    dt = (kwargs['te_max'] - kwargs['te_min'])/20 + (kwargs['ts_max'] - kwargs['ts_min'])/20
    # try the best:
    for p in points:
        try:
            opt = minimizer(F,p,dt) #type:ignore
            break # found one
        except (ValueError, ArithmeticError):
            # this one didn't work
            continue
    else:
        raise ArithmeticError("no trajectory could be found")

    s_opt = opt[0]
    t_opt = opt[1]
    

    # compute properties:
    r1,v1 = origin.t2vectors(s_opt)
    r2,v2 = destination.t2vectors(s_opt+t_opt)
    vl1,vl2 = lambert(r1,r2,t_opt,origin.mu)
    return {
        "ts": s_opt,
        "te": s_opt+t_opt,
        "dv1": np.linalg.norm(vl1-v1),
        "dv2": np.linalg.norm(vl2-v2),
        'r': np.linalg.norm(r2)
    }

def _pois(origin:Orbit,
          destination:Orbit,
          bounds:tuple[float,float,float,float],
          gridsize:int = 8)->np.ndarray:
    '''helper functions to generate times of interest for the trajectory optimizer,
    gets all apses, node-crossings, and potential hohmann transfer points'''

    
    bottom_s = bounds[0]
    top_s = bounds[1]
    bottom_e = bounds[2]
    top_e = bounds[3]
    bottom_s, top_s = sorted((bottom_s, top_s))
    bottom_e, top_e = sorted((bottom_e, top_e))


    # apses and nodes
    ori_cross = origin.plane_crossing(destination)
    dest_cross = destination.plane_crossing(origin)
    if origin.e < 1:
        ori_nodes = [0,m.pi, ori_cross, ori_cross + m.pi]
    else:
        ori_nodes = [0.0]
        if abs(ori_cross) < origin.finf: ori_nodes.append(ori_cross) 

        if abs(ori_cross - m.pi) < origin.finf: ori_nodes.append(ori_cross - m.pi) 

    if destination.e < 1:
        dest_nodes = [0,m.pi, dest_cross, dest_cross + m.pi]
    else:
        dest_nodes = [0.0]
        if abs(dest_cross) < destination.finf: dest_nodes.append(dest_cross) 

        if abs(dest_cross - m.pi) < destination.finf: dest_nodes.append(dest_cross - m.pi) 
    

    starts = []
    ends = []
    for n in ori_nodes:
        starts.extend(mod_bounds(bottom_s, origin.t(n), top_s, origin.period))
    for n in dest_nodes:
        ends.extend(mod_bounds(bottom_e, destination.t(n), top_e, destination.period))

    pois = np.array(np.meshgrid(starts,ends)).T.reshape(-1,2)
    # mesh of nodes/apses

    # hohmann:
    if (origin.e < 1 and destination.e < 1):

        # TODO figure out eccentric hohmann
        # pe->ap, ap->pe, ap->ap, pe->pe ?
        # values need to be shifted a bit




        t = origin.hohmann_time(destination)
        t = np.array(mod_bounds(bottom_s,t,top_s, origin.synodic_period(destination)))
        t2 = t + origin.hohmann(destination)[2]
        poi = np.column_stack((t,t2))
        pois = np.vstack((pois,poi))

    # add edges of range:


    starts = np.linspace(bottom_s, top_s, gridsize+2)[1:-1]
    ends = np.linspace(bottom_e, top_e, gridsize+2)[1:-1]

    grid = np.array(np.meshgrid(starts,ends)).T.reshape(-1,2)

    pois = np.vstack((pois,grid)) # extra grid of points

    # pois = np.vstack((pois, np.array([ # square of the bounds a distance of dt from the edge
    #     [bottom_s, bottom_e],
    #     [bottom_s, top_e],
    #     [top_s, bottom_e],
    #     [top_s, top_e],
    # ])))

    # plt.scatter(pois[:,0],pois[:,1])

    return pois





def apse_line_rotating_direct_transfer(
        origin:Orbit,
        destination:Orbit,
        t_i_change:float,
        i_bounds:tuple[float,float]=(-m.pi,m.pi),
        bounds:tuple[float,float,float,float]|None = None,
        **kwargs)->dict:
    ''' 
    Similar to direct transfer but allows a change in rotation around the apse line before optimizer is run.
    
    '''
    raise NotImplementedError()
