''' 
Direct transfer from one orbit to another
'''
from ..kep import Orbit, lambert
import math as m
import numpy as np
from ..utils import *


__all__ = [
    'direct_transfer'
]



# origin, destination, start bound, end bound, **kwargs (weights and min-max)

# variables to consider (each has _min, _max, _w)
# dv1, dv2, ts, tt, te, r
lim_vars = ['dv1','dv2','ts','tt','te','r']
suffix_vars = ['_min','_max','_w']

def direct_transfer(
        origin:Orbit,
        destination:Orbit,
        **kwargs,
)->tuple[float,float,float,float,float]:
    
    
    # set **kwargs defaults
    kwargs.setdefault('dv1_min', 0)
    kwargs.setdefault('dv1_max', m.inf)
    kwargs.setdefault('dv1_weight', 1) # intercept dV
    kwargs.setdefault('dv2_min', 0)
    kwargs.setdefault('dv2_max', m.inf)
    kwargs.setdefault('dv2_weight', 1) # rendezvous dV (rendezvous by default)

    default_max__time = max(origin.canonical_time_period, destination.canonical_time_period)
    kwargs.setdefault('ts_min',0)
    kwargs.setdefault('ts_max', default_max__time)
    kwargs.setdefault('ts_w', 0)
    kwargs.setdefault('te_min',0)
    kwargs.setdefault('te_max', default_max__time)
    kwargs.setdefault('te_w', 0)
    kwargs.setdefault('tt_min',0)
    kwargs.setdefault('tt_max', default_max__time)
    kwargs.setdefault('tt_w', 0)

    kwargs.setdefault('r_min',0)
    kwargs.setdefault('r_max', m.inf)
    kwargs.setdefault('r_w', 0)
        








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
            vl1,vl2 = lambert(r1,r2,t,origin.mu)
        except (ArithmeticError, ValueError): return m.inf # trajectories doesn't work

        
        dv1 = np.linalg.norm(vl1-v1)
        dv2 = np.linalg.norm(vl2-v2)
        r = np.linalg.norm(r2)/AU

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
    
    
    # # get points?
    # points = _pois(origin,destination,kwargs['ts_min'],kwargs['te_max'])
    # # try the best:
    # for p in best:
    #     try:
    #         s_opt,t_opt = nelder_mead_2d(F,p,-dt/2, 1e-6, max_iter=1000) #type:ignore
    #         break # found one
    #     except (ValueError, ArithmeticError):
    #         # this one didn't work
    #         continue
    # else:
    #     raise ArithmeticError("no trajectory could be found")
    bounds = [(kwargs['ts_min'],kwargs['ts_max']),
                (kwargs['ts_min'],kwargs['ts_max'])]
    s_opt, t_opt = bboptim(F,
                           np.array([np.average(bounds[0]),np.average(bounds[1])]),
                           bounds,) # type:ignore
    

    # compute properties:
    r1,v1 = origin.t2vectors(s_opt)
    r2,v2 = destination.t2vectors(s_opt+t_opt)
    vl1,vl2 = lambert(r1,r2,t_opt,origin.mu)
    return np.linalg.norm(vl1-v1), np.linalg.norm(vl2-v2), s_opt, s_opt+t_opt, np.linalg.norm(r2) # type: ignore




def _pois(origin:Orbit, destination:Orbit, lower_time:float, upper_time:float)->np.ndarray:
    '''helper functions to generate times of interest for the trajectory optimizer'''
    raise NotImplementedError()

    # apses and nodes
    ori_cross = origin.relative_node_crossing(destination)
    dest_cross = destination.relative_node_crossing(origin)
    if origin.e < 1:
        ori_nodes = [0,m.pi, ori_cross, ori_cross + m.pi]
    else:
        ori_nodes = [0.0]
        if abs(ori_cross) < origin.asymptote_angle(): ori_nodes.append(ori_cross) 

        if abs(ori_cross - m.pi) < origin.asymptote_angle(): ori_nodes.append(ori_cross - m.pi) 

    if destination.e < 1:
        dest_nodes = [0,m.pi, dest_cross, dest_cross + m.pi]
    else:
        dest_nodes = [0.0]
        if abs(dest_cross) < destination.asymptote_angle(): dest_nodes.append(dest_cross) 

        if abs(dest_cross - m.pi) < destination.asymptote_angle(): dest_nodes.append(dest_cross - m.pi) 
    

    starts = []
    ends = []
    for n in ori_nodes:
        starts.extend(inside_modulo_bounds(lower_time, origin.theta_to_time(n), upper_time, origin.period))
    for n in dest_nodes:
        ends.extend(inside_modulo_bounds(lower_time, destination.theta_to_time(n), upper_time, destination.period))

    pois = np.array(np.meshgrid(starts,ends)).T.reshape(-1,2)
    # mesh of nodes/apses

    # hohmann:
    if (origin.e < 1 and destination.e < 1):

        # figure out eccentric hohmann
        # pe->ap, ap->pe, ap->ap, pe->pe ?




        t = origin.hohmann_time(destination)
        t = np.array(inside_modulo_bounds(lower_time,t,upper_time, origin.synodic_period(destination)))
        t2 = t + origin.hohmann_travel_time(destination)
        poi = np.column_stack((t,t2))
        pois = np.vstack((pois,poi))
    # presort invalids:
    pois = pois[pois[:,1] < upper_time]
    pois = pois[pois[:,0] < pois[:,1]]
    return pois