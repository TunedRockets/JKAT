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

    kwargs.setdefault('ts_min',0)
    default_max_start_time = max(origin.canonical_time_period)
    kwargs.setdefault('ts_max', max(origin.T,destination.T))

        








    # first define optimizer function:
    def F(s:float,t:float)->float: # start + travel time
        if not (start_bounds[0]<s < start_bounds[1]): return m.inf
        if not (end_bounds[0] < s + t < end_bounds[1]): return m.inf # ensure we're not outside bounds
        r1,v1 = origin.t2vectors(s)
        r2,v2 = destination.t2vectors(s+t)
        try:
            vl1,vl2 = lambert(r1,r2,t,origin.mu)
        except (ArithmeticError, ValueError): return m.inf # trajectories doesn't work

        # exclude maximums:
        if np.linalg.norm(vl1-v1) > max_insertion: return m.inf
        if np.linalg.norm(vl2-v2) > max_relv: return m.inf
        if np.linalg.norm(r2)/AU > max_intercept_distance: return m.inf

        weight = float(
            np.linalg.norm(vl1-v1) * w_insertion +
            np.linalg.norm(vl2-v2) * w_relv + 
            t/DAY * w_travel_time +
            np.linalg.norm(r2)/AU * w_intercept_distance +
            (s+t)/DAY * w_intercept_time
        )
        return weight
    
    low_time = min(time_bounds)
    high_time = max(time_bounds)

    dts = (time_bounds[1] - time_bounds[0])/10
    dte = (time_bounds[3] - time_bounds[2])/10
    dt = (0.5*(dte+dts))

    pois = _times_of_interest(origin,destination,low_time,high_time)
    extremes = np.array([ # to fill out in case pois is empty (slightly inside the bounds)
        [time_bounds[0]+dts, time_bounds[2]+dte],
        [time_bounds[1]-dts, time_bounds[2]+dte],
        [time_bounds[1]-dts, time_bounds[3]-dte],
        [time_bounds[0]+dts, time_bounds[3]-dte],
        [time_bounds[0], time_bounds[2]+dte], # extra layer for when bounds are same
        [time_bounds[1]-dts, time_bounds[3]-dte],
        [time_bounds[0], time_bounds[3]-dte],
    ])
    pois = np.vstack((pois, extremes))

    # exclude based on bounds:
    pois = pois[pois[:,0] > time_bounds[0]]
    pois = pois[pois[:,0] < time_bounds[1]]
    pois = pois[pois[:,1] > time_bounds[2]]
    pois = pois[pois[:,1] < time_bounds[3]]


    # exclude times based on max:
    pois = pois[pois[:,1] < max_intercept_time]

    pois[:,1] -= pois[:,0] # make travel time

    # exclude times based on max:
    pois = pois[pois[:,1] < max_travel_time]

    
    # pick best of pois:
    dv_pois = []
    for poi in pois:
        try:
            dv = F(poi[0], poi[1])
        except (ArithmeticError, ValueError): 
            dv = m.inf
        dv_pois.append(dv)
    idx = np.argsort(dv_pois)
    dv_pois = np.array(dv_pois)[idx]
    pois = pois[idx]
    best = pois[0:5] if len(pois) >= 5 else pois
    dv_best = dv_pois[0:5] if len(dv_pois) >= 5 else dv_pois


    # prestudy those close to optimal:
    best = best[dv_best < dv_best[0]*1.5]
    bestopt = [simple_hill_descent_2d(F,x, dt/4, 10) for x in best]
    bof = []
    for bo in bestopt:
        bof.append(F(bo[0],bo[1]))
    idx = np.argsort(bof)
    bestopt = np.array(bestopt)[idx]

    # try the best:
    for p in bestopt:
        try:
            s_opt,t_opt = nelder_mead_2d(F,p,-dt/2, 1e-6, max_iter=1000) #type:ignore
            break # found one
        except (ValueError, ArithmeticError):
            # this one didn't work
            continue
    else:
        raise ArithmeticError("no trajectory could be found")

    # compute properties:
    r1,v1 = origin.time_to_rv(s_opt)
    r2,v2 = destination.time_to_rv(s_opt+t_opt)
    vl1,vl2 = lambert_vectors(r1,r2,t_opt,origin.sgp)
    return np.linalg.norm(vl1-v1), np.linalg.norm(vl2-v2), s_opt, s_opt+t_opt, np.linalg.norm(r2) # type: ignore




def _times_of_interest(origin:Orbit, destination:Orbit, lower_time:float, upper_time:float)->np.ndarray:
    '''helper functions to generate times of interest for the trajectory optimizer'''

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