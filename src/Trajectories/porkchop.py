''' 
Generate a porkchop plot for a given maneuver
'''
import src.kep as kp
from .direct import direct_transfer, lambert
import numpy as np
from numpy.linalg import norm
import math as m

__all__ = [
    'porkchop_plot'
]


def porkchop_plot(ob1:kp.Orbit, ob2:kp.Orbit,start_range:list[float]|np.ndarray, end_range:list[float]|np.ndarray,
                        prograde:bool = True, rendezvous = True, min_alt=0)->np.ndarray:
    '''calculates the porkchop plot between two orbits, assumes sgp based on the first orbit
    returns a 2d array of all Dv values, and index of the lowest one.\n
    if rendezvous is true then breaking dv is included, otherwise not\n
    if min_alt is selected, all orbits below that altitude between r1 and r2 are discarded (value inf)\n
    use lambert vector on the winning times to obtain the winning orbit'''


    

    def eval_fn(r1,v1,r2,v2,s,e)->np.floating|float:
        if s >= e: return np.nan
        try: tv1, tv2 = lambert(r1,r2,(e-s), ob1.mu, prograde=prograde)
        except ArithmeticError: return np.nan # failed to converge
        if cross_alt(r1,tv1): return np.nan
        dv1 = norm(tv1-v1)
        dv2 = norm(tv2-v2)

        return dv1 + (dv2 if rendezvous else 0)


    def cross_alt(r1:np.ndarray, tv1:np.ndarray)->bool: # does not take into account r1->ap->pe->r2
        if min_alt == 0 or r1.dot(tv1) > 0: return False
        hvec = np.cross(r1,tv1)
        evec = np.cross(tv1,hvec)/ob1.mu - r1/norm(r1)
        pe = hvec.dot(hvec)/ob1.mu / (1+norm(evec))
        return pe < min_alt
    
    # vectorization for convenience
    rr1,vv1 = np.vectorize(ob1.t2vectors,otypes=[np.ndarray,np.ndarray])(start_range)
    rr2,vv2 = np.vectorize(ob2.t2vectors,otypes=[np.ndarray,np.ndarray])(end_range)
    
    DDVV = [] # can't vectorize with 1d and 2d array (without meshgrid shenanigans)
    for r1,v1,s in zip(rr1,vv1,start_range):
        row = []
        for r2,v2,e in zip(rr2,vv2,end_range):
            row.append(eval_fn(r1,v1,r2,v2,s,e))
        DDVV.append(row)
    
    return np.array(DDVV)
