''' 
mathematical optimizers, all put in one place so that they (being computationally heavy)
can be referred to another package or language if a speedup is needed

also a place to defer to scipy when relevant
'''
from typing import Callable, Sequence
import numpy as np
import math as m
# import scipy.optimize as o

# for star importing:
__all__ = [
    "root_finder_bisection",
    "root_finder_newton",
    "root_finder_fallback",
    'minimizer'
]



# === root finders ===

def root_finder_bisection(f:Callable[[float],float], lower:float, upper:float, tolerance:float = 1e-10)->float:
    '''takes a univariate function and finds the root of that function
    through recursive bisection.
    converges on a root between bounds, provided bounds are of different sign

    :param f: Function to find root of
    :type f: Callable
    :param lower: lower bound on x
    :type lower: float
    :param upper: upper bound on x
    :type upper: float
    :param tolerance: maximum allowed difference between actual root and returned x, defaults to 1e-8
    :type tolerance: float, optional
    :return: x_r such that x_r approx x_bar where f(x_bar)=0
    :rtype: float
    '''

    if not ( f(lower) * f(upper) < 0): # check the initial interval contains a root
        raise ValueError("bounds have same sign")   
    middle = (lower + upper)/2
    while 0.5*abs(upper-lower) > tolerance:  # check that we're not converged
        middle = (lower + upper)/2                  # midpoint of current interval
        if f(lower) * f(middle) < 0:           # select which 1/2 interval to continue with
            upper = middle
        else:
            lower = middle
    return middle

def root_finder_newton(f:Callable[[float],float], df:Callable[[float],float],x0:float, max_iter:int = 100, precision:float=1e-10)->float:
    '''Newton-Raphson root finding algorithm in 1D

    :param f: function to find root of
    :type f: Callable[[float],float]
    :param df: derivative of said function
    :type df: Callable[[float],float]
    :param x0: initial guess 
    :type x0: float
    :param max_iter: maximum allowed iterations, defaults to 100
    :type max_iter: int, optional
    :param precision: value of f(x) below which method exists, defaults to 1e-6
    :type precision: float, optional

    :return: x s.t. f(x) approx 0
    :rtype: float
    '''

    for i in range(max_iter):
        fx = f(x0)
        if abs(fx) < precision: return x0
        dx = fx/df(x0)
        x0 = x0 - dx
    else: raise ArithmeticError("Newton's method failed to converge")

def root_finder_fallback(f:Callable[[float],float],
                         df:Callable[[float],float],
                         lower:float, 
                         upper:float,
                         max_iter:int = 100, 
                         precision:float=1e-10):
    '''Root finder that tries to use newton, but falls back on bisection
    if newton fails to converge'''
    try:
        x0 = root_finder_newton(f,df,(upper+lower)/2,max_iter,precision)
        if not (lower < x0 < upper): raise ArithmeticError("Converged wrong")
    except ArithmeticError:
        x0 = root_finder_bisection(f,lower,upper,precision)
    return x0

def minimizer_old(f:Callable[[np.ndarray],float],
              x0:np.ndarray, 
              initial_step:np.ndarray,
              precision:float = 1e-10, 
              max_iter:int=1000, 
              allow_nonconvergence:bool=False)->np.ndarray:
    '''function to find a minimum value over a R^n->R function.
    currently implements Nelder-Mead as the algorithm of choice'''

    a = 1 # default nelder mead coefficients
    b = 0.5
    c = 2
    d = 0.5

    # currently only 2d is implemented TODO
    if x0.size != 2: raise NotImplementedError('minimizer only works for 2d domains')
    
    p1 = [0,x0] # initial points
    p2 = [0,x0 - np.array([initial_step[0],0])]
    p3 = [0, x0 - np.array([0,initial_step[1]])]
    p1[0] = f(p1[1])
    p2[0] = f(p2[1])
    p3[0] = f(p3[1])
    avg_point = lambda p1,p2,p3: ((p1[1][0] + p2[1][0] + p3[1][0])/3, (p1[1][1] + p2[1][1] + p3[1][1])/3)

    for _ in range(max_iter):
            
        if p1[0] > p2[0]: p1,p2 = p2,p1 # ordering (lowest first)
        if p2[0] > p3[0]: p2,p3 = p3,p2
        if p1[0] > p2[0]: p1,p2 = p2,p1

        
        m = (p1[0] + p2[0] + p3[0])/3
        var = (((p1[0]-m)**2 + (p2[0]-m)**2 + (p3[0]-m)**2)/2)
        if var < precision**2: # termination (based on ssd of function values)
            return np.array(avg_point(p1,p2,p3))

        cent = 0.5*(p1[1] + p2[1]) # centroid

        reflect_p = cent + a*(cent - p3[1]) # transform
        reflect = [f(reflect_p), reflect_p]
        if p1[0] <= reflect[0] < p2[0]:
            p3 = reflect
            continue

        elif reflect[0] < p1[0]: # expand
            expand_p = cent + c*(reflect_p - cent)
            expand = [f(expand_p), expand_p]
            if expand[0] < reflect[0]:
                p3 = expand
                continue
            else:
                p3 = reflect
                continue
        elif reflect[0] < p2[0]: # contract
            if reflect[0] < p3[0]: 
                contract_p = cent + b*(reflect_p - cent)
            else:
                contract_p = cent + b*(p3[1] - cent)
            contract = [f(contract_p), contract_p]
            if contract[0] < p3[0]:
                p3 = contract
                continue
        else: # shrink
            p3_p = p1[1] + d*(p3[1]-p1[1])
            p3 = [f(p3_p), p3_p]
            p2_p = p1[1] + d*(p2[1]-p1[1])
            p2 = [f(p2_p), p2_p]
            continue
    else:
        if allow_nonconvergence: return np.array(avg_point(p1,p2,p3))
        else: raise ArithmeticError("Nelder-mead failed to converge")
    
    
    raise NotImplementedError("optimizer not implemented, sorry :(")
    

def minimizer(f:Callable[[np.ndarray],float],
              x0:np.ndarray, 
              initial_step:np.ndarray,
              precision:float = 1e-10, 
              max_iter:int=1000, 
              allow_nonconvergence:bool=False)->np.ndarray:
    '''function to find a minimum value over a R^n->R function.
    currently implements Nelder-Mead as the algorithm of choice'''


    a = 1 # default nelder mead coefficients
    b = 0.5
    c = 2
    d = 0.5

    parr = [x0]
    for i in range(len(x0)): # square parr
        x1 = x0.copy()
        x1[i] = x1[i] +  initial_step[i]
        parr.append(x1)
    parr = np.array(parr)

    Farr = []
    for p in parr:
        Farr.append(f(p))
    Farr = np.array(Farr)
    
    for _ in range(max_iter):

        idx = Farr.argsort() # sort ascending
        Farr = Farr[idx]; parr = parr[idx]

        std = np.std(Farr)
        if std < precision: # termination (based on ssd of function values)
            return np.average(parr, axis=0)

        po = np.average(parr[0:-1], axis=0) # centroid of all except last

        pr = po + a*(po - parr[-1]) # reflect
        Fr = f(pr)
        if Farr[0] <= Fr < Farr[-2]:
            parr[-1] = pr
            Farr[-1] = Fr
            continue

        if Fr < Farr[0]: # expand
            pe = po + c*(pr - po)
            Fe = f(pe)
            if Fe < Fr:
                parr[-1] = pe
                Farr[-1] = Fe
                continue
            else:
                parr[-1] = pr
                Farr[-1] = Fr
                continue
        if Fr < Farr[-1]: # contract setup
            pc = po + b*(pr - po)
            Fc = f(pc)
            Fcomp = Fr
        else:
            pc = po + b*(parr[-1] - po)
            Fc = f(pc)
            Fcomp = Farr[-1]
        if Fc < Fcomp: # contract comparison
            parr[-1] = pc
            Farr[-1] = Fc
        else: # shrink
            for i in range(1,len(parr)):
                parr[i] = parr[0] + d*(parr[i] - parr[0])
                Farr[i] = f(parr[i])
            continue
    else:
        if allow_nonconvergence: return np.average(parr, axis=0)
        else: raise ArithmeticError("Nelder-mead failed to converge")