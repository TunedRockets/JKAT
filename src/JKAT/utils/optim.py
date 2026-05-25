''' 
mathematical optimizers,

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
    'bb_glob_optim'
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

def bb_glob_optim(f:Callable[[np.ndarray],float], bounds:Sequence[tuple[float,float]], **kwargs):
    '''wrapper for global optimization black box algorithm, returns list of local optima'''
    raise NotImplementedError("optimizer not implemented, sorry :(")
    x = o.shgo(f,bounds,**kwargs, workers=-1)
    return x.xl