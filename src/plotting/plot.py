''' 
the pyplot-like plotter
'''
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import math as m
import numpy as np

from ..kep import Orbit
__all__ = [
    'plot',
    'center'
]



ax:Axes = None # type:ignore


def init():
    global ax
    if not ax is None: return;
    ax = plt.figure().add_subplot(projection='3d')
    


# to make:

def plot(
        ob:Orbit,
        f_bounds:tuple[float,float] = (0,2*m.pi),
        t_bounds:tuple[float,float] = (-m.inf,m.inf),
        *,
        stilts:bool = True,
        stilt_separation:float = m.radians(30),
        f:float|None = None,
        t:float|None = None,
        forward_predict:bool = False,
        backward_predict:bool = False,
        **kwargs
)->None:
    '''plot an orbit in the 3d plane, along with a point if supplied
    '''

    # set limits:
    if ob.e > 1: f_bounds = (
        max(f_bounds[0], -ob.finf),
        min(f_bounds[1], ob.finf)
    )
    if t_bounds[0] != -m.inf:
        f_low = max(f_bounds[0], ob.f(t_bounds[0]))
    else: f_low = f_bounds[0]
    if t_bounds[1] != m.inf:
        f_high = max(f_bounds[1], ob.f(t_bounds[1]))
    else: f_high = f_bounds[1]

    # get locus:
    pp = ob.point_locus(f_low, f_high).T

    # generate colors if not supplied
    kwargs.setdefault('color',np.random.random(3))

    # plot the orbit:
    init()
    ax.plot(pp[0],pp[1],pp[2], **kwargs)

    kwargs.pop('label',None) # to not duplicate labels

    if (not f is None) or (not t is None): # plot point
        p = ob.t2rvec(t) if t is not None else ob.rvec(f) # type: ignore
        ax.scatter(p[0],p[1],p[2], **kwargs)
    
    # add stilts:
    if stilts:
        ss =  np.arange(f_low, f_high+1e-6, stilt_separation)
        ss = ob.rvec(ss)
        for s in ss:
            ax.plot((s[0],s[0]),
                    (s[1],s[1]),
                    (s[2],0), **kwargs)
        # ax.stem(ss[0],ss[1],ss[2], linefmt='--', markerfmt='', **kwargs)
    
    # add predictions:
    if forward_predict or backward_predict: raise NotImplementedError()
    return;



def center(**kwargs):
    '''add center to the plot'''
    init()
    kwargs.setdefault('lw', 3)
    kwargs.setdefault('color', 'red')
    ax.scatter(0,0,0, **kwargs)
    return;


# =============================

def plot_orbit(ax,ob:Orbit,time:float|tuple[float,float]=0.0,trail:float=2*m.pi, ThreeDee:bool=True,hyper_predict:bool=False, max_alt:float=m.inf, **kwargs)->None:
    '''Plot orbit in the given axis, trail determines how far behind the orbit is plotted (defaults to entire orbit)
    hyperbolic orbits will only be plotted up to the current point, if hyper_predict is true, a dashed line will be plotted ahead,
    max_size determines how far out to plot,
    if time is a tuple, trail is ignored and instead using the two times'''

    if isinstance(time, float):
        theta = ob.time_to_theta(time)

        if trail != 2*m.pi:
            # add trail:
            end_theta = theta
            start_theta = theta-trail
        else: # no trail, 
            end_theta = m.pi
            start_theta = -m.pi
        if ob.e >= 1: end_theta = theta
    else:
        start_theta = ob.time_to_theta(time[0]) #type:ignore
        end_theta = ob.time_to_theta(time[1]) # type:ignore
        theta = end_theta


    cross = ob.crosses_altitude(max_alt)
    if cross is None and ob.periapsis > max_alt: return # don't render anything
    elif not cross is None:
        start_theta = bounds(-cross,start_theta,cross)
        end_theta = bounds(-cross,end_theta,cross)


    locus = ob.point_locus(start_theta,end_theta)
    if cross is None or abs(theta) < cross:
        point = ob.theta_to_rv(theta)[0]
    else: point = np.array([m.inf])

    if not 'color' in kwargs:
        kwargs['color'] = np.random.random(3)
    if ThreeDee: ax.plot(locus[:,0],locus[:,1],locus[:,2], **kwargs)
    else: ax.plot(locus[:,0],locus[:,1], **kwargs)
        
    kwargs.pop('label',None) # to not duplicate labels
    if np.linalg.norm(point) <= max_alt:
        if ThreeDee: ax.scatter(point[0],point[1],point[2], **kwargs)
        else: ax.scatter(point[0],point[1], **kwargs)
    
    if hyper_predict and ob.e >= 1:
        start_theta = end_theta
        end_theta = cross
        locus = ob.point_locus(start_theta,end_theta) # type:ignore
        if ThreeDee: ax.plot(locus[:,0],locus[:,1],locus[:,2], ls='--', **kwargs)
        else: ax.plot(locus[:,0],locus[:,1], **kwargs)
    return;

def get_solar_system_ax():
    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(0,0,0, lw=3, color="red")
    return ax