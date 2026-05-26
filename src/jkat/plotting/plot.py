''' 
the pyplot-like plotter
'''
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib as mpl
import math as m
import numpy as np

from ..kep import Orbit
from ..ephemeris import *
__all__ = [
    'plot',
    'center',
    'add_solar_system'
]



ax:Axes = None # type:ignore


def init():
    global ax
    if not ax is None: return;
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'},
                            gridspec_kw=dict(top=1.07, bottom=0.02, left=0, right=1))
    ax.grid(False)
    ax.set_axis_off()
    mpl.rcParams['axes3d.mouserotationstyle'] = 'azel'



# to make:

def plot(
        ob:Orbit,
        f_bounds:tuple[float,float]|None = None,
        t_bounds:tuple[float,float]|None = None,
        *,
        stilts:bool|float = m.radians(5),
        stilt_separation:float = m.radians(10),
        f:float|None = None,
        t:float|None = None,
        forward_predict:bool = False,
        backward_predict:bool = False,
        max_distance:float = m.inf,
        point_label:str|None = None,
        **kwargs
)->None:
    '''plot an orbit in the 3d plane, along with a point if supplied.
    bounds in time of true anomaly limit what is shown,
    stilts are added if true or if inclination is above given value.
    forward and backward predict plots a dashed line for the parts of the orbit 
    not inside the bounds
    '''
    if not t_bounds is None:
        ft_low = ob.f(t_bounds[0]) if m.isfinite(t_bounds[0]) else -2*m.pi
        ft_high = ob.f(t_bounds[1]) if m.isfinite(t_bounds[1]) else 2*m.pi
        while ft_low > ft_high: ft_high += 2*m.pi

    if f_bounds is None and t_bounds is None:
        # full circle of -finf to finf
        if ob.e > 1:  f_low = -ob.finf; f_high = ob.finf
        else: f_low = -m.pi; f_high = m.pi
    elif f_bounds is None: # t_bounds given
        if ob.e > 1:
            f_low = max(-ob.finf, ft_low) # type:ignore
            f_high = min(ob.finf, ft_high) # type:ignore
        else: 
            f_low = ob.f(t_bounds[0]); f_high = ft_high # type: ignore
    elif t_bounds is None: # fbounds given
        if ob.e > 1:
            f_low = max(-ob.finf, f_bounds[0]) # type:ignore
            f_high = min(ob.finf, t_bounds[1]) # type:ignore
        else: f_low = f_bounds[0]; f_high = t_bounds[1] # type: ignore
    else:
        if ob.e > 1:
            f_low = max(-ob.finf, ft_low, f_bounds[0]) # type:ignore
            f_high = min(ob.finf, ft_high, f_bounds[1]) # type:ignore
        else: 
            f_low = max(ft_low, f_bounds[1])
            f_high = min(ft_high, f_bounds[1]) # type: ignore

    # ensure no double wrap:
    if f_high - f_low > 2*m.pi: f_high = f_low + 2*m.pi

    # get locus:
    pp = ob.point_locus(f_low, f_high).T

    # get rid of things beyond max distance:
    pp = pp[:,np.linalg.norm(pp,axis=0) < max_distance]

    # generate colors if not supplied
    kwargs.setdefault('color',np.random.random(3))
    pw = kwargs.pop('pw',3)

    # plot the orbit:
    init()
    ax.plot(pp[0],pp[1],pp[2], **kwargs)

    kwargs.pop('label',None) # to not duplicate labels

    if (not f is None) or (not t is None): # plot point
        
        p = ob.t2rvec(t) if t is not None else ob.rvec(f) # type: ignore
        ax.scatter(p[0],p[1],p[2], **kwargs, lw=pw)

        # add label:
        if point_label:
            ax.text(p[0],p[1],p[2], point_label, ha="center", va="center", weight='bold', zorder=99) # type:ignore
    
    # add stilts:
    if stilts is True or (isinstance(stilts,float) and ob.i > stilts):
        ss =  np.arange(f_low, f_high+1e-6, stilt_separation)
        ss = ob.rvec(ss)
        ss = ss[np.linalg.norm(ss,axis=1) < max_distance]
        for s in ss:
            ax.plot((s[0],s[0]),
                    (s[1],s[1]),
                    (s[2],0), **kwargs)
        # ax.stem(ss[0],ss[1],ss[2], linefmt='--', markerfmt='', **kwargs)
    
    # add predictions:
    if forward_predict or backward_predict: raise NotImplementedError()


    

    ax.set_aspect('equal') # to keep scale aspect 1:1:1
    return;



def center(**kwargs):
    '''add center to the plot'''
    init()
    kwargs.setdefault('lw', 14)
    kwargs.setdefault('color', 'gold')
    ax.scatter(0,0,0, **kwargs)
    return;

def add_solar_system(t:float = 0, planets:str='11111000',symbols:bool=False, initials:bool=False, **kwargs):
    '''adds the solar system. planet string determines which planets to add'''

    pla = [Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
    colors = ['silver','wheat','cornflowerblue','orangered','darkorange','khaki','dodgerblue','royalblue']
    sym = ['☿', '♀', '♁', '♂', '♃', '♄', '⛢', '♆'] # 🜨, ♁
    initial = ['H','V','E','M','J','S','U','N']

    for i in range(len(planets)):
        if planets[i] == '0': continue
        plot(pla[i], t=t, color = colors[i], **kwargs, 
             point_label= (sym[i] if symbols else (initial[i] if initials else None)),
             pw = (9 if (symbols or initials) else None)
             )
    center(color='gold', lw = (9 if (symbols or initials) else 4))

