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
from ..utils import r2f
__all__ = [
    'plot',
    'center',
    'add_solar_system',
    'clf',
    'show',
    'set_view_angle',
]



ax:Axes = None # type:ignore


def init():
    global ax
    if not ax is None: return;
    else: return clf()
    

def clf():
    '''clear plot so a new plot can be drawn'''
    global ax
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'},
                            gridspec_kw=dict(top=1.07, bottom=0.02, left=0, right=1))
    ax.grid(False)
    ax.set_axis_off()
    mpl.rcParams['axes3d.mouserotationstyle'] = 'azel'

    # default view:
    ax.azim = 0 # type:ignore
    ax.elev = 90 # type:ignore

    return;

def set_view_angle(az:float,el:float,r:float=10):
    '''viewing angle of next plot, angles in degrees, r seems not to work
    due to pyplot'''
    init()
    global ax
    ax.azim = az # type:ignore
    ax.elev = el # type:ignore
    ax.dist = r # type:ignore
    return;

# to make:
def plot(
        ob:Orbit,
        f_bounds:tuple[float,float]|None = None,
        t_bounds:tuple[float,float]|None = None,
        *,
        stilts:bool|float = m.radians(5), # min inclination to add stilts
        stilt_number:int = 20, # number of stilts
        stilt_spacing:str = 'angle',
        f:float|None = None,
        t:float|None = None,
        forward_predict:bool = False,
        backward_predict:bool = False,
        max_distance:float = m.inf,
        point_label:str|None = None,
        resolution:int = 360,
        **kwargs
)->None:
    '''plot an orbit in the 3d plane, along with a point if supplied.
    bounds in time of true anomaly limit what is shown,
    stilts are added if true or if inclination is above given value.
    forward and backward predict plots a dashed line for the parts of the orbit 
    not inside the bounds
    '''
    if not t_bounds is None:
        try:
            ft_low = ob.f(t_bounds[0]) if m.isfinite(t_bounds[0]) else -2*m.pi
        except (OverflowError):
            ft_low = -2*m.pi
        try:
            ft_high = ob.f(t_bounds[1]) if m.isfinite(t_bounds[1]) else 2*m.pi
        except (OverflowError):
            ft_high = 2*m.pi
        while ft_low > ft_high: ft_high += 2*m.pi

        if t_bounds[1]-t_bounds[0] > ob.T: ft_low = -2*m.pi; ft_high = 2*m.pi

    if f_bounds is None and t_bounds is None:
        # full circle of -finf to finf
        if ob.e > 1:  f_low = -ob.finf; f_high = ob.finf
        else: f_low = -m.pi; f_high = m.pi
    elif f_bounds is None: # t_bounds given
        if ob.e > 1:
            f_low = max(-ob.finf, ft_low) # type:ignore
            f_high = min(ob.finf, ft_high) # type:ignore
        else: 
            f_low = ft_low; f_high = ft_high # type: ignore
    elif t_bounds is None: # fbounds given
        if ob.e > 1:
            f_low = max(-ob.finf, f_bounds[0]) # type:ignore
            f_high = min(ob.finf, f_bounds[1]) # type:ignore
        else: f_low = f_bounds[0]; f_high = f_bounds[1] # type: ignore
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
    pp = ob.point_locus(f_low, f_high, resolution).T

    # get rid of things beyond max distance:
    pp = pp[:,np.linalg.norm(pp,axis=0) < max_distance]

    # generate colors if not supplied
    kwargs.setdefault('color',np.random.random(3))
    point_size = kwargs.pop('s', 40)
    if point_size is None: point_size = 40

    # plot the orbit:
    init()
    ax.plot(pp[0],pp[1],pp[2], **kwargs)

    kwargs.pop('label',None) # to not duplicate labels

    if (not f is None) or (not t is None): # plot point
        
        p = ob.t2rvec(t) if t is not None else ob.rvec(f) # type: ignore
        if np.linalg.norm(p) < max_distance:
            ax.scatter(p[0],p[1],p[2], **kwargs, s = float(point_size)) # type:ignore

            # add label:
            if point_label:
                ax.text(p[0],p[1],p[2], point_label, ha="center", va="center", weight='bold', zorder=99) # type:ignore
    
    # add stilts:
    if stilts is True or (isinstance(stilts,float) and abs(ob.i) > stilts):
        match stilt_spacing:
            case "range": # based on distance from barycenter
                r = np.linspace(ob.pe, min(max_distance, ob.ap), stilt_number)
                ss = r2f(r, ob.p, ob.e)
                ss = np.concat((ss,[0], -ss))
            # case "distance": # based on distance travelled. TODO

            case _: #angle
                ss =  np.linspace(max(-np.pi, -ob.finf), min(np.pi, ob.finf), stilt_number) # straight angle arange
        ss = ss[f_low <= ss]; ss = ss[ss <= f_high]
        ss = ob.rvec(ss)
        ss = ss[np.linalg.norm(ss,axis=1) <= max_distance]
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
    kwargs.setdefault('s', 70)
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
             s = (90 if (symbols or initials) else None)
             )
    center(color='gold', s = (90 if (symbols or initials) else 50))

def show(**kwargs):
    '''wrapper of plt.show that ensures the frame is cleared for future renders'''
    plt.show(**kwargs)

    global ax # erase for next show
    ax = None # type:ignore
