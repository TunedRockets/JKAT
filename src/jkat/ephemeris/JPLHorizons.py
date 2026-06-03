'''
sub-module for interfacing with the JPL hoizons system of Ephemeris

'''
import numpy as np
import requests as req
import re
import datetime as dt
import math
from pathlib import Path
from typing import Callable
import time

from ..utils.consts import JULIAN_YEAR
from .realtime import to_JD, from_JD, to_time
from ..utils import a2p
from ..kep import Orbit
import pickle
__all__ = [
    'Ephemerical_Orbit',
    'horizons_request'
]


URL = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text"


# =====

URI_INVALID = r'[^a-zA-Z0-9%._~-]|(%(?:(?![0-9]{2})))'
def URI_sanitize(link:str):
    '''sanitize a link by replacing invalid characters with % codes'''

    while not (m := re.match(URI_INVALID, link)) is None:
        target = m[0]
        sub = ''.join([(f'%{b:X}') for b in target.encode()])
        link = re.sub(m[0], sub, link)
    return link

class Ephemerical_Orbit(Orbit):
    '''Modification of the orbit to take into account the precession of 
    the orbital elements.
    for access, set current_time to wanted time.
    will be automatically set for any function that uses time'''

    def __init__(self,
                 pp: np.ndarray,
                 ee: np.ndarray,
                 ii: np.ndarray,
                 rraan: np.ndarray,
                 aargp: np.ndarray,
                 ttp: np.ndarray,
                 tt:np.ndarray,
                 mu: float,
                 name:str,
                 jplid:int) -> None:
        
        self.pp = pp; self.ee = ee; self.ii = ii; self.rraan = rraan; self.tt = tt
        self.aargp = aargp; self.ttp = ttp
        self.current_time = 0
        self.mu = mu
        self.name = name
        self.jplid = jplid
    @property
    def p(self): return np.interp(self.current_time,self.tt,self.pp)
    @property
    def e(self): return np.interp(self.current_time,self.tt,self.ee)
    @property
    def i(self): return np.interp(self.current_time,self.tt,self.ii)
    @property
    def raan(self): return np.interp(self.current_time,self.tt,self.rraan)
    @property
    def argp(self): return np.interp(self.current_time,self.tt,self.aargp)
    @property
    def tp(self): return np.interp(self.current_time,self.tt,self.ttp)

    def t2rvec(self, t: float) -> np.ndarray:
        self.current_time = t
        return super().t2rvec(t)
    def t2vvec(self, t: float) -> np.ndarray:
        self.current_time = t
        return super().t2vvec(t)
    def t2vectors(self, t: float) -> tuple[np.ndarray, np.ndarray]:
        self.current_time = t
        return super().t2vectors(t)
    
    def f(self, t: float = 0, after_periapsis: bool = False) -> float:
        self.current_time = t
        return super().f(t, after_periapsis)
    def link_tf(self, t: float, f: float) -> float:
        self.current_time = t
        return super().link_tf(t, f)
    
    def t2M(self, t: float, after_periapsis: bool = False) -> float:
        self.current_time = t
        return super().t2M(t, after_periapsis)
    def t2L(self, t: float, after_periapsis: bool = False) -> float:
        self.current_time = t
        return super().t2L(t, after_periapsis)
    def t2X(self, t: float, after_periapsis: bool = False) -> float:
        self.current_time = t
        return super().t2X(t, after_periapsis)
    
    def osculating_orbit(self,t):
        self.current_time = t
        return Orbit(self.p,self.e,self.i,self.raan,self.argp,self.tp,self.mu)
    
    def __repr__(self) -> str:
        return f"{self.name} ({self.jplid}):\n{self.current_time=}\n {self.p=}\n {self.e=}\n {self.i=}\n {self.raan=}\n {self.argp=}\n {self.tp=}\n {self.mu=}"
    

# ====
CACHE = Path(__file__).parent / 'JPLHorizons.cache'
def _request_cache(request:str, grab_first:bool)->Ephemerical_Orbit:
    '''hold a cache of requests to not need to contact jpl all the time'''

    try:
        with open(CACHE, 'rb') as file:
            cache = pickle.load(file)
    except FileNotFoundError: cache = {}

    result = cache.get(request, None)
    if result is None:
        result = _request_to_orbit(request, grab_first)
        cache[request] = result
        with open(CACHE, 'wb') as file:
            pickle.dump(cache,file)
    return result





def horizons_request(name:str|int,
                     center:str|int=10,
                     ecliptic:bool=False,
                     t_start:float|dt.datetime = 0,
                     t_end:float|dt.datetime = JULIAN_YEAR,
                     steps:int = 50,
                     grab_first:bool=True)->Ephemerical_Orbit:
    '''request ephemeris from horizons and get mu and the osculating elements as a list along with time.
    if several items match the name, it will pick the first that matches.
    Stores all requests and their returns in a cache'''

    if isinstance(name,int): name = str(name)
    if isinstance(center,int): center = str(center)
    if center.lower() in ('0', '10', 'sun', 'ssb'): ecliptic = True # nothing else for the sun

    # create request
    request = URL

    request += f"&COMMAND='{URI_sanitize(name)}'"
    request += f"&CENTER='{URI_sanitize(center)}'"
    ecl = 'ECLIPTIC' if ecliptic else 'BODY_EQUATOR'
    request += f"&REF_PLANE='{ecl}'"
    request += f"&EPHEM_TYPE='ELEMENTS'"
    request += f"&OBJ_DATA='{"NO"}'" # don't need that/not implemented yet
    request += f"&STEP_SIZE='{steps}'"

    # times (in JD)
    if isinstance(t_start, dt.datetime): t_start = to_time(t_start)
    if isinstance(t_end, dt.datetime): t_end = to_time(t_end)
    str_start = f'JD {to_JD(t_start, MJD2000=False)}'
    str_end = f'JD {to_JD(t_end, MJD2000=False)}'
    request += f"&START_TIME='{str_start}'"
    request += f"&STOP_TIME='{str_end}'"

    ob = _request_cache(request, grab_first=grab_first)
    return ob

def _request_to_orbit(request, grab_first:bool)->Ephemerical_Orbit:
    '''internal function incorporating a cache'''

    re.M = True
    result = req.get(request)
    if not result.ok: raise LookupError(f'Bad return from JPL-Horizons: {result.text}')
    result = result.text

    if re.search(r'(Multiple major-bodies match string)|(Small-body Index Search Results)', result):
        # found many:
        if not grab_first: raise LookupError("Several items found with given name")
        name = re.search(r"COMMAND='([^']*)'",request)[1] # type:ignore

        # grab all occurances:
        tgt_id = None
        matches = re.finditer(r'(\n|^)[\s]*([-0-9]+)[\s]*(([\w]+ ?)+)', result)
        for m in matches:
            
            jplname = m[3]
            jplid = m[2]
            
            if name.strip().lower() == jplname.strip().lower():
                # perfect match, use this:
                tgt_id = jplid
                break
            if tgt_id == None:
                tgt_id = jplid # else use first
        # replace search:
        request = re.sub(r"COMMAND='([^']*)'",f"COMMAND='{tgt_id}'", request)
        return _request_to_orbit(request, False)


    mu = re.search(r'Keplerian GM[\s]*:[\s]*([\-+0-9E\.]+)', result)
    if mu is None: raise LookupError("Couldn't find value: Keplerian GM")
    mu = float(mu[1])

    # get names:
    m = re.search(r'Target body name:[\s]*([0-9]*)?[\s]*(.*)[\s]+\((.*)\)[\s]*{source:', result)
    if m is None: raise LookupError("Couldn't find value: Target body name")
    # why JPL do you make the name scheme inconsistent?!?!?!
    if m[1] != '': # small body
        jplidx = int(m[1])
    else: # main body
        jplidx = int(m[3])
    jplname = m[2]

    # get ephemeris values
    content = re.search(r'\$\$SOE(?s:(.*))\$\$EOE', result)
    if content is None: raise LookupError("Couldn't find value: ephemeris table")
    content = content[1]

    # extract times:
    times = []
    while not (m := re.search(r'([0-9.]+) = A.D.', content)) is None:
        content = content.replace(m[0], '')
        t = float(m[1])
        times.append(from_JD(t))

    # extract eccentricity:
    ee = []
    while not (m := re.search(r'EC= ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        e = float(m[1])
        ee.append(e)
    
    # extract inclination:
    ii = []
    while not (m := re.search(r'IN= ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        i = float(m[1])
        ii.append(math.radians(i))
    
    # extract argument of periapsis:
    ww = []
    while not (m := re.search(r'W = ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        w = float(m[1])
        ww.append(math.radians(w))
    
    # extract raan:
    oo = []
    while not (m := re.search(r'OM= ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        o = float(m[1])
        oo.append(math.radians(o))
    
    # extract parameter (by way of semi-major axis):
    pp = []
    idx = 0
    while not (m := re.search(r'A = ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        a = float(m[1])
        pp.append(a2p(a, ee[idx]))
        idx += 1

    # extract tp:
    ttp = []
    while not (m := re.search(r'Tp=  ([0-9.E\-+]+)', content)) is None:
        content = content.replace(m[0], '')
        tp = float(m[1])
        ttp.append(from_JD(tp))

    if not (len(times) == len(ee) == len(ii) == len(ww) == len(oo) == len(pp) == len(ttp)):
        raise LookupError('Malformed result, (element length mismatch)')
    pp = np.array(pp); ee = np.array(ee); ii = np.array(ii)
    oo = np.array(oo); ww = np.array(ww); ttp = np.array(ttp)
    times = np.array(times)  


    # create orbit:
    ob = Ephemerical_Orbit(
        pp,ee,ii,oo,ww,ttp, times,mu,
        jplname,jplidx
    )
    ob.current_time = ob.tt[0]
    return ob