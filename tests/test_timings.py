''' 
test all the time conversions

'''

from src.utils.anomalies import *
from src.utils import EARTH_RADIUS, EARTH_MU, a2T, h2a, finf
import math as m
import numpy as np

import pytest
from pytest import approx

REL = 2e-3

def test_kepler_universal_time_elliptical():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    for e in [0.0,0.3,0.5,0.7,0.9]:
        period = a2T(h2a(h,e,mu),mu)
        tt = np.linspace(0,period,endpoint=False)
        for t in tt:
            kep = t2f_kep(t,e,h,mu)
            uni = t2f(t,e,h,mu)
            assert uni == approx(kep, rel=REL)

def test_kepler_universal_time_hyper():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    period = 90*60
    for e in [1.0,1.5,2.0,3.0,4.0]:
        tt = np.linspace(-period,period,endpoint=False)
        for t in tt:
            try:
                kep = t2f_kep(t,e,h,mu)
                uni = t2f(t,e,h,mu)
                assert uni == approx(kep, rel=REL)
            except ArithmeticError: continue

def test_kepler_universal_theta_elliptical():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    for e in [0.0,0.3,0.5,0.7,0.9]:
        ff = np.linspace(0.01,2*m.pi,endpoint=False)
        for f in ff:
            kep = f2t_kep(f,e,h,mu)
            uni = f2t(f,e,h,mu)
            assert uni == approx(kep, rel=REL)

def test_kepler_universal_theta_hyper():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    for e in [1.0,1.5,2.0,3.0,4.0]:
        asymp_ang = finf(e)
        if m.isnan(asymp_ang): asymp_ang = m.pi
        ff = np.linspace(-asymp_ang + 0.05,asymp_ang - 0.05,endpoint=False)
        for f in ff:
            try:
                kep = f2t_kep(f,e,h,mu)
                uni = f2t(f,e,h,mu)
                assert uni == approx(kep, rel=REL)
            except ArithmeticError: continue


def test_uni_round_trip_elliptical():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    for e in [0.0,0.3,0.5,0.7,0.9]:
        ff = np.linspace(0,2*m.pi, 49,endpoint=False)
        for f in ff:
            t = f2t(f,e,h,mu)
            f2 = t2f(t,e,h,mu)
            t2 = f2t(f2,e,h,mu)
            assert f == approx(f2, rel=REL)
            assert t == approx(t2, rel=REL)

def test_uni_round_trip_hyper():
    h = EARTH_RADIUS*8  
    mu = EARTH_MU
    for e in [1.0,1.5,2.0,3.0,4.0]:
        asymp_ang = finf(e)
        if m.isnan(asymp_ang): asymp_ang = m.pi
        ff = np.linspace(-asymp_ang + 0.05,asymp_ang - 0.05,endpoint=False)
        for f in ff:
            try:
                t = f2t(f,e,h,mu)
                f2 = t2f(t,e,h,mu)
                t2 = f2t(f2,e,h,mu)
                assert f == approx(f2, rel=REL)
                assert t == approx(t2, rel=REL)
            except ArithmeticError: continue