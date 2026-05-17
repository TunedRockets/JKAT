'''
Tests that run the examples in curtis' book.
'''

import pytest
from pytest import approx
import math as m
import numpy as np

from src.kep import Orbit, orbit_from_keplerian, orbit_from_rv, orbit_from_lambert
from src.utils import EARTH_RADIUS, EARTH_MU, DAY, elazr2vec




@pytest.fixture
def fixt_ob()->Orbit:
    '''standard earth orbit'''
    return orbit_from_keplerian(EARTH_RADIUS,0,0,0,0,0,EARTH_MU)


def test_curtis_2_5(fixt_ob:Orbit):
    ob = fixt_ob
    ob.T = DAY
    assert ob.apoapsis == approx(EARTH_RADIUS + 35_786)
    assert ob.vt(0) == approx(3.075)

def test_curtis_2_7( fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(4000+6378,400+6378)
    assert ob.e == approx(0.2098)
    assert ob.h == approx( 57_172)
    assert ob.a == approx(8578)
    assert ob.period == approx(2.196*(60*60))

def test_curtis_2_9(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 1
    ob.pe = 7000
    theta1 = ob.cross_radius(8000)
    assert m.isfinite(theta1)
    theta2 = ob.cross_radius(16000)
    assert  m.isfinite(theta2)
    p1 = ob.rvec(theta1)
    p2 = ob.rvec(theta2)
    d = np.linalg.norm(p1-p2)
    assert d == approx(13_270)


def test_curtis_2_10():
    
    r = np.array([14_600,0,0])
    v = elazr2vec(0,m.radians(90-50), 8.6)
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(80_708)
    assert ob.e == approx(1.3393)
    assert ob.periapsis == approx(6986)
    assert ob.a == approx(-20_590)
    assert ob.c3 == approx(19.36)
    assert ob.aiming_radius == approx(18_340)



def test_curtis_2_11(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 0.3
    ob.h = 60_000
    r,v = ob.vectors(np.radians(120))
    assert r == approx(np.array([-5312.7, 9201.9, 0]))
    assert v == approx(np.array([-5.7533, -1.3287, 0]))

def test_curtis_2_12():
    r = np.array([7000,9000,0])
    v = np.array([-3.3472,9.1251,0])
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(94_000)
    assert ob.f(0) == approx(m.radians(52.125))
    assert ob.e == approx(1.538)

def test_curtis_2_13_and_14():
    r = np.array([8182.4,-6865.9,0])
    v = np.array([0.47572,8.8116,0])
    ob = orbit_from_rv(r,v,EARTH_MU)
    t = ob.f(0)
    t2 = t + m.radians(120)
    r,v = ob.t2vectors(t2)
    r_facit = np.array([1454.9, 8251.6,0])
    v_facit = np.array([-8.1323,5.6785, 0])
    assert r == approx(r_facit)
    assert v == approx(v_facit)
    assert ob.e == approx(1.0563)
    assert t % (2*m.pi) == approx(m.radians(288.44))


def test_curtis_3_1(fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(9600,21000)
    T = ob.t(np.radians(120))
    assert T == approx(4077)

def test_curtis_3_2(fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(9600,21000)
    assert ob.e == approx(0.37255)
    theta = ob.f(3*60*60)
    assert theta == approx(np.radians(193.2))

def test_curtis_3_4(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 1
    # not testing the energy eq: TODO change
    ob.pe = 7972
    ob.tp = 0
    f = ob.f(6*60*60)
    r = ob.r(f)
    assert r == approx(86_899)

def test_curtis_3_5():
    ob = orbit_from_rv(np.array([6378+300,0,0]), np.array([0,15,0]),EARTH_MU)
    T = ob.t(np.radians(100))
    r = ob.r(np.radians(100))
    assert T == approx(4141.4)
    assert r == approx(48_497)
    r,v = ob.t2vectors(T + 3*60*60)
    assert np.linalg.norm(r) == approx(163_180)
    assert np.linalg.norm(v) == approx(10.51)

def test_curtis_3_6(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 1.4682
    ob.h = 95_154
    ob.link_tf(0,m.radians(30))
    assert ob.r(m.radians(30)) == approx(10_000)
    assert ob.v(m.radians(30)) == approx(10)
    t = ob.f(1*60*60)
    assert t % (2*m.pi) == approx(m.radians(100.04))

def test_curtis_3_7():
    ob = orbit_from_rv(np.array([7000,-12_124, 0]), np.array([2.6679, 4.6210, 0]), EARTH_MU)
    r,v = ob.t2vectors(60*60)
    assert r == approx(np.array([-3296.8,7413.9,0]))
    assert v == approx(np.array([-8.2977,-0.96309,0]))


def test_curtis_4_3():
    r = np.array([-6045, -3490, 2500])
    v = np.array([-3.457, 6.618, 2.533])
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(58_310)
    assert ob.i == approx(np.radians(153.2))
    assert ob.raan == approx(np.radians(255.3))
    assert ob.e == approx(0.1712)
    assert ob.argp == approx(np.radians(20.07))
    assert ob.f(0) == approx(np.radians(28.45))

def test_curtis_4_7():
    
    ob = orbit_from_keplerian(-80_000,
                                1.4,
                                m.radians(30),
                                m.radians(40),
                                m.radians(60),
                                m.radians(30),
                                EARTH_MU)
    ob.h = 80_000
    ob.link_tf(0,m.radians(30))
    r,v = ob.t2vectors(0)
    r_facit = np.array([-4040,4815,3629])
    v_facit = np.array([-10.39, -4.772, 1.744])
    assert r == approx(r_facit)
    assert v == approx(v_facit)

def test_curtis_5_2():
    r1 = np.array([5000,10_000,2100])
    r2 = np.array([-14_600,2500,7000])
    dt = 60*60
    ob = orbit_from_lambert(r1,r2,0,dt,EARTH_MU)
    assert ob.h == approx(80_470)
    assert ob.a == approx(20_000)
    assert ob.e == approx(0.4335)
    assert ob.raan == approx(np.radians(44.60))
    assert ob.argp == approx(np.radians(30.71))
    assert ob.i == approx(np.radians(30.19))
    assert ob.r(0) == approx(4952+6378)
    assert ob.tp == approx(256.1)
    assert ob.f(0) == approx(np.radians(350.8))

    
def test_curtis_5_3():
    r1 = np.array([273_378,0,0])
    r2 = np.array([145_820,12_758,0])
    dt = 48_600
    ob = orbit_from_lambert(r1,r2,0,dt,EARTH_MU)
    p_alt = ob.r(0)
    t_p2 = ob.tp - dt
    assert p_alt == approx(160.2+6378)
    assert t_p2 == approx(38_396)

