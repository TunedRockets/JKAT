'''
Tests that run the examples in curtis' book.
'''

import pytest
from pytest import approx
import math as m
import numpy as np

from src.jkat.kep import Orbit, orbit_from_keplerian, orbit_from_rv, orbit_from_lambert
from src.jkat.utils import *
    
from src.jkat.trajectories import *

REL = 2e-3

@pytest.fixture
def fixt_ob()->Orbit:
    '''standard earth orbit'''
    return orbit_from_keplerian(EARTH_RADIUS,0,0,0,0,0,EARTH_MU)


def test_curtis_2_5(fixt_ob:Orbit):
    ob = fixt_ob
    ob.T = SIDEREAL_DAY
    assert ob.apoapsis == approx(42_164, rel=REL)
    assert ob.vt(0) == approx(3.075, rel=REL)

def test_curtis_2_7(fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(4000+6378,400+6378)
    assert ob.e == approx(0.2098, rel=REL)
    assert ob.h == approx( 57_172, rel=REL)
    assert ob.a == approx(8578, rel=REL)
    assert ob.period == approx(2.196*(60*60), rel=REL)

def test_curtis_2_8():
    z1 = 1545
    f1 = m.radians(126)
    z2 = 852
    f2 = m.radians(58)
    z1 += 6378
    z2 += 6378 # earth radius used by Curtis
    mu = EARTH_MU
    e = 0.08164 # manually given since we can't solve non-linear equations
    h = 54_830
    # ensure it works:
    p = h2p(h, mu)
    assert r2f(z1,p,e) == approx(f1, rel=REL)
    assert f2r(f2, p,e) == approx(z2, rel=REL)

    # and the other answers:
    zp = f2r(0,p,e)
    za = f2r(m.pi,p,e)
    assert zp - 6378 == approx(595.5, rel=REL)
    a,_ = apse2ae(za,zp)
    assert a == approx(7593, rel=REL)
    T = a2T(a,mu)
    assert T == approx(6585,rel=REL)
    

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
    assert d == approx(13_270, rel=REL)


def test_curtis_2_10():
    
    r = np.array([14_600,0,0])
    v = elazr2vec(0,m.radians(90-50), 8.6)
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(80_708, rel=REL)
    assert ob.e == approx(1.3393, rel=REL)
    assert ob.periapsis == approx(6986, rel=REL)
    assert ob.a == approx(-20_590, rel=REL)
    assert ob.c3 == approx(19.36, rel=REL)
    assert ob.aiming_radius == approx(18_340, rel=REL)



def test_curtis_2_11(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 0.3
    ob.h = 60_000
    r,v = ob.vectors(np.radians(120))
    assert r == approx(np.array([-5312.7, 9201.9, 0]), rel=REL)
    assert v == approx(np.array([-5.7533, -1.3287, 0]), rel=REL)

def test_curtis_2_12():
    r = np.array([7000,9000,0])
    v = np.array([-3.3472,9.1251,0])
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(94_000, rel=REL)
    assert ob.f(0) == approx(m.radians(52.125), rel=REL)
    assert ob.e == approx(1.538, rel=REL)

def test_curtis_2_13_and_14():
    r = np.array([8182.4,-6865.9,0])
    v = np.array([0.47572,8.8116,0])
    ob = orbit_from_rv(r,v,EARTH_MU)
    f = ob.f(0)
    f2 = f + m.radians(120)
    r,v = ob.vectors(f2)
    r_facit = np.array([1454.9, 8251.6,0])
    v_facit = np.array([-8.1323,5.6785, 0])
    assert r == approx(r_facit, rel=REL)
    assert v == approx(v_facit, rel=REL)
    assert ob.e == approx(1.0563, rel=REL)
    assert f % (2*m.pi) == approx(m.radians(288.44), rel=REL)


def test_curtis_3_1(fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(9600,21000)
    T = ob.t(np.radians(120))
    assert T == approx(4077, rel=REL)

def test_curtis_3_2(fixt_ob:Orbit):
    ob = fixt_ob
    ob.set_apses(9600,21000)
    assert ob.e == approx(0.37255, rel=REL)
    theta = ob.f(3*60*60)
    assert theta == approx(np.radians(193.2), rel=REL)

def test_curtis_3_4(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 1
    # not testing the energy eq: TODO change
    ob.pe = 7972
    ob.tp = 0
    f = ob.f(6*60*60)
    r = ob.r(f)
    assert r == approx(86_899, rel=REL)

def test_curtis_3_5():
    ob = orbit_from_rv(np.array([6378+300,0,0]), np.array([0,15,0]),EARTH_MU)
    T = ob.t(np.radians(100))
    r = ob.r(np.radians(100))
    assert T == approx(4141.4, rel=REL)
    assert r == approx(48_497, rel=REL)
    r,v = ob.t2vectors(T + 3*60*60)
    assert np.linalg.norm(r) == approx(163_180, rel=REL)
    assert np.linalg.norm(v) == approx(10.51, rel=REL)

def test_curtis_3_6(fixt_ob:Orbit):
    ob = fixt_ob
    ob.e = 1.4682
    ob.h = 95_154
    ob.link_tf(0,m.radians(30))
    assert ob.r(m.radians(30)) == approx(10_000, rel=REL)
    assert ob.v(m.radians(30)) == approx(10, rel=REL)
    t = ob.f(1*60*60)
    assert t % (2*m.pi) == approx(m.radians(100.04), rel=REL)

def test_curtis_3_7():
    ob = orbit_from_rv(np.array([7000,-12_124, 0]), np.array([2.6679, 4.6210, 0]), EARTH_MU)
    r,v = ob.t2vectors(60*60)
    assert r == approx(np.array([-3296.8,7413.9,0]), rel=REL)
    assert v == approx(np.array([-8.2977,-0.96309,0]), rel=REL)

def test_curtis_3_7_alt():
    r = np.array([7000,-12_124, 0])
    v = np.array([2.6679, 4.6210, 0])
    dt = 60*60
    r1,v1 = propagate_vectors(r,v,dt,EARTH_MU)
    assert r1 == approx(np.array([-3296.8,7413.9,0]), rel=REL)
    assert v1 == approx(np.array([-8.2977,-0.96309,0]), rel=REL)


def test_curtis_4_3():
    r = np.array([-6045, -3490, 2500])
    v = np.array([-3.457, 6.618, 2.533])
    ob = orbit_from_rv(r,v,EARTH_MU)
    assert ob.h == approx(58_310, rel=REL)
    assert ob.i == approx(np.radians(153.2), rel=REL)
    assert ob.raan == approx(np.radians(255.3), rel=REL)
    assert ob.e == approx(0.1712, rel=REL)
    assert ob.argp == approx(np.radians(20.07), rel=REL)
    assert ob.f(0) == approx(np.radians(28.45), rel=REL)

def test_curtis_4_7():
    
    ob = orbit_from_keplerian(-80_000,
                                1.4,
                                m.radians(30),
                                m.radians(40),
                                m.radians(60),
                                m.radians(30),
                                EARTH_MU)
    assert ob.f(0) == approx(m.radians(30))
    ob.h = 80_000
    ob.link_tf(0,m.radians(30))
    assert ob.f(0) == approx(m.radians(30))
    
    r,v = ob.t2vectors(0)
    r_facit = np.array([-4040,4815,3629])
    v_facit = np.array([-10.39, -4.772, 1.744])
    assert r == approx(r_facit, rel=REL)
    assert v == approx(v_facit, rel=REL)

def test_curtis_5_2():
    r1 = np.array([5000,10_000,2100])
    r2 = np.array([-14_600,2500,7000])
    dt = 60*60
    ob = orbit_from_lambert(r1,r2,0,dt,EARTH_MU)
    assert ob.h == approx(80_470, rel=REL)
    assert ob.a == approx(20_000, rel=REL)
    assert ob.e == approx(0.4335, rel=REL)
    assert ob.raan == approx(np.radians(44.60), rel=REL)
    assert ob.argp == approx(np.radians(30.71), rel=REL)
    assert ob.i == approx(np.radians(30.19), rel=REL)
    assert ob.r(0) == approx(4952+6378, rel=REL)
    assert ob.f(0) == approx(np.radians(350.8), rel=REL)
    assert ob.tp % ob.T == approx(256.1, rel=REL)
    

    
def test_curtis_5_3():
    r1 = np.array([273_378,0,0])
    r2 = np.array([145_820,12_758,0])
    dt = 48_600
    ob = orbit_from_lambert(r1,r2,0,dt,EARTH_MU)
    p_alt = ob.r(0)
    t_p2 = ob.tp - dt
    assert p_alt == approx(160.2+6378, rel=REL)
    assert t_p2 == approx(38_396, rel=REL)

def test_curtis_6_1():
        # using trajectory optimizer
        origin = orbit_from_keplerian(1,0,0,0,0,0,EARTH_MU)
        origin.set_apses(800+6378, 480+6378)
        destination = orbit_from_keplerian(16_000+6378,0,0,0,0,0,EARTH_MU)


        res = orbit_transfer(origin,
                            destination,
                            prograde=True, dv2_w=1)
        
        assert res['dv1'] == approx(1.7225, rel=REL)
        assert res['dv2'] == approx(1.3297, rel=REL)

def test_curtis_6_1_alt():
    dv1,dv2,T = hohmann_transfer(480+6378, 16000+6378,EARTH_MU)
    
    # adjust for the elliptic orbit
    a,e = apse2ae(480+6378, 800+6378)
    p = a2p(a,e)
    h = p2h(p,EARTH_MU)
    dv3 = f2v(0,p,e,h) - vcirc(0,p,e,h)
    dv1 -= dv3


    assert dv1 == approx(1.7225, rel=REL)
    assert dv2 == approx(1.3297, rel=REL)


# need time-independent transfer between orbits

def test_curtis_6_2():
    # not right test for orbit transfer

    ob = Orbit(500+6378, 0,0,0,0,0,EARTH_MU)
    dt = a2T(apse2ae(5000+6378,500+6378)[0],EARTH_MU)/2
    f = ob.f(-dt)
    assert f == approx(m.radians(-275.2), rel=REL)
    
def test_curtis_6_5():

    dt, T = circular_phasing(-m.radians(12), 3, 42_164,EARTH_MU)
    assert 2*dt == approx(0.0227, rel=REL)

def test_curtis_6_6():
    raise NotImplementedError()
    # fpa and velocity diff

def test_curtis_6_7():
    raise NotImplementedError()
    # fpa and velocity diff

def test_curtis_6_8():
    raise NotImplementedError()
    # arbitrary v maneuver

def test_curtis_6_9():
    raise NotImplementedError()
    # complex trajectory

def test_curtis_6_11():
    raise NotImplementedError()
    # complex trajectory