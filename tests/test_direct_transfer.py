


import math as m
import numpy as np

from src.jkat.utils.consts import YEAR
from src.jkat import Earth, Mars, Jupiter

import src.jkat.trajectories as tra

import matplotlib.pyplot as plt

def assert_optim_works(time_range, origin,destination):
    dvarr = tra.porkchop_plot(origin,destination,time_range,time_range, rendezvous=True, prograde=True)
    idx = np.unravel_index(np.nanargmin(dvarr),dvarr.shape)
    dv_cost = dvarr[idx]
    ts = time_range[idx[0]]
    te = time_range[idx[1]]

    res = tra.direct_transfer(origin,destination,(time_range[0],time_range[-1],time_range[0],time_range[-1]),dv1_w=1,dv2_w=1, prograde=True)
    
    if False:
        plt.imshow(dvarr.T,extent=(time_range[0],time_range[-1],time_range[0],time_range[-1]), origin="lower",vmax=30)
        plt.scatter(ts,te, label="prokchop min")
        plt.scatter(res['ts'],res['te'],label="optimizer")
        plt.legend()
        plt.show()
    
    dv_tot = res['dv1'] + res['dv2']
    assert dv_cost >= dv_tot

def test_opt_mars_2000():
    rang = np.linspace(0,3*YEAR)
    assert_optim_works(rang,Earth,Mars)

def test_opt_mars_2010():
    rang = np.linspace(10*YEAR,13*YEAR)
    assert_optim_works(rang,Earth,Mars)

def test_opt_jupiter_2000():
    rang = np.linspace(0,6*YEAR)
    assert_optim_works(rang,Earth,Jupiter)

def test_opt_jupiter_2010():
    rang = np.linspace(10*YEAR,16*YEAR)
    assert_optim_works(rang,Earth,Jupiter)