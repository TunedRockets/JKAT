


from src.utils import *
from src.kep import *
import matplotlib.pyplot as plt
import numpy as np

def _yzval(d):
    '''calculate approximate root of lambert y function
    from values of d = r1+r2 / A'''
    yarr = [39.478, 9.869, 0,    -27.965, -44.651, -56.166, -64.706]
    darr = [-1.414, 0,     1.414, 10,     20,       30,      39.478]
    y = np.interp(d,darr,yarr)
    return y


def y(z):
    return (z*stumpff_s(z)-1) / m.sqrt(stumpff_c(z))

yy = []
dd = np.linspace(0.01-m.sqrt(2),39,1000)
for d in dd:
    z_guess = _yzval(d)
    y_guess = y(z_guess) # want y > -d
    yy.append(y_guess)

plt.plot(dd,yy, color='red')
plt.plot(dd,-dd, color='blue')
plt.show()
    