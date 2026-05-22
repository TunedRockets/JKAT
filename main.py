


from src.utils import *
from src.kep import *
import src.plotting as oplt
import matplotlib.pyplot as plt
import numpy as np
from tests.test_timings import test_uni_round_trip_hyper
from tests.test_curtis import test_curtis_4_7



ob = Orbit(100,0.2,0.1,0.3,0.5,200,1000)
oplt.plot(ob, f=0)
oplt.center()
plt.axis('scaled')
oplt.show()