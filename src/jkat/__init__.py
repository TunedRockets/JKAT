''' 
global inits for accessing the basic functions like orbit and such
'''


# always want the kep stuff
from .kep import *

# also good to have some planets
from .ephemeris.examples import Earth, Mars, Jupiter
from .ephemeris import examples

# and some basic consts
from .utils.consts import *


# and basic plotting:
from .plotting.plot import add_solar_system, plot
from .plotting import show

from .trajectories import direct_transfer