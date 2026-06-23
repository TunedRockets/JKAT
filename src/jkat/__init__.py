''' 
global inits for accessing the basic functions like orbit and such
'''

# also good to have some planets
from .ephemeris.examples import \
    Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
from .ephemeris import from_time, to_time, horizons_request

# always want the kep stuff
from .kep import *



# and some basic consts
from .utils.consts import *


# and basic plotting:
from .plotting import show, add_solar_system, plot, clf

from .trajectories import direct_transfer