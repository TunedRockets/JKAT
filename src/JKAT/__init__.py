''' 
global inits for accessing the basic functions like orbit and such
'''


# always want the kep stuff
from .kep import *

# also good to have some planets
from .ephemeris.examples import Earth, Mars, Jupiter


# and some basic consts
from .utils.consts import *