''' 
Mathematical, physical and engineering constants.

All values given in the units:
length - km
mass   - kg
time   - s

'''
import math as m

# source: (https://ssd.jpl.nasa.gov/astro_par.html)
# physical constants:

G = 6.6743e-20 # [m^3 kg^-1 s^-2]
c = 299_792.458 # [km/s]

# solar system bodies constants:
EARTH_MU = 3.986004e5 # [km^3/s^2]
SUN_MU =   1.32712e11 # [km^3/s^2]
LUNA_MU =  4.9028e3   # [km^3/s^2]
MARS_MU =  4.2828e4   # [km^3/s^2]

# source: (https://ssd.jpl.nasa.gov/planets/phys_par.html)
EARTH_RADIUS = 6371 # (mean) [km]
MARS_RADIUS = 3389.5 # (mean) [km]

# time:
DAY = 86400 # [s]
SIDEREAL_DAY = 86164.0905 # [s]
JULIAN_YEAR = DAY*365.25 # [s]
JULIAN_CENTURY = JULIAN_YEAR*100 # [s]
MJD_OFFSET = 2_451_554.5 # [days]
# time between julian date epoch and J2000 epoch
# lengths:
AU = 149_597_870.7 # [km]
