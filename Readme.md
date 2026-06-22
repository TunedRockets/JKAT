<img src="Dont%20Ask.png" width="400" alt="IMGTEXT: A giant cat perched atop Saint Basil's Cathedral reaching for a rocket.
Taken from the cover art of Frank Hayes' Don't Ask">

# Johannes-Kepler Astrodynamical Toolbox
JKAT is a python package that makes astrodynamical calculations easy. Many astrodynamical problems are simple to formulate, but require lengthy calculation processes to solve. Several tools such as GMAT, Pykep, and TUDAT exist to make this easier, but those programs often has a certain barrier of entry (and aren't written by me). So as a fun project, JKAT aims to solve many of the same problems, but with an hopefully simpler API/interface, and quicker learning curve.

JKAT started development in the summer 2025, and was made public in april 2026. current features include:
- creating and modifying keplerian orbits
- solving for position and velocity on an orbit given a time
- creating transfers between orbits using a lambert solver
- optimizing said transfer for different parameters
- interfacing with the JPL Horizons database
- plot resulting orbits

# installation
JKAT is available on the [Python Package Index](https://pypi.org/project/jkat/). To install it with pip, simply type `python -m pip install jkat`.
The source code is also available on [github](https://github.com/TunedRockets/JKAT).

# guide

> Warning!
This package is very much in alpha, expect regular changes.


Guide coming soon...





# TODO:
 - [ ] flesh out tests for curtis chapter 6
 - [ ] remove duplication in transfer.direct

 - [ ] gravity assist optimizer
 - [ ] perturbations
 - [ ] gravity assists and assist planning
 - [ ] figure out best way to represent horizons in examples (and figure out wierdness at t=0)
 - [ ] reset plot on window close
 - [ ] proper docstrings
 - [ ] make imports like `from x import y`/`from x import *` consistent.

