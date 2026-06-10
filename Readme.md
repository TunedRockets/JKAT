<img src="Dont%20Ask.png" width="400" alt="IMGTEXT: A giant cat perched atop Saint Basil's Cathedral reaching for a rocket.
Taken from the cover art of Frank Hayes' Don't Ask">

# Johannes-Kepler Astrodynamical Toolbox
JKAT is a python package that makes astrodynamical calculations easy. Many astrodynamical problems are simple to formulate, but require lengthy calculation processes to solve. Several tools such as GMAT, Pykep, and TUDAT exist to make this easier, but those programs often has a certain barrier of entry (and aren't written by me). So as a fun project, JKAT aims to solve many of the same problems, but with an hopefully simpler API/interface, and quicker learning curve.

# installation
JKAT is available on the [Python Package Index](https://pypi.org/project/jkat/). To install it with pip, simply type `python -m pip install jkat`.
The source code is also available on [github](https://github.com/TunedRockets/JKAT).

# guide

> Warning!
This package is very much in alpha, expect regular changes, and even
total remakes of the API. do NOT use this code for anything except playing around
this is just published for fun, and will be tidied up at a later date


Guide coming soon...





# TODO:
 - [ ] flesh out tests for curtis chapter 6
 - [ ] remove duplication in transfer.direct

 - [ ] gravity assist optimizer
 - [ ] perturbations
 - [ ] gravity assists and assist planning
 - [x] hook into jpl-horizons for ephemerides (more accurate probably not needed)
 - [ ] figure out best way to represent horizons in examples (and figure out wierdness at t=0)

 - [ ] proper docstrings
 - [x] PyPI release
 - [ ] make imports like `from x import y`/`from x import *` consistent.

