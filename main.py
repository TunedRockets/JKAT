


from src.utils import *
import matplotlib.pyplot as plt
import numpy as np


rvec,vvec = kep2vectors(
    0.5,
    0.3,
    1000,
    1000,
    0.4,
    0,
    0.3
)
print(vectors2kep(rvec,vvec,1000))