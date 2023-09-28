from geometry import Naca
from models import ThinAirfoilTheory
import matplotlib.pyplot as plt
import numpy as np

"""
First example
"""

myfoil = Naca("2412", chord=15, points=50)
myfoil.addFlap(10, 0.80)
myfoil.characterize()
myfoil.compute()
ThinAirfoilTheory.solve(myfoil)
print(myfoil.aL_0 * 180 / np.pi)
