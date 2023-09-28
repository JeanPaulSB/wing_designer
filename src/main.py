from airfoil.naca import Naca
from models.tat import ThinAirfoilTheory
import matplotlib.pyplot as plt
import numpy as np

"""
myfoil = Naca("0015", 0.23, 100)
myfoil.compute()
print(myfoil.thinAirfoilTheory(0))
"""
from geometry import Naca

myfoil = Naca("2412", chord=15, points=50)
myfoil.addFlap(10, 0.80)

myfoil.characterize()
myfoil.compute()

ThinAirfoilTheory.solve(myfoil)

print(myfoil.aL_0 * 180 / np.pi)
