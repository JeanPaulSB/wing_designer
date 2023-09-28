from airfoil.naca import Naca
import matplotlib.pyplot as plt
import numpy as np

"""
myfoil = Naca("0015", 0.23, 100)
myfoil.compute()
print(myfoil.thinAirfoilTheory(0))
"""
from geometry import Naca

myfoil = Naca("0012", chord=15, points=50)
myfoil.addFlap(1, 0.80)
myfoil.characterize()
myfoil.compute()
myfoil.plot()
