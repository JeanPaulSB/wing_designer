from airfoil.naca import Naca
import matplotlib.pyplot as plt
import numpy as np

"""
myfoil = Naca("0015", 0.23, 100)
myfoil.compute()
print(myfoil.thinAirfoilTheory(0))
"""
from geometry import Naca


myfoil = Naca("2412", chord=15, points=50)

print(myfoil.addFlap(30, 80))
