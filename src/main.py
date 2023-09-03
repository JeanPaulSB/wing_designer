from airfoil.naca import Naca
import matplotlib.pyplot as plt
import numpy as np


myfoil = Naca("4415", 0.23, 100)
myfoil.compute()
print(myfoil.thinAirfoilTheory(0))