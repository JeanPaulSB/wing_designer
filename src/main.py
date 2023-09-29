from geometry import Naca
from models import ThinAirfoilTheory, LiftingLineTheory
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

"""
First example
"""


from geometry import Naca
from models import ThinAirfoilTheory
import numpy as np

"""airfoil = Naca(
    "2412", 15, 50
)  # NACA 2412, 15m of chord and 50 points for its definition.
airfoil.addFlap(20, 0.85)
ThinAirfoilTheory.solve(airfoil)
# cl @ 5deg
cl = 2 * np.pi * (np.deg2rad(5) - airfoil.aL_0)
"""
