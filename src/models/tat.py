from typing import Union
from sympy import symbols, integrate, cos
from sympy.abc import x
import numpy as np


"""
ThinAirfoilThoery Solver class for any airfoil, it computes the aL_0 and assigns it to the airfoil.

"""


class ThinAirfoilTheory:
    @staticmethod
    def solve(airfoil):
        airfoil.characterize()
        airfoil.compute()
        # creating the theta symbol
        theta = symbols("theta")
        airfoil.a0 = 2 * np.pi
        # making the corresponding change of variable
        if airfoil.symmetrical != True:
            airfoil.slope_sym = airfoil.slope_sym.subs(
                x, airfoil.chord / 2 * (1 - cos(theta))
            )
            airfoil.aL_0 = (
                -1
                / np.pi
                * integrate(airfoil.slope_sym * (cos(theta) - 1), (theta, 0, np.pi))
            )

        if airfoil.symmetrical == True:
            if airfoil.flapped:
                phi = np.arccos(-2 * airfoil.flap_position + 1)
                airfoil.aL_0 = airfoil.flap_slope * (
                    1 - phi / np.pi + np.sin(phi) / np.pi
                )
            else:
                airfoil.aL_0 = 0
