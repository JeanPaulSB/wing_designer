from typing import Union

from sympy.abc import x
from sympy import (
    Eq,
    Piecewise,
    GreaterThan,
    lambdify,
    cos,
    symbols,
    integrate,
    AccumBounds,
    simplify,
)
import numpy as np


class Naca:
    # four digits naca
    def __init__(
        self, digits: Union[list, str], chord: Union[float, int], points: int
    ) -> None:
        self.digits = digits
        self.symmetrical = False
        self.chord = chord
        self.points = points
        self.set_type()
        self.x_coordinates = ()
        self.y_coordinates = ()
        self.flapped = False
        self.flap_angle = 0
        self.flap_position = 0

    """
    Exports the airfoil coordinates in a csv file
    """

    def export(self, filename=None):
        if not filename:
            filename = str(self.digits) + ".csv"
        if filename:
            print(filename)
        print("Saving txt file with all the coordinates")
        with open(filename, "w") as outputfile:
            for x, y in zip(self.x_coordinates[0][::-1], self.y_coordinates[0][::-1]):
                print(f"{x},{y}", file=outputfile)
            for x, y in zip(self.x_coordinates[1][1::], self.y_coordinates[1][1::]):
                print(f"{x},{y}", file=outputfile)

    # source: https://archive.aoe.vt.edu/mason/Mason_f/CAtxtAppA.pdf
    """
    Performs the corresponding geometrical computations for the NACA airfoil
    using the mean camber line, its slope and the thickness distribution whether the airfoil
    is symmetrical or not.
    """

    def compute(self):
        if self.family == 4:
            # generating x values
            self.x = np.linspace(0, self.chord, self.points)
            # computing thickness distribution
            self.yt = np.array(list(map(self.thickness_distribution, self.x)))

            if self.symmetrical:
                self.x_coordinates = (self.x, self.x)
                self.y_coordinates = (self.yt, -self.yt)

                # NOTE: symmetrical airfoils do not have camber, but when a flap its added, it creates a cambered section with the form
                # of a straigth line from the position of the hing point of the flap to the trailing edge.
                if self.flapped:
                    self.yc = self.camber(self.x)

                    self.dyc = self.slope(self.x)

                    theta = np.arctan(self.dyc)

                    xu = self.x - self.yt * np.sin(theta)
                    xl = self.x + self.yt * np.sin(theta)

                    yl = self.yc - self.yt * np.cos(theta)
                    yu = self.yc + self.yt * np.cos(theta)

                    self.x_coordinates = (xl, xu)

                    self.y_coordinates = (yl, yu)

            if self.symmetrical == False:
                # computing camber line
                self.yc = self.camber(self.x)
                # camber slope
                self.dyc = self.slope(self.x)

                theta = np.arctan(self.dyc)

                # airfoil surface coordinates
                xu = self.x - self.yt * np.sin(theta)
                xl = self.x + self.yt * np.sin(theta)

                yl = self.yc - self.yt * np.cos(theta)
                yu = self.yc + self.yt * np.cos(theta)

                self.x_coordinates = (xl, xu)

                self.y_coordinates = (yl, yu)

        if self.family == 5:
            # TODO: pending
            pass

    def thickness_distribution(self, x):
        term1 = 0.2969 * (np.sqrt(x / self.chord))
        term2 = -0.1260 * (x / self.chord)
        term3 = -0.3516 * np.power(x / self.chord, 2)
        term4 = 0.2843 * np.power(x / self.chord, 3)
        term5 = -0.1015 * np.power(x / self.chord, 4)
        return self.thickness / 0.2 * (term1 + term2 + term3 + term4 + term5)

    """
    adds a flap to the airfoil
    angle: deg
    position: in percent of the cord   
    """

    def addFlap(self, angle: float, position: float):
        self.flap_angle = self.deg2rad(angle)
        self.flap_position = position
        self.flapped = True

    """
    Set the corresponding equations for the mean camber line and its slope
    for the thin airfoil theory
    """

    def set_equations(self):
        if self.family == 4:
            if not self.symmetrical:
                # defining symbolic camber and camber slope equations
                self.camber_eq1 = (
                    self.maximum_camber
                    / (self.camber_location**2)
                    * (
                        (
                            2 * self.camber_location * (x / self.chord)
                            - (x / self.chord) ** 2
                        )
                    )
                )
                self.camber_eq2 = (
                    self.maximum_camber
                    * (
                        (1 - 2 * self.camber_location)
                        + 2 * self.camber_location * (x / self.chord)
                        - (x / self.chord) ** 2
                    )
                    / (1 - self.camber_location) ** 2
                )
                self.camber_slope_eq1 = (
                    2
                    * self.maximum_camber
                    / np.power(self.camber_location, 2)
                    * (self.camber_location - (x / self.chord))
                )
                self.camber_slope_eq2 = (
                    2
                    * self.maximum_camber
                    / np.power(1 - self.camber_location, 2)
                    * (self.camber_location - (x / self.chord))
                )
                # integrating the previous equations into a piecewise function that describes the airfoil camber and its slope taking into account
                # the position of the maximum camber
                self.camber_sym = Piecewise(
                    (self.camber_eq1, x <= self.chord * self.camber_location),
                    (self.camber_eq2, x > self.chord * self.camber_location),
                )

                self.slope_sym = Piecewise(
                    (self.camber_slope_eq1, x <= self.chord * self.camber_location),
                    (self.camber_slope_eq2, x > self.chord * self.camber_location),
                )
                # making those symbolic equations actually evaluables
                self.camber = lambdify(x, self.camber_sym)
                self.slope = lambdify(x, self.slope_sym)

                if self.flapped:
                    # the slope it's actually negative
                    self.flap_slope = -(self.flap_angle)
                    # defining graph for the flap eq in order to intersect the corresponding camber line
                    # y = mx+b
                    # b = y -mx
                    b = self.camber(
                        self.flap_position * self.chord
                    ) - self.flap_slope * (self.flap_position * self.chord)
                    flap_eq = self.flap_slope * x + b
                    #  now we have to redifine the the camber equations knowing the position of the flap and the same with the slope sym eq
                    self.camber_sym = Piecewise(
                        (self.camber_eq1, x <= self.chord * self.camber_location),
                        (
                            self.camber_eq2,
                            (x > self.chord * self.camber_location)
                            & (x < self.chord * self.flap_position),
                        ),
                        (flap_eq, (x >= self.chord * self.flap_position)),
                    )
                    self.slope_sym = Piecewise(
                        (
                            self.camber_slope_eq1,
                            (x >= 0) & (x <= self.chord * self.camber_location),
                        ),
                        (
                            self.camber_slope_eq2,
                            (x >= self.chord * self.camber_location)
                            & (x < self.chord * self.flap_position),
                        ),
                        (
                            self.flap_slope,
                            (x >= self.chord * self.flap_position) & (x <= self.chord),
                        ),
                    )
                    # repeating this step since our piecewise functions have been redefined with the addition of the high lift device
                    self.camber = lambdify(x, self.camber_sym)
                    self.slope = lambdify(x, self.slope_sym)

            if self.symmetrical:
                if self.flapped:
                    # NOTE: this step is just an "extra" since the model for symmetrical airfoils do not require all these things.
                    # we are finding the flap eq just to plot it but only the slope is necessary for the TaT computations.
                    # we still need to find the incercept of the line with the y axis
                    self.flap_slope = -(self.flap_angle)
                    b = -self.flap_slope * self.chord * self.flap_position
                    self.flap_eq = self.flap_slope * x + b
                    # just redefining it
                    self.camber_sym = Piecewise(
                        (0, x <= self.chord * self.flap_position),
                        (self.flap_eq, x > self.chord * self.flap_position),
                    )
                    self.slope_sym = Piecewise(
                        ((0, x < self.chord * self.flap_position)),
                        (self.flap_slope, x >= self.chord * self.flap_position),
                    )
                    self.camber = lambdify(x, self.camber_sym)
                    self.slope = lambdify(x, self.slope_sym)

    def thinAirfoilTheory(self, alpha: Union[float, tuple]):
        # creating the theta symbol
        theta = symbols("theta")
        # making the corresponding change of variable
        if self.symmetrical != True:
            if self.flapped:
                # making the substitution in order to continue integrating
                self.slope_sym = self.slope_sym.subs(
                    x, self.chord / 2 * (1 - cos(theta))
                )
                # computing alpha @ L=0
                self.aL_0 = (
                    -1
                    / np.pi
                    * integrate(self.slope_sym * (cos(theta) - 1), (theta, 0, np.pi))
                )
            else:
                self.slope_sym = self.slope_sym.subs(
                    x, self.chord / 2 * (1 - cos(theta))
                )
                self.aL_0 = (
                    -1
                    / np.pi
                    * integrate(self.slope_sym * (cos(theta) - 1), (theta, 0, np.pi))
                )

        if self.symmetrical == True:
            if self.flapped:
                phi = np.arccos(-2 * self.flap_position + 1)
                self.aL_0 = self.flap_slope * (1 - phi / np.pi + np.sin(phi) / np.pi)
            else:
                self.aL_0 = 0

        if isinstance(alpha, float) or isinstance(alpha, int):
            return (alpha, 2 * np.pi * (np.deg2rad(alpha) - self.aL_0))
        if isinstance(alpha, tuple):
            alphas = np.arange(alpha[0], alpha[-1])
            cl_list = [2 * np.pi * (np.deg2rad(elem) - self.aL_0) for elem in alphas]
            return (alphas, cl_list)

    """
    Identifies the corresponding NACA parameters, its family and characterizes the airfoil
    """

    def set_type(self):
        # four digit series
        if len(self.digits) == 4:
            self.family = 4
            self.maximum_camber = int(self.digits[0]) / 100
            self.camber_location = int(self.digits[1]) / 10
            self.thickness = int(self.digits[2:]) / 100
            if self.maximum_camber == 0.0 or self.camber_location == 0.0:
                self.symmetrical = True

        # five digit series
        elif len(self.digits) == 5:
            self.family = 5
            self.design_cl = (3 / 2) * int(self.digits[0]) / 10
            self.maximum_camber = int(self.digits[1:3]) / 200
            self.thickness = int(self.digits[-2:]) / 100

    def deg2rad(self, deg):
        # 180 ° -> pi rad
        # x °
        return deg * np.pi / 180
