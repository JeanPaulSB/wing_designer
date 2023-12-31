import numpy as np

from .airfoil import Airfoil
from typing import Union

from sympy.abc import x
from sympy import Eq, Piecewise, lambdify


class Naca(Airfoil):
    def __init__(self, digits: str, chord: Union[float, int], points: int) -> None:
        self.digits = digits
        self.chord = chord
        self.points = points
        self.x_coordinates = ()
        self.y_coordinates = ()
        self.flapped = False
        self.symmetrical = False
        self.flap_angle = 0
        self.flap_position = 0
        self.aL_0 = None

    """
    Sets the following parameters:
    x: list of x values from 0 to c
    yt: thickness distribution
    yc: camber
    dyc: camber slope
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

    """
    Returns the thickness distribution y(t) at a specific point
    """
    # source: https://archive.aoe.vt.edu/mason/Mason_f/CAtxtAppA.pdf

    def thickness_distribution(self, x) -> float:
        term1 = 0.2969 * (np.sqrt(x / self.chord))
        term2 = -0.1260 * (x / self.chord)
        term3 = -0.3516 * np.power(x / self.chord, 2)
        term4 = 0.2843 * np.power(x / self.chord, 3)
        term5 = -0.1015 * np.power(x / self.chord, 4)
        return self.thickness / 0.2 * (term1 + term2 + term3 + term4 + term5)

    """
    Set the corresponding equations in terms of the family (4-series or 5-series)
    """

    def characterize(self):
        if len(self.digits) == 4:
            self.family = 4
            self.maximum_camber = int(self.digits[0]) / 100
            self.camber_location = int(self.digits[1]) / 10
            self.thickness = int(self.digits[2:]) / 100
            if self.maximum_camber == 0.0 or self.camber_location == 0.0:
                self.symmetrical = True

        if len(self.digits) == 5:
            self.family = 5
            self.design_cl = (3 / 2) * int(self.digits[0]) / 10
            self.maximum_camber = int(self.digits[1:3]) / 200
            self.thickness = int(self.digits[-2:]) / 100

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
