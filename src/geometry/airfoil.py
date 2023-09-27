from abc import ABC, abstractmethod
from typing import Union
import numpy as np


class Airfoil(ABC):
    @abstractmethod
    def thickness_distribution(self, x):
        pass

    """
    Adds a flap to the airfoil at a specified position
    in terms of x/c.
    """

    def addFlap(self, angle: float, position: float):
        self.flap_angle = np.deg2rad(angle)
        self.flap_position = position
        self.flapped = True

    def export(self, filename: str):
        if not filename:
            if self.digits:
                filename = self.digits + ".csv"
            if self.name:
                filename = self.name + ".csv"
            else:
                raise ValueError(
                    "Either provide a name to your airfoil or assign digits if it is a NACA one."
                )
        with open(filename, "w") as output_file:
            for x, y in zip(self.x_coordinates[0][::-1], self.y_coordinates[0][::-1]):
                print(f"{x},{y}", file=output_file)
            for x, y in zip(self.x_coordinates[1][1::], self.y_coordinates[1][1::]):
                print(f"{x},{y}", file=output_file)
