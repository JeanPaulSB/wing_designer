import numpy as np


class LiftingLineTheory:
    @staticmethod
    def solve(airfoil_root, airfoil_tip, b, c_root, c_tip, btip, alphageo, N):
        lamda = c_root / c_tip
        y = [-b / 2 * (1 - (2 * k - 1) / 2 * N) for k in range(N)]
        theta = [np.arccos(-2 * y / b)]
        c = [c_root * (1 - 2 * (lamda - 1) / b * y)]
