from geometry import Airfoil, Naca
from models import LiftingLineTheory

naca0015 = Naca("0015", 1.84, 50)
naca0012 = Naca("0012", 0.83, 50)

naca0015.a0 = 6.436
naca0012.a0 = 6.363
naca0015.aL_0 = 0
naca0012.aL_0 = 0

LiftingLineTheory.solve(naca0015, naca0012, 8, 1.84, 0.83, 0, 2, 4,10.7,317 * 10/36)
