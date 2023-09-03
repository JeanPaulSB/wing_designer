from naca import Naca
import matplotlib.pyplot as plt
import numpy as np
profile = Naca("2412",0.43,200)

profile.addFlap(5,0.90)
profile.set_equations()
profile.compute()
profile.thinAirfoilTheory(0)



plt.plot(profile.x_coordinates[0],profile.y_coordinates[0],'b--',label = "airfoil surface")
plt.plot(profile.x_coordinates[1],profile.y_coordinates[1],'b--')
plt.plot(profile.x_coordinates[1],profile.yc,'b--')
plt.show()
