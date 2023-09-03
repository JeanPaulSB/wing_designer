from naca import Naca
import matplotlib.pyplot as plt
import numpy as np
profile = Naca("0012",0.43,1000)

profile.addFlap(55,0.75)
profile.set_equations()
profile.compute()
profile.thinAirfoilTheory(0)



plt.plot(profile.x_coordinates[0],profile.y_coordinates[0],'b--',label = "airfoil surface")
plt.plot(profile.x_coordinates[1],profile.y_coordinates[1],'b--')
plt.plot(profile.x_coordinates[1],profile.camber(profile.x),'b--')
plt.show()
