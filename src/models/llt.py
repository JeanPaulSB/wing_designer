from geometry import Airfoil
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt


class LiftingLineTheory:
    @staticmethod
    def solve(
        airfoil_root: Airfoil,
        airfoil_tip: Airfoil,
        b: float,
        c_root: float,
        c_tip: float,
        b_tip: float,
        alphageo: float,
        N: int,
        S: float,
        V: float
    ):
        lamda = c_tip / c_root
        y = np.zeros(N)
        c = np.zeros(N)
        beta = np.zeros(N)
        a0 = np.zeros(N)
        aL_0 = np.zeros(N)
        theta = np.zeros(N)
        D = np.zeros(N)
        C = np.zeros((N, N))

        h_tip = c_tip * np.sin(b_tip)
        AR = b**2 / S

        for i in range(0, N):
            i_dummy = i + 1
            y[i] = (-b / 2) * (1 - ((2 * i_dummy - 1) / (2 * N)))

        c = c_root * (1 - 2 * ((lamda - 1) / b) * y)
        beta = np.arcsin(-2 * y / (b * c) * h_tip)
        a0 = airfoil_root.a0 - 2 * ((airfoil_tip.a0 - airfoil_root.a0) / b) * y
        aL_0 = airfoil_root.aL_0 - 2 * ((airfoil_tip.aL_0 - airfoil_root.aL_0) / b) * y
        theta = np.arccos(-2 * y / b) 
        D = (alphageo - aL_0 + beta) * np.pi / 180


 

        for k_row in range(N):
            a = (4*b) / (a0[k_row]*c[k_row])
            for n_col in range(N):
                b_aux = (2*(n_col +1)-1) / np.sin(theta[k_row])
                C[k_row][n_col] = (a+b_aux) * np.sin((2*(n_col +1)-1)*theta[k_row])
                
        A = solve(C,D)
        CL = np.pi * AR * A[0]
        

        # computing drag
        delta = 0
        for i in range(1,N):
            
            delta  =(i+1) * (A[i]/A[0])**2

        e = 1 / (1+delta)

        # estimating CDi
        CDi = CL**2 / (np.pi *e * AR)

        # now getting stations in order to analyze the lift distribution
        y = np.insert(y,0,-b/2)
        y = np.append(y,0)
    

        # recomputing theta
        print(V)
        theta = np.arccos((-2*y)/b)  
        strength = 2*b*V*(A[0] * np.sin(theta) + A[1]*np.sin(3*theta))
        normalized_strength = strength / strength[-1]
        
        # redefining c 
        c = c_root * (1 - 2 * ((lamda - 1) / b) * y)
        CL_w = 2 * strength / (V*c)
        qinf = 0.5*1.225*V**2

        L = qinf*c*(CL_w/CL)

        plt.plot(y/(b/2),L)
        plt.show()

        
         