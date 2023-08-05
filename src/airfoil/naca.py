from typing import Union
from scipy.integrate import simps,trapz,cumtrapz
from sympy.abc import x
from sympy import Eq, Piecewise, GreaterThan,lambdify,cos,symbols,integrate,AccumBounds
import matplotlib.pyplot as plt
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
            self.set_equations()
            # computing thickness distribution
            self.yt = np.array(list(map(self.thickness_distribution, self.x)))

            if self.symmetrical:
                self.x_coordinates = (self.x, self.x)
                self.y_coordinates = (self.yt, -self.yt)

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
    
    def thickness_distribution(self,x):
        term1 =  0.2969 * (np.sqrt(x/self.chord))
        term2 = -0.1260 * (x / self.chord) 
        term3 = -0.3516 * np.power(x/self.chord,2)
        term4 =  0.2843 * np.power(x/self.chord,3)
        term5 = -0.1015 * np.power(x/self.chord,4)
        return self.thickness / 0.2 *  (term1 + term2 + term3 + term4 + term5)

  




       
   
    """
    Set the corresponding equations for the mean camber line and its slope
    """
    def set_equations(self):
        if self.family == 4:
            if not self.symmetrical:
                camber_eq1 = self.maximum_camber / (self.camber_location**2) * ( ( 2 * self.camber_location * (x/self.chord) - ( x / self.chord)**2))
                camber_eq2 = self.maximum_camber* ((1 - 2 * self.camber_location)+ 2 * self.camber_location * (x / self.chord) - (x / self.chord) ** 2) / (1 - self.camber_location) ** 2
                camber_slope_eq1 = 2 * self.maximum_camber / np.power(self.camber_location,2) * (self.camber_location - (x / self.chord))
                camber_slope_eq2 = 2 * self.maximum_camber / np.power(1-self.camber_location,2) * (self.camber_location - (x / self.chord))
                
                self.camber_sym = Piecewise((camber_eq1, x <= self.chord * self.camber_location),(camber_eq2,x > self.chord * self.camber_location))
                self.slope_sym = Piecewise((camber_slope_eq1, x <= self.chord * self.camber_location),(camber_slope_eq2,x > self.chord * self.camber_location))
                
                self.camber = lambdify(x,self.camber_sym)
                self.slope = lambdify(x,self.slope_sym)
                
                

    def thinAirfoilTheory(self,interval: tuple):
        # creating the theta symbol
        theta = symbols("theta")
        # making the corresponding change of variable
        self.camber_sym = self.camber_sym.subs(x,self.chord / 2 * (1 - cos(theta)))
        self.slope_sym = self.slope_sym.subs(x,self.chord / 2 * (1- cos(theta)))




        location_maximum_camber = np.arccos(abs(self.camber_location * 2  - 1 ))

    
        A0 = lambda alpha: alpha  - ( integrate(self.slope_sym,(theta,0,np.pi)) /(np.pi))
        A1 = 2 / np.pi * (integrate(self.slope_sym * cos(theta),(theta,0,np.pi)))
        A2 = 2 / np.pi * (integrate(self.slope_sym * cos(theta) * cos(theta),(theta,0,np.pi)))

        # now computing cl

        angles = np.linspace(self.deg2rad(interval[0]),self.deg2rad(interval[1]),20)


        cl_list = np.array(list(map(lambda rad: 2 * np.pi*(A0(rad)+A1/2),angles)))

        cm_ac = np.pi / 4 * (A2-A1)
        print(cm_ac)

         
 
  
      

       
        """
        # x in terms of theta
        theta_values = np.array(list(map(lambda x: np.arccos((2*x/self.chord) - 1),self.x)))
        
        derivative = self.slope_theta(theta_values)

        plt.plot(theta_values,derivative)
        print(trapz(derivative))
        plt.show()
        """
        
        
   

        """
        print(A0)
        A0 = lambdify(theta,A0)
        
        # f(b) - f(a) 
        print( - A0(0) + A0(np.pi) )
        """
        



    
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

            if self.maximum_camber == 0.0 and self.camber_location == 0.0:
                self.symmetrical = True

        # five digit series
        elif len(self.digits) == 5:
            self.family = 5
            self.design_cl = (3 / 2) * int(self.digits[0]) / 10
            self.maximum_camber = int(self.digits[1:3]) / 200
            self.thickness = int(self.digits[-2:]) / 100


    def deg2rad(self,deg):
        # 180 ° -> pi rad
        # x °
        return deg * np.pi / 180
    
    def summary(self):
        if self.family == 5:
            print(
                f"""
            {self.family} digits NACA {self.digits}
            chord: {self.chord}
            maximum camber: {self.maximum_camber}
            thickness: {self.thickness}
            design lift coefficient: {self.design_cl}
            points: {self.points}

                """
            )
        elif self.family == 4:
            print(
                f"""
            {self.family} digits NACA {self.digits}
            chord: {self.chord}
            maximum camber: {self.maximum_camber}
            thickness: {self.thickness}
            camber location: {self.thickness}
            points: {self.points}



                """
            )
