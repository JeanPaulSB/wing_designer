from typing import Union

from sympy.abc import x
from sympy import Eq, Piecewise, GreaterThan,lambdify,cos,symbols,integrate,AccumBounds,simplify
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
    adds a flap to the airfoil
    angle: deg
    position: in percent of the cord   
    """
    def addFlap(self,angle:float,position: float):
        self.flap_angle = angle
        self.flap_position = position
        self.flapped  = True
       
   
    """
    Set the corresponding equations for the mean camber line and its slope
    """
    def set_equations(self):
        if self.family == 4:
            if not self.symmetrical:
                self.camber_eq1 = self.maximum_camber / (self.camber_location**2) * ( ( 2 * self.camber_location * (x/self.chord) - ( x / self.chord)**2))
                self.camber_eq2 = self.maximum_camber* ((1 - 2 * self.camber_location)+ 2 * self.camber_location * (x / self.chord) - (x / self.chord) ** 2) / (1 - self.camber_location) ** 2


                self.camber_slope_eq1 = 2 * self.maximum_camber / np.power(self.camber_location,2) * (self.camber_location - (x / self.chord))
                self.camber_slope_eq2 = 2 * self.maximum_camber / np.power(1-self.camber_location,2) * (self.camber_location - (x / self.chord))

               
                
                self.camber_sym = Piecewise((self.camber_eq1, x <= self.chord * self.camber_location),
                                            (self.camber_eq2,x > self.chord * self.camber_location))
                
                self.slope_sym = Piecewise((self.camber_slope_eq1, x <= self.chord * self.camber_location),
                                           (self.camber_slope_eq2,x > self.chord * self.camber_location))

                
                self.camber = lambdify(x,self.camber_sym)
                self.slope = lambdify(x,self.slope_sym)

                if self.flapped:
                    # getting phi in radians
                    phi = np.arccos(2*self.flap_position-1)
                    
                    self.flap_slope = -(self.flap_angle)
                   



                    # defining graph for the flap eq
                    # y = mx+b
                    # b = y -mx
                    
                    b = self.camber(self.flap_position * self.chord) - self.flap_slope * (self.flap_position * self.chord) 
                    flap_eq = self.flap_slope *  x + b


                    #  now we have to redifine the the camber equations knowing the position of the flap and the same with the slope sym eq
                    self.camber_sym = Piecewise((self.camber_eq1, x <= self.chord * self.camber_location),
                                                (self.camber_eq2,(x > self.chord * self.camber_location) & (x < self.chord * self.flap_position))  ,
                                                (flap_eq, (x >= self.chord * self.flap_position)))
                    
                    self.slope_sym = Piecewise((self.camber_slope_eq1, (x >= 0 ) & (x <= self.chord * self.camber_location)),
                                               (self.camber_slope_eq2,(x >= self.chord * self.camber_location) & (x < self.chord * self.flap_position)),
                                               (self.deg2rad(self.flap_slope),(x >= self.chord * self.flap_position) & (x <= self.chord)))
                    

            
                self.camber = lambdify(x,self.camber_sym)
                self.slope = lambdify(x,self.slope_sym)
                  




                
                

    def thinAirfoilTheory(self,alpha):
        # creating the theta symbol
        theta = symbols("theta")
        # making the corresponding change of variable

        convertorad = lambda x: np.arccos(-2*x + 1) 
        

        if self.flapped:
            integral_limits = (
            convertorad(self.camber_location),
            convertorad(self.flap_position),
            np.pi )
            print(integral_limits)
            self.slope_sym = self.slope_sym.subs(x,self.chord / 2 * (1- cos(theta)))
            aL_0 = - 1 / np.pi * integrate(self.slope_sym * (cos(theta) - 1 ),(theta,0,np.pi))
            print(simplify(aL_0 * 180 / np.pi).evalf())



        else:
            self.slope_sym = self.slope_sym.subs(x,self.chord / 2 * (1- cos(theta)))
            print(self.slope_sym)
            aL_0 = - 1 / np.pi * integrate(self.slope_sym * (cos(theta) - 1 ),(theta,0,np.pi))
            print((aL_0 * 180 / np.pi))

        
        

    


      





    
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
    