from typing import Union
from scipy.integrate import simps
import matplotlib.pyplot as plt
import numpy as np


class Naca:
    # four digits naca
    def __init__(self,digits: Union[list,str], chord: Union[float,int],points: int) -> None:
        self.digits = digits
        self.symmetrical = False
        self.chord = chord
        self.points = points
        self.set_type()


        self.x_coordinates = ()
        self.y_coordinates = ()
        
    

    def export(self,filename = None):
        if not filename:
            filename = str(self.digits) +".csv"
        if filename:
            print(filename)
        print("Saving txt file with all the coordinates")
        with open(filename,"w") as outputfile:
            for x,y in zip(self.x_coordinates[0][::-1],self.y_coordinates[0][::-1]):
                print(f"{x},{y}",file = outputfile)     
            for x,y in zip(self.x_coordinates[1][1::],self.y_coordinates[1][1::]):
                print(f"{x},{y}",file = outputfile) 
                
    # source: https://archive.aoe.vt.edu/mason/Mason_f/CAtxtAppA.pdf
    def compute(self):
        if self.family == 4:
            # generating x values
            self.x = np.linspace(0,self.chord,self.points)

            self.yt = np.array(list(map(self.thickness_distribution,self.x)))

         


            if self.symmetrical:

                self.x_coordinates = (self.x,self.x)
                self.y_coordinates = (self.yt,-self.yt)
                
            if self.symmetrical == False:
                self.yc = np.array(list(map(self.camber_line,self.x)))

                self.dyc = np.array(list(map(self.camber_slope,self.x)))

                theta = np.arctan(self.dyc)

               

                xu = self.x - self.yt * np.sin(theta)
                xl = self.x + self.yt * np.sin(theta)

                yl = self.yc - self.yt * np.cos(theta)
                yu = self.yc + self.yt * np.cos(theta)


                self.x_coordinates = (xl,xu)

                self.y_coordinates = (yl,yu)
            
        if self.family == 5:
            pass
    
    def thickness_distribution(self,x):
        term1 =  0.2969 * (np.sqrt(x/self.chord))
        term2 = -0.1260 * (x / self.chord) 
        term3 = -0.3516 * np.power(x/self.chord,2)
        term4 =  0.2843 * np.power(x/self.chord,3)
        term5 = -0.1015 * np.power(x/self.chord,4)
        return self.thickness / 0.2 *  (term1 + term2 + term3 + term4 + term5)

  




       
    def camber_line(self, x ):
        if x >= 0 and x <=  self.camber_location * self.chord:
            return self.maximum_camber / (self.camber_location ** 2) * ((2 * self.camber_location * (x/self.chord) - (x/self.chord)**2))
        else:
            return  (self.maximum_camber * ((1 - 2 * self.camber_location) + 2 * self.camber_location * (x/self.chord) - (x/self.chord) ** 2)) / (1 - self.camber_location)**2

    def camber_slope(self,x):
        if x >=0 and x <=  self.camber_location * self.chord:
            return (2.0 * self.maximum_camber) / np.power(self.camber_location,2) * (self.camber_location - (x / self.chord))
        else:
            return (2.0 * self.maximum_camber) / np.power(1-self.camber_location,2) * (self.camber_location - (x / self.chord))




   
    # identifies the corresponding naca family
    def set_type(self):
        
        # four digit series
        if len(self.digits) == 4:
            self.family = 4
            self.maximum_camber = (int(self.digits[0]) / 100 ) 
            self.camber_location = (int(self.digits[1])/10) 
            self.thickness = (int(self.digits[2:])/ 100) 


      

            if self.maximum_camber == 0.0 and self.camber_location == 0.0:
             
                self.symmetrical =  True
            


        # five digit series
        elif len(self.digits) == 5:
            self.family = 5
            self.design_cl = (3 / 2) * int(self.digits[0]) / 10
            self.maximum_camber =  int(self.digits[1:3]) / 200
            self.thickness =  int(self.digits[-2:]) / 100
            

    

    def summary(self):

        if self.family == 5:
            print(f"""
            {self.family} digits NACA {self.digits}
            chord: {self.chord}
            maximum camber: {self.maximum_camber}
            thickness: {self.thickness}
            design lift coefficient: {self.design_cl}
            points: {self.points}

                """)
        elif self.family == 4:
            print(f"""
            {self.family} digits NACA {self.digits}
            chord: {self.chord}
            maximum camber: {self.maximum_camber}
            thickness: {self.thickness}
            camber location: {self.thickness}
            points: {self.points}



                """
                )
            




