import typing
import random
import matplotlib.pyplot as plt

from PyQt5 import QtCore, QtWidgets,uic
from PyQt5.QtCore import Qt,QPointF
from PyQt5.QtGui import QPainter,QColor,QFont,QPen
from PyQt5.QtWidgets import QWidget,QVBoxLayout
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg 
from matplotlib.figure import Figure

class DrawingWindow(QWidget):
    def __init__(self, parent):
        QWidget.__init__(self,parent)
        self.airfoil = None
        fig = Figure(figsize=(5, 5))
        self.can = FigureCanvasQTAgg(fig)
        self.toolbar = NavigationToolbar(self.can, self)
        

        layout = QVBoxLayout(self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.can)

        self.ax = self.can.figure.add_subplot(111)
        self.ax.set_ylabel("y")
        self.ax.set_ylabel("x")


    def plot(self):
        # random data
        data = [random.random() for i in range(10)]
  

        if self.airfoil:
            self.airfoil.compute()

            # clearing old figure
            self.ax.cla()

            self.ax.plot(self.airfoil.x_coordinates[0],self.airfoil.y_coordinates[0],'b--',label = "airfoil surface")
            self.ax.plot(self.airfoil.x_coordinates[1],self.airfoil.y_coordinates[1],'b--')

            if not self.airfoil.symmetrical:
         
                self.ax.plot(self.airfoil.x,self.airfoil.yt,'-.y',label = "thickness distribution")
                self.ax.plot(self.airfoil.x,self.airfoil.yc,'-r',label = "mean camber line")
           
            

            self.ax.set_xlim(0,self.airfoil.chord )
            self.ax.set_ylim(-self.airfoil.thickness * 1.5,self.airfoil.thickness * 1.5)

            self.ax.set_title(f"NACA {self.airfoil.digits}")

            self.ax.legend()

            

            self.can.figure.tight_layout()
            self.can.draw()