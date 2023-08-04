import sys
import typing
import numpy as np
import random
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPainter, QColor, QFont, QPen
from PyQt5.QtWidgets import QWidget, QFileDialog

from naca import Naca


class Window(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self, None)
        uic.loadUi("gui.ui", self)
        self.setWindowTitle("Airfoil generator")  #

        self.plotButton.clicked.connect(self.validateInputs)
        self.exportButton.clicked.connect(self.getFile)

    def validateInputs(self):
        naca = self.nacaInput
        chord = self.chordInput
        points = self.pointsInput

        nacaisValid = False
        chordisValid = False
        pointsisValid = False

        if naca.text() != "":
            if len(naca.text()) == 4 or len(naca.text()) == 5:
                nacaisValid = True
            else:
                print("unaccepted naca")
            if self.containsLetters(naca.text()):
                print("just acepting numbers")

        if chord.text() != "":
            if self.containsLetters(chord.text()):
                print("Chord must be a number")

            else:
                chordisValid = True

        if points.text() != "":
            if self.containsLetters(points.text()):
                print("Points must be numerical")
            else:
                pointsisValid = True

        if all([nacaisValid, chordisValid, pointsisValid]):
            self.widget.airfoil = Naca(
                naca.text(), float(chord.text()), int(points.text())
            )
            self.widget.plot()

    # function that checks if an input contains any letter
    def containsLetters(self, string: str) -> bool:
        return any([c.isalpha() for c in string])

    def getFile(self):
        if self.widget.airfoil:
            options = QFileDialog.Options()
            filename = QFileDialog.getSaveFileName(
                self,
                "Save file",
                f"NACA{self.widget.airfoil.digits}.csv",
                options=options,
            )
            self.widget.airfoil.export(filename[0])


app = QtWidgets.QApplication(sys.argv)
# se crea la instancia de la ventana
miVentana = Window()
# se muestra la ventana
miVentana.show()
# se entrega el control al sistema operativo
sys.exit(app.exec_())
