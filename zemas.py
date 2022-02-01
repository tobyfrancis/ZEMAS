# This Python file uses the following encoding: utf-8
import os
from pathlib import Path
import sys
import strictyaml as syaml

import numpy as np
import cv2, imutils

from include.dm3_lib import DM3
from include.crystals.parsers import CIFParser
from include.workers import *

from PyQt6 import QtGui, uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow, QFileDialog
from PyQt6.QtCore import QFile, QThreadPool


class ZEMAS(QMainWindow):
    def __init__(self,args=None):
        super(ZEMAS, self).__init__(args)

        # Global variables
        self.cif_data      = None
        self.current_image = None
        self.current_hkl   = None
        self.current_uvw   = None
        self.current_ormat = None

        self.alpha   = 0
        self.beta    = 0
        self.imrot   = 0
        self.diffrot = 0

        # Thread pool for paralllelization
        self.threadpool = QThreadPool.globalInstance()
        
        # Load the UI from file, providing access to named elements
        self.load_ui()

        # ------ Configuration tab ------
        self.browseCIF.clicked.connect(self.browse_cif)

        self.imRotDial.valueChanged.connect(self.set_im_rotation)
        self.diffRotDial.valueChanged.connect(self.set_diff_rotation)
        self.imRotSpinBox.valueChanged.connect(self.set_im_rotation)
        self.diffRotSpinBox.valueChanged.connect(self.set_diff_rotation)

        # ------ Tilting and acquisition tab ------
        self.acquireImage.clicked.connect(self.acquire_image)
        self.calcZoneAxis.clicked.connect(self.calc_zone_axis)
        self.moveToZoneAxis.clicked.connect(self.calc_alpha_beta)

        self.alphaDial.valueChanged.connect(self.set_alpha_spinbox)
        self.betaDial.valueChanged.connect(self.set_beta_spinbox)
        self.alphaDial.sliderReleased.connect(self.set_alpha_released)
        self.betaDial.sliderReleased.connect(self.set_beta_released)
        self.alphaSpinBox.valueChanged.connect(self.set_alpha_released)
        self.betaSpinBox.valueChanged.connect(self.set_beta_released)

    def load_ui(self):
        path = os.fspath(Path(__file__).resolve().parent / "form.ui")
        uic.loadUi(path, self)

    # ------ Configuration tab ------
    def browse_cif(self):
        path = os.fspath(Path(__file__).resolve().parent / "conf/")
        #open the file dialog for the
        filename, _ = QFileDialog.getOpenFileName(
                      self,
                      "Open file",
                      "",
                      "CIF Files (*.cif)",
        )
        self.cif_data = CIFParser(filename)

    def set_im_rotation(self, value):
        self.imrot   = value
        self.imRotDial.setValue(value)
        self.imRotSpinBox.setValue(self.imrot)

    def set_diff_rotation(self, value):
        self.diffrot = value
        self.diffRotDial.setValue(value)
        self.diffRotSpinBox.setValue(self.diffrot)

    # ------ Tilting and acquisition tab ------

    def set_alpha_spinbox(self,value):
        self.alphaSpinBox.setValue(value)

    def set_beta_spinbox(self,value):
        self.betaSpinBox.setValue(value)

    def set_alpha_released(self):
        self.alpha = self.alphaSpinBox.value()
        self.alphaDial.setValue(self.alpha)

    def set_beta_released(self):
        self.beta  = self.betaSpinBox.value()
        self.betaDial.setValue(self.beta)

    def acquire_image(self):
        #open the file dialog for the
        filename, _ = QFileDialog.getOpenFileName(
                      self,
                      "Open file",
                      "",
                      "Tiff File (*.tiff);;DM3 File (*.dm3)",
        )

        if filename is not None:
            if filename[-3:] == 'dm3':
                dm3 = DM3(filename)
                image = dm3.imagedata
                self.current_image = image

                #display code
                frame = imutils.resize(image, width=480)
                height, width  = frame.shape
                bytes_per_line = width
                qt_image = QtGui.QImage(frame, width, height, bytes_per_line, QtGui.QImage.Format.Format_Indexed8)
                self.currentImage.setPixmap(QtGui.QPixmap.fromImage(qt_image))

            if filename[-4:] == 'tiff':
                image = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
                self.current_image = image

                #display code
                frame = imutils.resize(image, width=480)
                height, width  = frame.shape
                bytes_per_line = width
                qt_image = QtGui.QImage(frame, width, height, bytes_per_line, QtGui.QImage.Format.Format_Indexed8)
                self.currentImage.setPixmap(QtGui.QPixmap.fromImage(qt_image))

    def calc_zone_axis(self):
        image_analyzer = SpotDiffractionImageAnalyzer(self.current_image, self.diffrot, self.alpha, self.beta, self.cif_data)
        self.threadpool.start(image_analyzer)
        image_analyzer.signals.result.connect(self.set_hkl_uvw_ormat)

    def set_hkl_uvw_ormat(self,results):
        hkl,uvw,ormat = results
        self.current_hkl   = hkl
        self.current_uvw   = uvw
        self.current_ormat = ormat

        self.uBox1.setValue(uvw[0])
        self.vBox1.setValue(uvw[1])
        self.wBox1.setValue(uvw[2])

    def calc_alpha_beta(self):
        u1 = self.uBox1.value()
        v1 = self.vBox1.value()
        w1 = self.wBox1.value()
        uvw1 = np.array([u1, v1, w1])
        u2 = self.uBox2.value()
        v2 = self.vBox2.value()
        w2 = self.wBox2.value()
        uvw2 = np.array([u2, v2, w2])
        calculator = AlphaBetaTilter(self.alpha, self.beta, uvw1, uvw2, self.current_ormat, self.diffrot, self.cif_data)
        self.threadpool.start(calculator)
        calculator.signals.result.connect(self.set_alpha_beta)

    def set_alpha_beta(self,results):
        alpha, beta = results
        self.alphaDial.setValue(alpha)
        self.betaDial.setValue(beta)
        self.alphaSpinBox.setValue(alpha)
        self.betaSpinBox.setValue(beta)
        self.alpha = alpha
        self.beta  = beta

def run_zemas():
    #Initialize QT application
    app = QApplication([])
    #Set application icon
    app.setWindowIcon(QtGui.QIcon('assets/zonexus.png'))

    #TODO: show splash screen
    #TODO: during splash screen, load YAML files from conf folder and initialize Microscopes/Holders

    #load and set the stylesheet for the application
    path = os.fspath(Path(__file__).resolve().parent / "include/style.qss")
    with open(path, 'r') as f:
        style = f.read()
        app.setStyleSheet(style)

    #start the main ZEMAS widget, subclassed from QMainWindow
    widget = ZEMAS()
    widget.setWindowTitle("ZEMAS")
    widget.setWindowIcon(QtGui.QIcon('assets/zonexus.png'))
    #show widget (using show function from QMainWindow)
    widget.show()
    #exit application
    sys.exit(app.exec())

if __name__ == "__main__":
    run_zemas()
