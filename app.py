# This Python file uses the following encoding: utf-8
import os
from pathlib import Path
import sys
import strictyaml as syaml

import numpy as np
import cv2, imutils

from include.dm3_lib import DM3
from include.workers import *

from PyQt6 import QtGui, uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow, QFileDialog
from PyQt6.QtCore import QFile, QThreadPool


class ZEMAS(QMainWindow):
    def __init__(self,args=None):
        super(ZEMAS, self).__init__(args)

        #private variables 
        self.current_image = None
        self.threadpool = QThreadPool.globalInstance()
        self.image_analyzer = ImageAnalyzer()
        
        #function calls
        self.load_ui()
        self.browseYaml.clicked.connect(self.browse_yaml)
        self.acquireImage.clicked.connect(self.acquire_image)
        self.calcZoneAxis.clicked.connect(self.calc_zone_axis)
        

    def browse_yaml(self):
        path = os.fspath(Path(__file__).resolve().parent / "conf/")
        #open the file dialog for the
        filename, _ = QFileDialog.getOpenFileName(
                      self,
                      "Open file",
                      "",
                      "YAML Configuration Files (*.yaml)",
        )
        #Todo: load with strictyaml

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
                frame = imutils.resize(image,width=480)
                frame = cv2.cvtColor(frame, cv2.COLOR_GRAY2RGB)
                height, width, channels  = frame.shape
                bytes_per_line = channels * width
                qt_image = QtGui.QImage(frame, width, height, bytes_per_line, QtGui.QImage.Format.Format_Indexed8)
                self.currentImage.setPixmap(QtGui.QPixmap.fromImage(qt_image))
                

            if filename[-4:] == 'tiff':
                image = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
                self.current_image = image

                #display code
                frame = imutils.resize(image,width=480)
                frame = cv2.cvtColor(frame, cv2.COLOR_GRAY2RGB)
                height, width, channels  = frame.shape
                bytes_per_line = channels * width
                qt_image = QtGui.QImage(frame, width, height, bytes_per_line, QtGui.QImage.Format.Format_Indexed8)
                self.currentImage.setPixmap(QtGui.QPixmap.fromImage(qt_image))

    def calc_zone_axis(self):
        if self.image_analyzer.image != self.current_image:
            self.image_analyzer.set_image(self.current_image)
            self.threadpool.start(self.image_analyzer)
            self.image_analyzer.signals.finished.connect(self.set_hkl)

    def set_hkl(self):
        self.hBox1.setValue(self.image_analyzer.hkl[0])
        self.kBox1.setValue(self.image_analyzer.hkl[1])
        self.lBox1.setValue(self.image_analyzer.hkl[2])

    def load_ui(self):
        path = os.fspath(Path(__file__).resolve().parent / "form.ui")
        uic.loadUi(path, self)

if __name__ == "__main__":
    #Initialize QT application
    app = QApplication([])
    #Set application icon
    app.setWindowIcon(QtGui.QIcon('assets/zonexus.png'))

    #TODO: show splash screen
    #TODO: during splash screen, load YAML files from conf folder and initialize Microscopes/Holders

    #load and set the stylesheet for the application
    path = os.fspath(Path(__file__).resolve().parent / "style.qss")
    with open(path, 'r') as f:
        style = f.read()
        app.setStyleSheet(style)

    #start the main ZEMAS widget, subclassed from QMainWindow
    widget = ZEMAS()
    widget.setWindowTitle("ZEMAS")
    #show widget (using show function from QMainWindow)
    widget.show()
    #exit application
    sys.exit(app.exec())
