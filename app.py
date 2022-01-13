# This Python file uses the following encoding: utf-8
import os
from pathlib import Path
import sys
import strictyaml as syaml

from PyQt6 import QtGui, uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow, QFileDialog
from PyQt6.QtCore import QFile

class ZEMAS(QMainWindow):
    def __init__(self,args=None):
        super(ZEMAS, self).__init__(args)
        self.load_ui()
        self.browseYaml.clicked.connect(self.browse_yaml)
        self.acquireImage.clicked.connect(self.acquire_image)

    def browse_yaml(self):
        path = os.fspath(Path(__file__).resolve().parent / "conf/")
        #open the file dialog for the
        dlg = QFileDialog()
        dlg.setDirectory(path)
        dlg.setFileMode(QFileDialog.FileMode.ExistingFiles)
        dlg.setNameFilter("YAML Configuration Files (*.yaml)")
        dlg.show()

    def acquire_image(self):
        path = os.fspath(Path(__file__).resolve().parent)
        #open the file dialog for the
        dlg = QFileDialog()
        dlg.setDirectory(path)
        dlg.setFileMode(QFileDialog.FileMode.ExistingFiles)
        dlg.setNameFilter("DM4 Files (*.dm4)")
        dlg.show()

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
