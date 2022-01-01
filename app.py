# This Python file uses the following encoding: utf-8
import os
from pathlib import Path
import sys

from PyQt6 import QtGui, uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow
from PyQt6.QtCore import QFile

class App(QMainWindow):
    def __init__(self,args=None):
        super(App, self).__init__(args)
        self.load_ui()

    def load_ui(self):
        path = os.fspath(Path(__file__).resolve().parent / "form.ui")
        uic.loadUi(path, self)

if __name__ == "__main__":

    app = QApplication([])
    app.setWindowIcon(QtGui.QIcon('assets/zonexus.png'))
    path = os.fspath(Path(__file__).resolve().parent / "style.qss")
    with open(path, 'r') as f:
        style = f.read()
        app.setStyleSheet(style)

    widget = App()
    widget.setWindowTitle("ZEMAS")

    widget.show()
    sys.exit(app.exec())
