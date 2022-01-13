# This Python file uses the following encoding: utf-8
from microscope import Microscope
from PyQt6.QtCore import QThread
class Worker:
    def __init__(self,QThread,Microscope):
        self.scope = Microscope
