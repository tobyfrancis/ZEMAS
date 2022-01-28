# This Python file uses the following encoding: utf-8
import traceback
import numpy as np

from PyQt6.QtCore import QObject, QRunnable, pyqtSignal
from include.image_analysis import *
from include.crystal_text_input import *
from include.tem_calcs import *

class WorkerSignals(QObject):
    finished = pyqtSignal()
    error = pyqtSignal(tuple)

class ImageAnalyzer(QRunnable):
    def __init__(self):
        super(ImageAnalyzer, self).__init__()

        self.image = None
        self.structure = Structure("arbitrary")

        self.alpha = 90
        self.beta  = 90
        self.gamma = 90

        self.ang = 0
        self.detector = "MSC"
        self.mode = "d"
        self.setting = 30

        self.hkl   = None
        self.uvw   = None
        self.ormat = None
        
        self.signals = WorkerSignals()
        

    def set_image(self,image):
        self.image = image

    def set_abg(self,alpha,beta,gamma): #set angles for analysis
        self.alpha = alpha 
        self.beta  = beta
        self.gamma = gamma

    def run(self):
        try:
            props = GetDiffractionSpotProperties(self.image) #from image_analysis
            df = props[0]
            annotated_im = props[1]
            df = GetAngleAndHKLDistances(df)

            df = df[df['Angle']>10] # threshold for angle > 10
            # these are all inputs for image_analysis
            Avoid_Crystallographically_Equivalent_Solutions = False
            row_index_of_pair_of_spots = 0 # arbitrary choice

            GV1 = df['Vector 1 1/dhkl'].values[row_index_of_pair_of_spots]
            GV2 = df['Vector 2 1/dhkl'].values[row_index_of_pair_of_spots] 
            RR  = GV2 / GV1
            AAA = df['Angle'].values[row_index_of_pair_of_spots]
            Ratio_of_two_g_vectors = RR 
            Acute_Angle_Between_Vectors = AAA 

            Camera_Constant_Known, Accuracy, Max_hkl_index = self.variables_for_HB_code()
            Crystal_System = 6# HB_BravaisLattice It's not reading in the variable name from the other if statement
            Crystal_System_Type = 1  #This needs to have other possibilities - Read it in from CIF file
            Ratio_of_two_g_vectors=RR 
            Acute_Angle_Between_Vectors = AAA 
            answers = call_crystal(Crystal_System, Crystal_System_Type, AP, BP, CP, self.alpha, self.beta, self.gamma, Camera_Constant_Known, Ratio_of_two_g_vectors, Acute_Angle_Between_Vectors, Accuracy, Max_hkl_index, Avoid_Crystallographically_Equivalent_Solutions)[0]
            
            h = answers[0]
            k = answers[1]
            l = answers[2]

            u = answers[8]
            v = answers[9]
            w = answers[10]

            self.hkl = np.array([h,k,l])
            self.uvw = np.array([u,v,w])
            #self.ormat = self.calcOrient(miller(u,v,w), miller(h,k,l),\
            #                             self.ang, self.detector, self.mode, self.setting)
            self.signals.finished.emit()

        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

    def variables_for_HB_code(self): #these  are currently hard-coded variables

        #Get these from microscope metadata
        Camera_Constant_Known= 'N'

        #optional choices
        Accuracy = 0.03
        Max_hkl_index=3
        return Camera_Constant_Known, Accuracy, Max_hkl_index


    def calc_ormat(self, za, ref, ang, detector, mode, setting, acur = 1e-9):
        """The crystal has a certain orientation with respect to the stage. The orientation is most easily found when studying the zone axis (real space) || screen z-axis and a visible reflection on the detector (recyprocal space) defined by an angle from the detector x-axis.
        The orientation must be simply a rotation matrix and hence defines how the cartesian system stuck to the crystal is rotated compared to the stage coordinate system."""
        #first check that za (real space) and ref (recyprocal space) are indeed perpendicular. This follows the normal h*u + k*v + l*w = 0 relationship valid for any crystal system.
        if abs(np.dot(za, ref))<acur:
            #turn angle from degrees to radians
            ang = ang/360*2*np.pi
            
            #calculate the cartesian equivalents of the vectors
            zaC = self.structure.millerToCartesian(za)
            refC = self.structure.millerToCartesian(ref, typ = "recyp")
            #normalize the vectors
            zaC = zaC/np.linalg.norm(zaC)
            refC = refC/np.linalg.norm(refC)
            depC = np.cross(zaC, refC)
            #the vectors of the crystal to be transformed
            mat1 = np.array([zaC, refC, depC]).T
            
            #the matrix of corresponding detector vectors
            c1 = np.array([0,0,1])
            c2 = np.array([np.cos(ang), np.sin(ang), 0])
            c3 = np.array([np.cos(ang+np.pi/2), np.sin(ang+np.pi/2), 0])
            mat2 = np.array([c1, c2, c3]).T
            

            #TODO: get this stuff from configuration

            dec=Detector(detector)  #Added
            if mode == 'd':
                calib_file = 'rotationCalibrationMSC-DIFF.txt'
            if mode == 'i':
                calib_file = 'rotationCalibrationMSC-IMG.txt'
            dec.setCalibration(calib_file, mode = mode) #Added  - change to input filename through load session
            #print(type(dec))   #added
            stage=Stage('Doubletilt')  #added  - change to input through load session

            realcords = dec.detectorToAbs(mat2, mode, setting) #change to absolute coordinates
            
            #stagecoords = self.stage.absToStage(realcords)
            stagecoords = stage.absToStage(realcords)
            
            #the rotation matrix needs to turn mat 1 (cartesian vectors stuck to crystal) into stagecoords (stage vectors). Therefore
            ormat = np.dot(stagecoords, np.linalg.inv(mat1))
            #multiplying by ormat goes from crystal cartesian vector to stage coordinates, ormat.T (inverse) goes from stage to cartesian.
            return ormat

        else:
            print("Warning: ZA vector and reflection vector are not perpendicular")
            return np.identity(3)
