# This Python file uses the following encoding: utf-8
import traceback
import numpy as np

from PyQt6.QtCore import QObject, QRunnable, pyqtSignal
from include.image_analysis import *
from include.crystal_text_input import *
from include.CifFile import get_number_with_esd

from scipy.optimize import minimize, differential_evolution, shgo, dual_annealing

class WorkerSignals(QObject):
    finished = pyqtSignal()
    result   = pyqtSignal(tuple) #make sure that you're always returning a tuple
    error    = pyqtSignal(tuple) 

class DiffractionRunnable(QRunnable):

    def __init__(self): 
        super(DiffractionRunnable, self).__init__()
        self.signals = WorkerSignals()
    
    def getRealMatrix(self, rnd=10):
        """returns in the general triclinic case the real space vectors in terms of cartesian axes
        The cartesian axes stuck to the crystal is defined as: a is on x-axis, a-b plane is on x-y plane"""
        a = float(self.crystal_a)
        b = float(self.crystal_b)
        c = float(self.crystal_c)
        alpha = float(self.crystal_alpha) / 360*2*np.pi
        beta  = float(self.crystal_beta)  / 360*2*np.pi
        gamma = float(self.crystal_gamma) / 360*2*np.pi
        av = np.array([a, 0, 0])
        bv = np.array([np.round(b*np.cos(gamma), rnd), np.round(b*np.sin(gamma), rnd), 0])
        c1 = np.round(c*np.cos(beta), rnd) 
        c2 = (c*b*np.cos(alpha) - c1*bv[0])/(bv[1])
        c3 = np.round(np.sqrt(c**2-c1**2-c2**2), rnd)
        cv = np.array([c1, c2, c3])
        matR = np.array([av, bv, cv]).T
        return matR
        
    def getReciprocalMatrix(self, rnd=10):
        """Matrix containing the basis vectors of the recyprocal lattice. Based on the real matrix."""
        a1 = self.RM[:, 0]
        a2 = self.RM[:, 1]
        a3 = self.RM[:, 2]
        vol = np.dot(a1, np.cross(a2, a3))
        astar = np.cross(a2, a3)/vol
        bstar = np.cross(a3, a1)/vol
        cstar = np.cross(a1, a2)/vol
        return np.array([astar, bstar, cstar]).T

    def millerToCartesian(self, vec, typ = "real"):
        """This function returns the coordinates in the cartesian coordinate system stuck to the crystal given a set of miller indices as columns or an array. Standard it will be assumed miller indices as defined by the real coordinate system, but recyp is also valid and must be supplied as typ = 'recyp' as argument"""
        if typ=="real":
            return np.dot(np.array(self.RM), np.array(vec))
        elif typ=="recip" or typ=="reciprocal":
            return np.dot(np.array(self.RRM), np.array(vec))
        else:
            return None
    
    def XR(self, d, rnd=10):
        d=float(d)/360*2*np.pi
        mat=np.array([[1, 0, 0],[0, np.round(np.cos(d), rnd), np.round(-np.sin(d), rnd)], [0, np.round(np.sin(d), rnd), np.round(np.cos(d), rnd)]])
        return mat
    
    def YR(self, d, rnd=10):
        d=float(d)/360*2*np.pi
        mat=np.array([[np.round(np.cos(d), rnd), 0, np.round(np.sin(d), rnd)], [0, 1, 0], [np.round(-np.sin(d), rnd), 0, np.round(np.cos(d), rnd)]])
        return mat

    def ZR(self, d, rnd=10):
        d=float(d)/360*2*np.pi
        mat=np.array([[np.round(np.cos(d), rnd), np.round(-np.sin(d), rnd), 0], [np.round(np.sin(d), rnd), np.round(np.cos(d), rnd), 0], [0, 0, 1]])
        return mat

    def detectorToAbs(self, vec):
        """Detector coordinates are transformed to absolute coordinates."""
        return np.dot(self.ZR(self.diffrot), vec)

    def absToStage(self, vecs):
        alpha=float(self.alpha)
        beta=float(self.beta)
        
        xx = self.XR(-alpha)
        yy = self.YR(-beta)
        
        #The inverse of the absolute to stage
        rot = np.dot(yy, xx)
        return np.dot(np.array(rot), np.array(vecs))

class SpotDiffractionImageAnalyzer(DiffractionRunnable):
    def __init__(self, image, diffrot, alpha, beta, cif_data): 

        super(SpotDiffractionImageAnalyzer, self).__init__()

        # Inputs
        self.image    = image
        self.diffrot  = diffrot
        self.alpha    = alpha
        self.beta     = beta
        self.cif_data = cif_data
        
        # Outputs
        self.hkl   = None
        self.uvw   = None
        self.ormat = None
        

    def run(self):
        try:
            if self.cif_data is None:
                raise ValueError("Structure needs to be loaded from CIF file in Configuration tab.")
            
            lattice_params = self.cif_data.lattice_parameters()
        
            self.crystal_a = lattice_params[0]
            self.crystal_b = lattice_params[1]
            self.crystal_c = lattice_params[2]

            self.crystal_alpha = lattice_params[3]
            self.crystal_beta  = lattice_params[4]
            self.crystal_gamma = lattice_params[5]

            self.RM  = self.getRealMatrix()
            self.RRM = self.getReciprocalMatrix()

            block = self.cif_data.structure_block

            try: 
                SpaceGroupITNumber, _ = get_number_with_esd(block["_space_group_IT_number"])
            except:
                SpaceGroupITNumber, _ = get_number_with_esd(block["_symmetry_Int_Tables_number"]) 

            HB_BravaisLattice = self.SpaceGroupItNumberToBravaiLattice(SpaceGroupITNumber)

            props = GetDiffractionSpotProperties(self.image) #from image_analysis.py
            df = props[0]
            annotated_im = props[1]
            df = GetAngleAndHKLDistances(df) #from image_analysis.py
            df = df[df['Angle']>10] # threshold for angle > 10
            
            row_index_of_pair_of_spots = 0 

            GV1 = df['Vector 1 1/dhkl'].values[row_index_of_pair_of_spots]
            GV2 = df['Vector 2 1/dhkl'].values[row_index_of_pair_of_spots] 
            Ratio_of_two_g_vectors  = GV2 / GV1
            Acute_Angle_Between_Vectors = df['Angle'].values[row_index_of_pair_of_spots]

            # hard-coded inputs for call_crystal
            Camera_Constant_Known = 'N'
            Accuracy = 0.03
            Max_hkl_index = 5
            Crystal_System_Type = 1 
            Avoid_Crystallographically_Equivalent_Solutions = False


            #call function crystal_text_input.py to get 
            results = call_crystal(HB_BravaisLattice, Crystal_System_Type,\
            self.crystal_a, self.crystal_b, self.crystal_c,\
            self.crystal_alpha, self.crystal_beta, self.crystal_gamma,\
            Camera_Constant_Known, Ratio_of_two_g_vectors, Acute_Angle_Between_Vectors, Accuracy,\
            Max_hkl_index, Avoid_Crystallographically_Equivalent_Solutions)[0] 
            self.hkl   = np.array([results[0],results[1],results[2]])
            self.uvw   = np.array([results[8],results[9],results[10]])

            self.ormat = self.calc_ormat(self.uvw, self.hkl, self.diffrot)

        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

        else:
            self.signals.result.emit((self.hkl,self.uvw,self.ormat))

        finally:
            self.signals.finished.emit()

    def calc_ormat(self, za, ref, ang, acur = 1e-9):
        #turn angle from degrees to radians
        ang = ang/360*2*np.pi
        
        #calculate the cartesian equivalents of the vectors
        zaC  = self.millerToCartesian(za)
        refC = self.millerToCartesian(ref, typ = "recip")
        #normalize the vectors
        zaC  = zaC/np.linalg.norm(zaC)
        refC = refC/np.linalg.norm(refC)
        depC = np.cross(zaC, refC)
        print("Constructing mat1")
        #the vectors of the crystal to be transformed
        mat1 = np.array([zaC, refC, depC])
        
        print("Getting mat2")
        #the matrix of corresponding detector vectors
        mat2 = self.ZR(ang)
        
        print("To absolute coordinates...")
        realcords = self.detectorToAbs(mat2) #change to absolute coordinates
        
        print("To stage coordinates...")
        #stagecoords = self.stage.absToStage(realcords)
        stagecoords = self.absToStage(realcords)
        
        print("Calculating ormat")
        #the rotation matrix needs to turn mat 1 (cartesian vectors stuck to crystal) into stagecoords (stage vectors). Therefore
        ormat = np.dot(stagecoords, mat1)
        #multiplying by ormat goes from crystal cartesian vector to stage coordinates, ormat.T (inverse) goes from stage to cartesian.
        return ormat

    def SpaceGroupItNumberToBravaiLattice(self, SpaceGroupITNumber):
        if 1 <= SpaceGroupITNumber <= 2:
            # BravaisLattice = 'Triclinic'
            HB_BravaisLattice = 6
        if 3 <= SpaceGroupITNumber <= 15:
            # BravaisLattice = 'Monoclinic'
            HB_BravaisLattice = 5
        if 16 <= SpaceGroupITNumber <= 74:
            # BravaisLattice = 'Orthorhombic'
            HB_BravaisLattice = 3
        if 75 <= SpaceGroupITNumber <= 142:
            # BravaisLattice = 'Tetragonal'
            HB_BravaisLattice = 2
        if 143 <= SpaceGroupITNumber <= 167:
            # BravaisLattice = 'Trigonal'
            HB_BravaisLattice = 4
        if 168 <= SpaceGroupITNumber <= 194:
            # BravaisLattice = 'Hexagonal'
            HB_BravaisLattice = 4
        if 195 <= SpaceGroupITNumber <= 230:
            # BravaisLattice = 'Cubic'
            HB_BravaisLattice = 1
        return HB_BravaisLattice
        
class AlphaBetaTilter(DiffractionRunnable):
    def __init__(self, alpha, beta, uvw1, uvw2, ormat, diffrot, cif_data): 
        super(AlphaBetaTilter, self).__init__()
        
        self.alpha    = alpha
        self.beta     = beta
        self.uvw1     = uvw1
        self.uvw2     = uvw2
        self.ormat    = ormat
        self.diffrot = diffrot
        self.cif_data = cif_data

    def calcAlphaBeta(self, ormat, uvw1, uvw2, rnd = 10):
        vec1 = np.dot(ormat, self.millerToCartesian(uvw1, typ = "reciprocal"))
        vec2 = np.dot(ormat, self.millerToCartesian(uvw2, typ = "reciprocal"))

        vec1 = vec1/np.linalg.norm(vec1)
        vec2 = vec2/np.linalg.norm(vec2)
        print(vec1,vec2)

        def equation(x): 
            a,b = x
            A = np.dot(self.XR(a), self.YR(b))
            vec2_comp = np.dot(A,vec1)
            vec2_comp = vec2_comp/np.linalg.norm(vec2_comp)
            return 1-np.dot(vec2,vec2_comp)
        
        result = differential_evolution(equation, bounds = ((0,360),(0,360)))
        a,b    = result.x
        y      = np.arccos(np.clip(-result.fun+1,-1,1)) * 360 / (2 * np.pi)

        if y < 1: #within a degree of the requested zone-axis.
            return a,b
        else:
            raise ValueError("New zone axis not possible with double-tilt stage.")

    def run(self):
        try:
            if self.cif_data is None:
                raise ValueError("Structure needs to be loaded from CIF file in Configuration tab.")
                
            if self.ormat is None:
                raise ValueError("Orientation matrix must be calculated.")
            lattice_params = self.cif_data.lattice_parameters()
        
            self.crystal_a = lattice_params[0]
            self.crystal_b = lattice_params[1]
            self.crystal_c = lattice_params[2]

            self.crystal_alpha = lattice_params[3]
            self.crystal_beta  = lattice_params[4]
            self.crystal_gamma = lattice_params[5]

            self.RM  = self.getRealMatrix()
            self.RRM = self.getReciprocalMatrix()

            alpha,beta = self.calcAlphaBeta(self.ormat, self.uvw1, self.uvw2)
            self.alpha = (((self.alpha + alpha) % 360) + 360) % 360
            self.beta  = (((self.beta  + beta ) % 360) + 360) % 360
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

        else:
            self.signals.result.emit((self.alpha,self.beta))

        finally:
            self.signals.finished.emit()
        
