from distutils.util import strtobool

import re
import cmath

import math
import numpy as np
import sys
import textwrap
from include import map


#Define Variables that are input in Fortran

#Get these from CIF FILE
Crystal_System = 1
Crystal_System_Type = 1
AP = 5.43
BP = 5.43
CP = 5.43
ALPHA = 90
BETA = 90
GAMMA = 90

#Get these from microscope metadata
Camera_Constant_Known= 'N'

#Get this from image processing
Ratio_of_two_g_vectors= 1.4 # 1.46
Acute_Angle_Between_Vectors = 45 # 44

#optional choices
Accuracy = 0.03
Max_hkl_index=3
Avoid_Crystallographically_Equivalent_Solutions = True


def LOGO(I, AVER):

    ###print(' *************************************************')
    ##print(' PRACTICAL CRYSTALLOGRAPHY  Version {0}'.format(AVER))

    if I == 2:
        
        print('      Calculation of the Metric Tensor    ')
        
    elif I == 4:
        
        print('   Spacings of Planes  ')
        
    elif I == 5:
        
        print(' Angles between Vectors')
        
    elif I == 6:
        
        print(' Conversion of Plane Normals to Directions')
        print('    or Directions to Plane Normals')
        
    elif I == 7:
        
        print('      Coordinate  Transformations     ')
        print('      (Orientation Relationships)     ')
        
    elif I == 8:
        
        print('  Axis-Angle Pairs for Cubic Crystals ')
        
    elif I == 9:
        
        print(' Hexagonal Lattices: 4 Index Notation ')
        print(' Miller-Bravais notation for Planes  ')
        print(' Weber notation for Directions     ')
        
    elif I == 10:
        
        print('   Magnitudes of Real Space Vectors   ')
        
    else:
        #print('  Analysis of Electron Diffraction Patterns')
        one=1
    ##print(' Copyright:                      1993 ')
    #input(' Press Return ...')

    return


#
# MAGNITUDE OF A LATTICE VECTOR
#
def MAGG(H, K, L, G):

   A, B, C = map.TRANS(G, H, K, L)
   M = H * A + K * B + L * C
   M = math.sqrt(M)
   M = 1.0/M if M != 0.0 else float("nan")
   
   return M


#
# Sorts an array arr(1:n) into ascending numerical order, by straight insertion,
# while making the corresponding rearrangement of the array brr(1:n) .
#
def PIKSR2(N, ARR, BRR):

    for J in range(0, N):
        A=ARR[J]
        B=BRR[J]
        for I in range(J - 1, 0, -1):
            if ARR[I] > A: 
                ARR[I+1]=ARR[I]
                BRR[I+1]=BRR[I]
            else:
                break
        I=0
        ARR[I+1]=A
        BRR[I+1]=B

    return (ARR, BRR)


#
#
#
def PUT(I, J):
    
    H, K, L = list_of_numbers(' Components of VECTOR {0} of CRYSTAL {1} ?'.format(I,J), 3)
    JR = REEDI(' Components REAL( =0) or RECIPROCAL( =1) ?')

    return [H, K, L, JR]


#
# Read a real number with error checks
#
def REED(message, min = sys.float_info.min, max = sys.float_info.max):

    while True:
        try:
            userInput = float(input(message))
            if userInput < min or userInput > max:
                #print("Not in range! Try again.")
                continue
        except ValueError:
            #print("Not an integer! Try again.")
            continue
        else:
            return userInput
            break


#
# Read an integer option from the command line with error checks
#
def REEDI(message, min = 0, max = 100):
    
    while True:
        try:
            userInput = int(input(message))
            if userInput < min or userInput > max:
                #print("Not in range! Try again.")
                continue
        except ValueError:
            #print("Not an integer! Try again.")
            continue
        else:
            return userInput
            break


#
# User input a list of numbers of specified type
#
def list_of_numbers(message, entries, type=int, min=0, max=180):
    
    while True:
        try:
            output = input(message)
            elements = re.split('[,\s]+', output)
            processed = []
            for x in elements:
                if x != None and x != '':
                    val = type(x)
                    if val >= min and val <= max:
                        processed.append(val)
                    else:
                        #print("One of your numbers is out of bounds.  Try again.")
                        continue
            if len(processed) < entries or len(processed) > entries:
                #print("You've entered either too many or too few values.  Try again.")
                continue
        except ValueError:
            #print("You tried entering a number of wrong type.  Try again.")
            continue
        else:
            return processed

        
#
#
#
def CRYSTAL(F, G, J2):
    
    #J = REEDI(''' Choose Crystal System for ** CRYSTAL {0} **
 #1 = Cubic            4 = Hexagonal or Trigonal
 #2 = Tetragonal       5 = Monoclinic
 #3 = Orthorhombic     6 = Triclinic
#> '''.format(J2))

    J=Crystal_System

    G = MET(J)

    F = map.MAP_UTIL_INVERS(G)
    return F


#
# Select crystal system and, if not cubic, lattice type as well
#
#  Call this to recover JTYPE even if you do not need to due
#  to the bundling of the return value
#
def QUES(JCHOS, JTYPE=0):

#    message = '''
#Choose Crystal System 
#1 = Cubic          4 = Hexagonal or Trigonal
#2 = Tetragonal     5 = Monoclinic
#3 = Orthorhombic   6 = Triclinic
#> '''
#    J = REEDI(message, 1, 6)

    J=Crystal_System
    result = (J, JTYPE)
    if JCHOS != 1:
#        message = '''Choose Lattice Type
#1 = Primitive
#2 = Body Centered (Cubic, Tetragonal, Orthorhombic)
#3 = Face Centered (Cubic, Orthorhombic)
#4 = A-Centered (Orthorhombic)
#5 = B-Centered (Orthorhombic)
#6 = C-Centered (Orthorhombic or Monoclinic)
#> '''
        #JTYPE = REEDI(message, 1, 6)
        JTYPE = Crystal_System_Type
        #print(JTYPE)
        result = (J, JTYPE)
        
    return result


# 
# Calculate the metric tensor for different crystal structures
#
def MET(J):
    global AP, BP, CP, ALPHA, BETA, GAMMA
    G = np.zeros(9)
    if J == 1: 

        #AP = REED('Lattice Parameter ?')
        
        #print(' CUBIC SYSTEM {0:10.4f} ANGSTROMS'.format(AP))
        G[0] = AP * AP
        G[4] = G[0]
        G[8] = G[0]
          
    if J == 2:
        AP, CP = list_of_numbers(' Lattice Parameters "a" & "c" ?', 2)
        BP=AP
        #print(" TETRAGONAL SYSTEM {0:10.4f} {1:10.4f} ANGSTROMS".format(AP,CP))
        G[0]=AP*AP
        G[4]=G[1]
        G[8]=CP*CP
          
    if J == 3:
        #AP,BP,CP = list_of_numbers(' Lattice Parameters "a", "b", "c" ?', 3)
        #print(" ORTHORHOMBIC SYSTEM {0:10.4f} {1:10.4f} {2:10.4f} ANGSTROMS".format(AP,BP,CP))
        G[0]=AP*AP
        G[4]=BP*BP
        G[8]=CP*CP

    if J == 4:
        AP, CP = list_of_numbers(' Lattice Parameters "a" & "c" ?', 2)
        BP=AP
        #print(" HEXAGONAL SYSTEM {0:10.4f} {1:10.4f} ANGSTROMS".format(AP,CP))
        G[0]=AP*AP
        G[1]=-0.5*G[0]
        G[3]=G[1]
        G[4]=G[0]
        G[8]=CP*CP

    if J == 5:
        #AP,BP,CP = list_of_numbers(' Lattice Parameters "a", "b", "c" ?', 3)
        BETA = list_of_numbers(' Angle Beta ?', 1)
        #print(" MONOCLINIC SYSTEM {0:10.4f} ANGSTROMS {1:10.4f} ANGSTROMS {2:10.4f} ANGSTROMS {3:10.4f} DEGREES".format(AP, BP, CP, BETA))
        G[0] = AP * AP
        G[2] = AP * CP * math.cos(BETA * 2.0 * 3.14159/360.0)
        G[4] = BP*BP
        G[6] = G[2]
        G[8] = CP*CP

    if J == 6:
        #AP, BP, CP = list_of_numbers(' Lattice Parameters "a", "b", "c" ?', 3)
        #ALPHA, BETA, GAMMA = list_of_numbers(' Angles Alpha, Beta, Gamma ?', 3, float)
        #print(" TRICLINIC SYSTEM {0:10.4f} {1:10.4f} {2:10.4f} ANGSTROMS ALPHA={3:8.2f} BETA={4:8.2f} GAMMA={5:8.2} DEGREES".format(AP, BP, CP, ALPHA, BETA, GAMMA))
        G[0] = AP * AP
        G[1] = AP * BP * math.cos(GAMMA * 2.0 * 3.14159 / 360.0)
        G[3] = AP * CP * math.cos(BETA * 2.0 * 3.14159 / 360.0)
        G[3] = G[2]
        G[4] = BP * BP
        G[5] = BP * CP * math.cos(ALPHA * 2.0 * 3.14159 / 360.0)
        G[6] = G[3]
        G[7] = G[6]
        G[8] = CP * CP

    return G


#
# Generic yes-no question inputs
#
def yes_or_no(question, default='no'):
    
    if default is None:
        prompt = " [y/n] "
    elif default == 'yes':
        prompt = " [Y/n] "
    elif default == 'no':
        prompt = " [y/N] "
    else:
        raise ValueError(f"Unknown setting '{default}' for default.")

    while True:
        try:
            resp = input(question + prompt).strip().lower()
            if default is not None and resp == '':
                return default == 'yes'
            else:
                return strtobool(resp)
        except ValueError:
            print("Please respond with 'yes' or 'no' (or 'y' or 'n').")


#
#
#
def dsp_trans(F, H, K, L):
    H1, K1, L1 = map.TRANS(F, H, K, L)
    A=map.MAG(H,K,L,H1,K1,L1)
    A=math.sqrt(A)
    AD=1.0/A
    return AD


#
# Obtains the spacing of planes, when the Miller indices are not those which are
# systematically absent.
#
def DSP(F, H, K, L, JTYPE, JLOG=-1):

    AD = None
    if JTYPE == 1:
        AD = dsp_trans(F, H, K, L)
        return (AD, JLOG)
        
    JLOG = map.TYPE(H, K, L, JTYPE)
    if JLOG == 1:
        return (AD, JLOG)

    AD = dsp_trans(F, H, K, L)
    return (AD, JLOG)


#
# This code has two distinct operations:
# + Perform electron diffraction analysis for cases where the camera constant is known.
#   An input of two d-spacings from a pair of reciprocal lattice vectors and the angle between the two vectors are required.
#   Calculations are made for a set of Miller indices which vary within a specified range. 
#
#
# https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/ed1-b.html
# https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/ed2-b.html
#
def electron_diffraction():

    G = np.empty(9)    #
    F = np.empty(9)    # converts components from reciprocal to real space.
    AH = np.zeros(9)   # Miller indices used for the first vector.
    AK = np.zeros(9)   # Miller indices used for the first vector.
    AL = np.zeros(9)   # Miller indices used for the first vector.
    AAD = np.zeros(9)  # d-spacing for AH(I), AK(I), AL(I).
    AHH = np.zeros(9)  # Miller indices used for the second vector.
    AKK = np.zeros(9)  # Miller indices used for the second vector.
    ALL = np.zeros(9)  # Miller indices used for the second vector.
    ABD = np.zeros(9)  # d-spacing for AHH(I), AKK(I), ALL(I).
    AU1 = np.zeros(9)  # normalised components of the zone axis.
    AV1 = np.zeros(9)  # normalised components of the zone axis.
    AW1 = np.zeros(9)  # normalised components of the zone axis.
    AAX = np.zeros(9)  # angles between the lattice vectors.
    
    LOGO(1, 1.0)
    GV1=1.0
    GV2=1.0
    JYES = 2
    
    while True:
        JI = 0
        if JYES == 2:
            G = np.zeros(9)

           # #print("{:^80}".format('*'*53))
            ##print("** {:^58} **".format('All dimensions in Angstroms, All angles in degrees'))
            ##print("{:^80}".format('If the camera constant is known, input data consist'))
            ##print("{:^80}".format('of d-spacings obtained from two reciprocal lattice'))
            ##print("{:^80}".format('vectors, and the acute angle between these vectors.'))
            ##print("{:^80}".format('If the camera constant is unknown, then the ratio'))
            ##print("{:^80}".format('of the  vectors may be used instead.'))
            ##print("{:^80}".format('*'*53))

            #J40 = yes_or_no('Is the Camera Constant known?')
            J40 = Camera_Constant_Known
            #J50 = yes_or_no('Do you wish to avoid crystallographically equivalent solutions?')
            J50 = Avoid_Crystallographically_Equivalent_Solutions
            ##print('J50', J50)
            message = '''
Specify maximum range of Miller Index. The typical value is 3  (avoid values > 4). 
(Larger values lead to longer run times PHX: use SETJD TIME=5mins to avoid errors)
> '''
            #J21 = REEDI(message)
            J21 = Max_hkl_index
            J21 = J21 + 1
            J20 = 2 * J21 - 1
            J21 = -J21
            J22 = J20 + J21 + 1

            #
            # Select the crystal system and lattice type
            #
            
            J, JTYPE = QUES(0)

            #
            # Collect lattice parameters corresponding to selected crystal system
            # and calculate corresponding metric tensor
            #
            
            G = MET(J)

            #
            # Calculate the inverse of the metric tensor
            #
            
            F = map.MAP_UTIL_INVERS(G)

            #ACCU = REED('Accuracy of measurement? (Typical value 0.03)', 0, 0.08)
            ACCU = Accuracy
        J60 = 1
        AEQ = np.zeros((5,3))

        #
        # If camera constant is known, collect d-spaces of vectors
        # then compute the ratio of the reciprocal lattice vectors
        # Otherwise, collect the ratio of the two reciprocal lattice
        # vectors
        #

        if J40 == True:
        
            GV1 = REED(' d-spacing for first vector ?')
            GV2 = REED(' d-spacing for second vector ?')
            RR = GV1/GV2
        
        else:
        
            #RR = REED('Ratio of two Reciprocal Lattice Vectors?')
            RR = Ratio_of_two_g_vectors

        #
        # Collect the angle between the vectors.
        # 
        #AAA = REED('Acute angle between vectors ?')
        AAA = Acute_Angle_Between_Vectors
        #
        
        # Assure that the ratio between the two reciprocal lattice vectors
        # is greater than one
        #
        try:
            if RR < 1.0:
                RR = 1.0/RR
        except ValueError:
            RR[RR < 1.0] = 1.0/RR[RR < 1.0]

        AA = (AAA * 2.0 * 3.14159) / 360.0
        
        answers = map.analyze(J21, J22, J20, ACCU, F, GV1, GV2, RR, AA, JTYPE, J40, J50, J60, AEQ)
        
        ##print('-----------------------------Electron Diffraction Analysis-------------------------------------')

        message = 'Continue analysis (1) Fresh analysis (2) Exit to main menu (0)?'
        #JYES = REEDI(message, 0, 2)
        JYES = 0  # Added 
        global I # Added 
        I = 11 # Added 
        if JYES == 0:
            break
    
    
    return answers


#
#
#
def MAP_UTIL_TRANS(LR, X1, X2, X3):

    X4 = LR[0] * X1 + LR[1] * X2 + LR[2] * X3
    X5 = LR[3] * X1 + LR[4] * X2 + LR[5] * X3
    X6 = LR[6] * X1 + LR[7] * X2 + LR[8] * X3

    return [X4, X5, X6]


#
# To multiply a single row matrix with a 3x3 martix
#
def MAP_UTIL_TRANS2(D,H,K,L):

    HH = H * D[0] + K * D[3] + L * D[6]
    KK = H * D[1] + K * D[4] + L * D[7]
    LL = H * D[2] + K * D[5] + L * D[8]

    return [HH, KK, LL]


# To convert the components of a vector from real to reciprocal space and vice-versa
# This subroutine works for any crystal system.
# 
# Reference
#  Worked Examples in the Geometry of Crystals
#  by H. K. D. H. Bhadeshia
#  Published by the Institute of Materials, London, 1987.
#  The metric tensor is explained on page 22
#
# https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/convert-b.html 
#
def MAP_CRYSTAL_CONVERT(F, G, H, K, L, JR1):

# JR1 is 0 if vector defined in real space, 1 if in reciprocal space
#
# The metric tensor F converts the components from reciprocal to real space. 
# The metric tensor G converts the components from real to reciprocal
#  space, and G is the inverse of F and vice-versa.

   if JR1 == 0:
       
       U, V, W = MAP_UTIL_TRANS(G, H, K, L)
       
   else:
       
       U, V, W = MAP_UTIL_TRANS(F, H, K, L)

   return [U, V, W]


#
#  Calculate the magnitude of a lattice vector given the metric tensor G
#
def MAP_UTIL_MAGG(H, K, L, G):

    A, B, C = MAP_UTIL_TRANS(G, H, K, L)
    M=H*A+K*B+L*C
    M=math.sqrt(M)
    M=1.0/M if M != 0.0 else float("nan")

    return M


#
#
#
def LIST(D,H,K,L):
    
#C H K L  FROM CRYSTAL 1
    
    HH, KK, LL = MAP_UTIL_TRANS2(D,H,K,L)
    #print('{0} {1} {2} {3} {4} {5}'.format(H,K,L,HH,KK,LL))

    return [HH, KK, LL]


#
#
#
def LIST2(C, D):

    H, K, L = list_of_numbers(' Components of VECTOR ?', 3 )
    JR = REEDI(' Components REAL( =0) or RECIPROCAL( =1) ?')
    J = REEDI(' Which CRYSTAL is the VECTOR from (1 or 2)?')
    if J == 1:
        J3 = 2
    else:
        J3 = 1

    if J == 1 and JR == 0:
        HH, KK, LL = MAP_UTIL_TRANS(C, H, K, L)
    if J == 1 and JR == 1:
        HH, KK, LL = MAP_UTIL_TRANS2(D, H, K, L)
    if J == 2 and JR == 0:
        HH, KK, LL = MAP_UTIL_TRANS(D, H, K, L)
    if J == 2 and JR == 1:
        HH, KK, LL = MAP_UTIL_TRANS2(C, H, K, L)

    if JR == 0:
        print(' Real space VECTOR: CRYSTAL {0} = {1} {2} {3} parallel to real space VECTOR: CRYSTAL {4} = {5} {6} {7}'.format(J,H,K,L,J3,HH,KK,LL))
    else:
        print(' Reciprocal space VECTOR: CRYSTAL {0} = {1} {2} {3}  parallel to reciprocal space VECTOR: CRYSTAL {4}  = {5} {6} {7}'.format(J, H, K, L, J3, HH, KK, LL))

    return


#
#
#
def LIST3(C, U, V, W):

# C H K L  FROM CRYSTAL 1
    
    UU, VV, WW = MAP_UTIL_TRANS(C, U, V, W)
    #print('{0} {1} {2} {3} {4} {5}'.format(U, V, W, UU, VV, WW))
    
    return [UU, VV, WW]


#
# 
#
def CORD():

    G1 = np.zeros(9)
    G2 = np.zeros(9)
    F1 = np.zeros(9)
    F2 = np.zeros(9)
    A = np.zeros(9)
    B = np.zeros(9)
    
    LOGO(7,1.0)
    JYES = 1

    while True:

        if JYES == 1:

            for I in range(0,9):
              G1[I]=0.0
              G2[I]=0.0

            print(''' Input data consist of four vectors (real or reciprocal):
 (Vector 1 of Crystal 1) parallel to (Vector 1 of Crystal 2)
 (Vector 2 of Crystal 1) parallel to (Vector 2 of Crystal 2)''')

            F1 = CRYSTAL(F1,G1,1)
            F2 = CRYSTAL(F2,G2,2)

        H1, K1, L1, JR1 = PUT(1, 1)
        H2, K2, L2, JR2 = PUT(2, 1)
        H4, K4, L4, JR4 = PUT(1, 2)
        H5, K5, L5, JR5 = PUT(2, 2)

        #print(' CRYSTAL 1 {0} {1} {2} parallel to {3} {4} {5} CRYSTAL 2'.format(H1,K1,L1,H4,K4,L4))
        #print(' CRYSTAL 1 {0} {1} {2} parallel to {3} {4} {5} CRYSTAL 2'.format(H2,K2,L2,H5,K5,L5))

        H3, K3, L3, JR3 = MAP_UTIL_CROSS(H1,K1,L1,JR1,H2,K2,L2,JR2,F1,G1)
        H6, K6, L6, JR6 = MAP_UTIL_CROSS(H4,K4,L4,JR4,H5,K5,L5,JR5,F2,G2)

        M1 = MAP_UTIL_MAGG(H1, K1, L1, G1)
        M2 = MAP_UTIL_MAGG(H2, K2, L2, G1)
        M3 = MAP_UTIL_MAGG(H3, K3, L3, G1)
        M4 = MAP_UTIL_MAGG(H4, K4, L4, G2)
        M5 = MAP_UTIL_MAGG(H5, K5, L5, G2)
        M6 = MAP_UTIL_MAGG(H6, K6, L6, G2)

        KP=M1/M4
        GP=M2/M5
        MP=M3/M6

        A[0]=H1*KP
        A[3]=K1*KP
        A[6]=L1*KP
        A[1]=H2*GP
        A[4]=K2*GP
        A[7]=L2*GP
        A[2]=H3*MP
        A[5]=K3*MP
        A[8]=L3*MP
        B[0]=H4
        B[1]=H5
        B[2]=H6
        B[3]=K4
        B[4]=K5
        B[5]=K6
        B[6]=L4
        B[7]=L5
        B[8]=L6

        C = map.MAP_UTIL_INVERS(B)
        D = map.MAP_UTIL_PROD(A,C)

        print('''** COORDINATE TRANSFORMATION MATRIX (1 J 2) **
    {0} {1} {2} {3} {4} {5} {6} {7} {8}'''.format(D[0], D[1], D[2], D[3], D[4], D[5], D[6], D[7], D[8]))

        C = map.MAP_UTIL_INVERS(D)

        print('''**   INVERSE TRANSFORMATION MATRIX  (2 J 1) **
    {0} {1} {2} {3} {4} {5} {6} {7} {8}'''.format(C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8]))

        JYES2 = yes_or_no(' List parallel vectors?')
        if JYES2 == True:
            I4 = REEDI(' Specify maximum value of Miller Index (e.g., 3)')
            I5 = 2 * I4 + 1
            #print(' (PLANE: CRYSTAL 1)           (PLANE: CRYSTAL 2)')
            for I1 in range(1, I5):
                H = I1 - 1 - I4
                for I2 in range(1, I5):
                    K=I2-1-I4
                    for I3 in range(1, I5):
                        L=I3-1-I4
                        if H == 0 and K == 0 and L == 0: #)GOTO74
                            break
                        HH, KK, LL = LIST(D, H, K, L)

            #print('  [DIRECTION: CRYSTAL 1]  [DIRECTION: CRYSTAL 2]')
            for I1 in range(1, I5):
                U=I1-1-I4
                for I2 in range(1, I5):
                    V=I2-1-I4
                    for I3 in range(1, I5):
                        W=I3-1-I4
                        if U == 0 and V == 0 and W == 0: #)GOTO 774
                            break
                        UU, VV, WW = LIST3(C, U, V, W)

            while True:
                JYES3 = yes_or_no(' Study any particular vector in detail?')
                if JYES3 == True:
                    
                    LIST2(C,D)
                    
                else:
                    
                    break

        JYES = REEDI(''' Continue with same crystals  = 0
 Start Fresh Analysis         = 1
 Return to main menu          = 2
 > ''', 0, 2)
            
        if JYES == 2:
            break

    return


# Multiply a 3x3 matrix with a single column matrix
#
#
# To calculate the coordinate transformation matrix relating 
#   two crystals of arbitrary structure, from data consisting
#   of a pair of vectors (real or reciprocal) from each crystal.
#   The pair from one crystal has to be parallel to that from 
#   the other crystal. This is useful in studying the orientation
#   relationships between crystals.
#
# This subroutine works for any crystal system.
# 
# Reference
#  Worked Examples in the Geometry of Crystals
#  by H. K. D. H. Bhadeshia
#  Published by the Institute of Materials, London, 1987.
#  See pages 12-18 including Examples 4 and 5.
# https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/cord-b.html
#
def MAP_CRYSTAL_CORD(H1,K1,L1,JR1,H2,K2,L2,JR2,H4,K4,L4,JR4,H5,K5,L5,JR5,G1,G2):
#
# F1 and F2 are described as F, and G1 and G2 as G in what follows.
# The metric tensor F converts the components from reciprocal to
#  real space. 
# The metric tensor G converts the components from real to reciprocal
#  space, and G is the inverse of F and vice-versa.
# 
      F1 = map.MAP_UTIL_INVERS(G1)
      F2 = map.MAP_UTIL_INVERS(G2)

# Find vector [H3 K3 L3] which is normal to [H1 K1 L1] and [H2 K2 L2]

      H3, K3, L3 = MAP_UTIL_CROSS(H1,K1,L1,JR1,H2,K2,L2,JR2,F1,G1)
      
# Find vector [H6 K6 L6] which is normal to [H4 K4 L4] and [H5 K5 L5]

      H4, K4, L4 = MAP_UTIL_CROSS(H4,K4,L4,JR4,H5,K5,L5,JR5,F2,G2)

# Find the magnitude of each of the six vectors

      M1 = MAP_UTIL_MAGG(H1,K1,L1,G1)
      M2 = MAP_UTIL_MAGG(H2,K2,L2,G1)
      M3 = MAP_UTIL_MAGG(H3,K3,L3,G1)
      M4 = MAP_UTIL_MAGG(H4,K4,L4,G2)
      M5 = MAP_UTIL_MAGG(H5,K5,L5,G2)
      M6 = MAP_UTIL_MAGG(H6,K6,L6,G2)

# Hence calculate the normalisation constants
#  see for example, the constants k, g, m on page 17

      KP=M1/M4
      GP=M2/M5
      MP=M3/M6
#
# Calculate the matrices which are manipulated to enable the
#   coordinate transformation matrix D, and its inverse C
#
      A[0]=H1*KP
      A[3]=K1*KP
      A[6]=L1*KP
      A[1]=H2*GP
      A[4]=K2*GP
      A[7]=L2*GP
      A[2]=H3*MP
      A[5]=K3*MP
      A[8]=L3*MP

      B[0]=H4
      B[1]=H5
      B[2]=H6
      B[3]=K4
      B[4]=K5
      B[5]=K6
      B[6]=L4
      B[7]=L5
      B[8]=L6

      C = map.MAP_UTIL_INVERS(B)
      D = map.MAP_UTIL_PROD(A,C)
      C = map.MAP_UTIL_INVERS(D)

      return [D, C]

#
# Calculate the determinant of a 3 x 3 matrix
#
def MAP_UTIL_DET(G):

    DEL=G[0]*(G[4]*G[8]-G[5]*G[7])+G[1]*(G[5]*G[6]-G[3]*G[8])+G[2]*(G[3]*G[7]-G[6]*G[4])

    return DEL

  
# Takes a (vector) cross product of vectors [H1 K1 L1] and [H2 K2 L2]
#   to generate vector [H3 K3 L3] 
# Converts all vectors to real space
#
# 
# Reference
#  Worked Examples in the Geometry of Crystals
#  by H. K. D. H. Bhadeshia
#  Published by the Institute of Materials, London, 1987.
#  See pages 23-24 including equation 10b for details.
#
def MAP_UTIL_CROSS(H1, K1, L1, JR1, H2, K2, L2, JR2, F, G):

# JR1 and JR2 have values of 0 or 1 if real or reciprocal
#   lattice vectors respectively
#
# The metric tensor F converts the components from reciprocal to
#  real space. 
# The metric tensor G converts the components from real to reciprocal
#  space, and G is the inverse of F and vice-versa.
#
# Calculate volume of cell from the determinant of G

    H3 = 0
    L3 = 0
    K3 = 0
    JR3 = 0

    DEL = MAP_UTIL_DET(G)
    VOL = math.sqrt(DEL)
#
# Refer all vectors to real space coordinates if not already so.
#
    if JR1 == 1:

        A, B, C = MAP_UTIL_TRANS(F, H1, K1, L1)
        H1=A
        K1=B
        L1=C
        
    else:

        if JR2 == 1: #GOTO 3

            A, B, C = MAP_UTIL_TRANS(F, H2, K2, L2)
            H2=A
            K2=B
            L2=C
            
        else:
#
# Compute the vector resulting from cross product

            H3=(K1*L2-L1*K2)*VOL
            K3=(L1*H2-H1*L2)*VOL
            L3=(H1*K2-H2*K1)*VOL
#
# Convert components of resultant vector to real space coordinates
#
            A, B, C = MAP_UTIL_TRANS(F, H3, K3, L3)
            H3=A
            K3=B
            L3=C
#
# Reset basis indicators to imply that all vectors are referred
#  to real space coordinates
#
            JR2=0
            JR1=0
            JR3=0

    return [H3, K3, L3, JR3]


#
# Function which calculates the scalar product of two 3-d vectors.
#
def MAP_UTIL_MAG(H, K, L, H1, K1, L1):

    result = H*H1 + K*K1 + L*L1

    return result


# To obtain the spacing of planes, when the Miller indices are not
#  those which are systematically absent
#
def MAP_CRYSTAL_DSP(F, H, K, L, JTYPE):
#
# AD is the spacing of the planes with Miller indices (HKL)
# The metric tensor F converts the components from reciprocal to
#  real space. 
# The metric tensor G converts the components from real to reciprocal
#  space, and G is the inverse of F and vice-versa.
#
# JTYPE represents the lattice type, with 1=primitive etc.
#   For further details on JTYPE see subroutine TYPE
# Note that there are no systematic absences for a primitive lattice
#   so that the calculation can jump straight away to statement 6.
# For all other lattice types, subroutine TYPE determines whether 
#   there is a systematic absence (i.e. JLOG=1) and the spacing is
#   calculated only if there is no absence

    if JTYPE != 1: 
        JLOG = MAP_CRYSTAL_TYPE(H, K, L)
        if JLOG == 1:
            return
# Calcuate spacing using a scalar product of the real and reciprocal
#   lattice components of the same vector
    else:
# 10
        MAP_UTIL_TRANS(F, H, K, L, H1, K1, L1)
        A=MAP_UTIL_MAG(H, K, L, H1, K1, L1)
        A=math.sqrt(A)
        AD=1.0/A
# 20
    return (AD, JLOG)


#
# 
#
def MAP_UTIL_ODD(A):

    IH2=A/2
    IH3=(A+1)/2
    if IH2 == IH3:
        JODD=0
    else:
        JODD=1
    return JODD


# This helps to determine systematic absences due to lattice type
#
# JTYPE defines the lattice type (face-centered cubic etc.)
#       JTYPE=1 is for primitive
#       JTYPE=2 is for body-centered
#       JTYPE=3 is for face-centered
#       JTYPE=4 is for A-centered
#       JTYPE=5 is for B-centered
#       JTYPE=6 is for C-centered
# Reference: "Essentials of Crystallography"
#            by D. McKie and C. McKie, Blackwell Scientific Publ.
#            London, 1986, page 197
#
#
def MAP_CRYSTAL_TYPE(H, K, L, JTYPE, JLOG):

    if JTYPE == 1:
        JLOG=0
    elif JTYPE == 2:

        A=H+K+L
        JODD = MAP_UTIL_ODD(A)
        if JODD == 1:
            JLOG=1

    elif JTYPE == 3:

        J1 = MAP_UTIL_ODD(H)
        J2 = MAP_UTIL_ODD(K)
        J3 = MAP_UTIL_ODD(L)
        if (J1 == 0 and J2 == 0 and J3 == 0) or (J1 == 1 and J2 == 1 and J3 == 1):
            JLOG = 0
        else:
            JLOG = 1
        
    elif JTYPE == 4:

        A=K+L
        JODD = MAP_UTIL_ODD(A)
        if JODD == 1:
           JLOG=1
        else:
            JLOG=0
        
    elif JTYPE == 5: #GOTO 6

        A=H+L
        JODD = MAP_UTIL_ODD(A)
        if JODD == 1:
            JLOG = 1
        else:
            JLOG = 0
        
    elif JTYPE == 6: #GOTO 8

        A=H+K
        JODD = MAP_UTIL_ODD(A)
        if JODD == 1:
            JLOG=1
        else:
            JLOG=0
#
# Now check for systematic absences and set JLOG=0 if the HKL is not absent, 1 if absent
#
    return JLOG


#
# JR1 HAS VALUE OF 0 IF VECTORS ARE REAL, 1 IF RECIPROCAL
#
def CONVERT(F,G,H,K,L,JR1):

    if JR1 == 0:
        U, V, W = map.TRANS(G, H, K, L)
    else:
        U, V, W = map.TRANS(F,H,K,L)
    
    return (U, V, W)


#
# 
#
def MAP_UTIL_NORM(X,Y,Z):

   C=math.sqrt(X*X + Y*Y + Z*Z)
   X=X/C
   Y=Y/C
   Z=Z/C

   return [ X, Y, Z ]


#
# Main menu
#
def logo3(AVER):

    message = '''
*************************************************
 CAMBRIDGE UNIVERSITY
 PRACTICAL CRYSTALLOGRAPHY Version {0}
 1.  Analysis of Electron Diffraction
 2.  Calculation of Metric Tensor
 3.  Interplanar Spacings
 4.  Angles between Vectors
 5.  Convert between Real & Reciprocal
 6.  Axis-Angle Pairs
 7.  Orientation Relations
 8.  Documentation           10. Vector Size
 9.  Four Index Notation     11. Exit
         by  H. K. D. H. Bhadeshia 
 Choose program number .......       '''.format(AVER)
    
    #I = REEDI(message, 1, 11)
    I = 1
    return I


#
# An array LST(I,J) is sorted so that the second column is ordered to give LST(IVAR,J) in ascending order, where IVAR is the dimension of the first column. 
#
def MAP_UTIL_SORT2(IVAR, IVAL, LST):

    ASTORE = np.zeros(100)
    for J in range(0, IVAL):
        for I5 in range(0, IVAR):
            ASTORE[I5]=LST[I5, J]
        for I in range(J-1, 0, -1):
            if LST[IVAR, I] > ASTORE[ivar]:
                for I3 in range(0, IVAR):
                    LST[I3, I+1] = LST[I3, I]
            else:
                break
        I=0
        for I4 in range(0, IVAR):
            LST[I4, I+1] = ASTORE[I4]

    return LST


#
# To convert from four index notation to three index notation for a hexagonal lattice.
# MAP: https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/notat1-b.html
#
#  MAP_CRYSTAL_NOTAT1 converts from four index notation to three index notation for the hexagonal lattice system.
# If the components in four index notation are U,V,T,W, then the components in three index notation are given by:-
#
# U1 = U - T, U2 = V - T, U3 = W, in real space, or
#
# U1 = U, U2 = V, U3 = W, in reciprocal space.
#
def MAP_CRYSTAL_NOTAT1(U, V, T, W, JR2):

    if JR2 != 1: #GOTO 5
        U1=U-T
        U2=V-T
        U3=W
    else:

        U1=U
        U2=V
        U3=W

    return [U1, U2, U3]    

      
# To convert from three to four index notation
#  for the hexagonal system
#
# JR=0 if vector is in real space, JR=1 if in reciprocal space.
#
# Reference:
#  Chapter on Crystallography, by H. K. D. H. Bhadeshia, 
#   in "Microstructural Characterisation" edited by E. Metcalfe,
#   Institute of Metals, London, 1988, Appendix 2, page 30.
#
def MAP_CRYSTAL_NOTAT2(U, V, W, JR2):

    if JR2 == 1: 
        U1=(2.0*U-V)/3.0
        U2=(2.0*V-U)/3.0
        U3=W
        U4=-1.0*(U+V)/3.0
    else:
        
        U1=U
        U2=V
        U3=W
        U4=-1.0*(U+V)

    return [U1, U2, U3, U4]


#
# METRIC
#
def calculate_metric_tensor():

    LOGO(2, 1.0)

    while True:

        G = np.zeros(9)

        print('''*****************************************************
  All dimensions in Angstroms
  All angles in Degrees
  Input data define the crystal system 
*****************************************************''')

        #
        # Select the crystal system
        #
            
        J, JTYPE = QUES(1)
        
        #
        # Collect lattice parameters corresponding to selected crystal system
        # and calculate corresponding metric tensor
        #
            
        G = MET(J)

        #
        # Calculate the inverse of the metric tensor
        #
            
        F = map.MAP_UTIL_INVERS(G)

        #
        # Calculate the square root of the determinant of the metric tensor
        #
        
        DEL = MAP_UTIL_DET(G)
        DEL = math.sqrt(DEL)

        #print('*********************************************')

        #print(F[0],F[1],F[2],F[3],F[4],F[5],F[6],F[7],F[8],G[0],G[1],G[2],G[3],G[4],G[5],G[6],G[7],G[8],DEL)
        JYES = yes_or_no('Continue Analysis?')
        if JYES == False:
            break
         
    return


#
# Computes the Metric Tensor which is in effect a coordinate transformation relating
# the real and reciprocal axes.
#
def DSPACE():
    AAA = np.zeros(3000)
    ARR = np.zeros(400)
    BRR = np.zeros(400)
    LOGO(4, 1.0)
    while True:

        G = np.zeros(9)
        JK=0
        II1 = REEDI(' Specify maximum value of Miller index (1-14)', 1, 14)
        II1 = II1 + 1

        #
        # Select the crystal system and lattice type
        #
        
        J, JTYPE = QUES(0)
        
        #
        # Collect lattice parameters corresponding to selected crystal system
        # and calculate corresponding metric tensor
        #
            
        G = MET(J)
        
        #
        # Calculate the inverse of the metric tensor
        #

        F = map.MAP_UTIL_INVERS(G)

        II2=0
        
        for I1 in range(0, II1):
            H=I1-1
            for I2 in range(1, II1):
                K = I2 - 1
                for I3 in range(0, II1):
                    L=I3-1
                    cond = H == 0 and K == 0 and L == 0
                    if not cond:
                        AD, JLOG = DSP(F, H, K, L, JTYPE)
                        #print(AD, JLOG)
                        if JLOG == 1:
                            break
                        else:
                            II2 = II2 + 1
                            AAA[II2] = AD
                            jump_out = False
                            for JI in range(0, II2 - 1):
                                if AD == AAA[JI]:
                                    jump_out = True
                                    break
                                
                            if jump_out == True:
                                break
                            
                            IN = int(H*H + K*K + L*L)
                            JK=JK+1
                            ARR[JK] = AD
                            BRR[JK] = IN
                                
        JYES = yes_or_no(' Do you want a sorted list ?')
        if JYES == True:
            #print(JK)
            ARR, BRR = PIKSR2(JK, ARR, BRR)
            for JK1 in range(JK, 0, -1):
                print(BRR[JK1], ARR[JK1])
            JK=0
        else:
            JYES = yes_or_no(' Study any particular plane in detail ?')
            if JYES == True:
                H, K, L = list_of_numbers(' Miller indices of plane ?', 3)
                AD, JLOG = DSP(F,H,K,L,1,1)
                print(H,K,L,AD)
                
        JYES = yes_or_no(' Continue analysis?')
        if JYES == False:
            break
        
    return


#
# To obtain the spacing of planes, when the Miller indices are not
# those which are systematically absent. The calculation is
# interactive, and ordered lists of d-spacings available
# MAP: https://www.phase-trans.msm.cam.ac.uk/map/crystal/subs/dspace-b.html
#
def MAP_CRYSTAL_DSPACE(G, JTYPE, II1, JYES):
#
# The metric tensor F converts the components from reciprocal to real space. 
# The metric tensor G converts the components from real to reciprocal
#  space, and G is the inverse of F and vice-versa.
# AAA is the matrix of d-spacings
#
    JK=0
    II1=II1+1

    F = map.MAP_UTIL_INVERS(G)
    
    #
    # Generate a list of spacings with systematically varying Miller indices
    #
    II2=0
    for I1 in range(0, II1):
        H=I1-1
        for I2 in range(0,II1):
            K=I2-1
            for I3 in range(0, II1):
                L=I3-1
                if H == 0 and K == 0 and L == 0:
                    break
                MAP_CRYSTAL_DSP(F, H, K, L, JTYPE, JLOG, AD)
                if JLOG == 1:
                    break
                II2=II2+1
                AAA[II2]=AD

                skip = False
                for JI in range(0, II2-1):
                    #
                    #  Eliminate planes with similar spacings
                    #
                    if AD == AAA[JI]:
                        skip = true
                        break 
                if skip == False:
                    IN=int(H*H+K*K+L*L)
                    JK=JK+1
                    AH[JK]=H
                    AK[JK]=K
                    AL[JK]=L
                    SPCE[JK]=AD
                    MAGN[JK]=IN

    if JYES != 0:
        MAP_UTIL_SORT2(JK, SPCE, MAGN)

    return [JK, AH, AK, AL, SPCE, MAGN]


#
# JR1,JR2 HAVE VALUES OF 0 IF VECTORS ARE REAL, 1 IF RECIPROCAL
#
def ANGLEE(F, G, H, K, L, JR1, HH, KK, LL, JR2):
    
    if JR1 != 0:
        H1, K1, L1 = map.TRANS(F,H,K,L)
        A=map.MAG(H,K,L,H1,K1,L1)
        A=math.sqrt(A)
        AD=1.0/A
    else:
        H1, K1, L1 = map.TRANS(G, H, K, L)
        A = map.MAG(H,K,L,H1,K1,L1)
        A = math.sqrt(A)
        AD = 1.0/A

    if JR2 != 0:
        HH1, KK1, LL1 = map.TRANS(F, HH, KK, LL)
        B=map.MAG(HH,KK,LL,HH1,KK1,LL1)
        B=math.sqrt(B)
        BD=1.0/B
    else:
        HH1, KK1, LL1 = map.TRANS(G, HH, KK, LL)
        B=map.MAG(HH, KK, LL, HH1, KK1, LL1)
        B=math.sqrt(B)
        BD=1.0/B

    if JR1 == 1 and JR2 == 0:
        C = map.MAG(H, K, L, HH, KK, LL)
    if JR1 == 1 and JR2 == 1:
        C=map.MAG(H, K, L, HH1, KK1, LL1)
    if JR1 == 0 and JR2 == 0:
        C=map.MAG(H1, K1, L1, HH, KK, LL)
    if JR1 == 0 and JR2 == 1:
        C=map.MAG(H, K, L, HH, KK, LL)
          
    ANG = C / (A * B)
    if ANG > 1.0:
        ANG = 1.0
    if ANG < -1.0:
        ANG = -1.0
    #print(ANG)
    ANGDEG = math.acos(ANG)
    ANGDEG = ANGDEG * 360.0 / (2 * 3.141592654)
    
    return (ANG, ANGDEG)


#
#
#
def AANG():

    LOGO(5,1.0)
    JYES = 2

    while True:

        if JYES == 2:
            G = np.zeros(9)
            print('''*****************************************************
                     All dimensions in Angstroms, All angles in degrees
                     Input data consist of two vectors which may be
                     real or reciprocal lattice vectors.
                  *****************************************************''')

            J, JTYPE = QUES(1)
            
            G = MET(J)
            
            F = map.MAP_UTIL_INVERS(G)

        H,K,L = list_of_numbers(' Components of Vector 1 ?', 3, float )
        JR1 = REEDI(' Vector 1:  REAL (= 0) or RECIPROCAL (= 1)?  ', 0, 1)
        HH,KK,LL = list_of_numbers(' Components of Vector 2 ?', 3, float)
        JR2 = REEDI(' Vector 2:  REAL  (= 0) or RECIPROCAL (= 1) ?', 0, 1)
        ANG, ANGDEG = ANGLEE(F, G, H, K, L, JR1, HH, KK, LL, JR2)

        if JR1 == 0 and JR2 == 0:
            msg = '''
Real space Vector 1 {0:10.4f} {1:10.4f} {2:10.4f} 
Real space Vector 2 {3:10.4f} {4:10.4f} {5:10.4f} 
Cosine of angle between Vectors 1 & 2 =  {6:10.7f}  
Angle in degrees = {7:10.2f}'''.format(H, K, L, HH, KK, LL, ANG, ANGDEG)

        if JR1 == 1 and JR2 == 1:
            msg = '''
Reciprocal space Vector 1 {0:10.4f} {1:10.4f} {2:10.4f}
Reciprocal space Vector 2 {3:10.4f} {4:10.4f} {5:10.4f}
Cosine of angle between Vectors 1 & 2 = {6:10.7f}
Angle in degrees = {7:10.2f}'''.format(H, K, L, HH, KK, LL, ANG, ANGDEG)

        if JR1 == 0 and JR2 == 1:
           msg = '''
Real space Vector 1 {0:10.4f} {1:10.4f} {2:10.4f}
Reciprocal space Vector 2 {3:10.4f} {4:10.4f} {5:10.4f}
Cosine of angle between Vectors 1 & 2 = {6:10.7f}
Angle in degrees = {7:10.2f}'''.format(H, K, L, HH, KK, LL, ANG, ANGDEG)

        if JR1 == 1 and JR2 == 0:

            msg = '''
Reciprocal space Vector 1 {0:10.4f} {1:10.4f} {2:10.4f}
Real space Vector 2 {3:10.4f} {4:10.4f} {5:10.4f}
Cosine of angle between Vectors 1 & 2 = {6:10.7f}
Angle in degrees = {7:10.2f}'''.format(H, K, L, HH, KK, LL, ANG, ANGDEG)


        JYES = REEDI('Continue Analysis ? 0 = No 1 = Yes same crystal 2 = Yes, new crystal', 0, 2)
        if JYES == 0:
            break

    return


#
#
#
def CONV():
    
    LOGO(6, 1.0)
    JYES = 1

    while True:
        
        if JYES == 1:
            
            G = np.zeros(9)

        print('''*****************************************************
  All dimensions in Angstroms
  All angles in Degrees
  Input data consist of crystal system, lattice type
  and a real or reciprocal lattice vector
*****************************************************''')
        
        J, JTYPE = QUES(1)
        
        G = MET(J)

        F = map.MAP_UTIL_INVERS(G)

        H, K, L = list_of_numbers(' Components of Vector ?', 3, float)
        JR1 = REEDI(' REAL (= 0) or RECIPROCAL (= 1) ?', 0 , 1)
        U, V, W = CONVERT(F, G, H, K, L, JR1)

        if JR1 != 1: 
            msg = ' THE DIRECTION {0:10.6} {1:10.6} {2:10.6} is parallel to THE PLANE NORMAL {3:10.6} {4:10.6} {5:10.6}'.format(H, K, L, U, V, W)
        else:
            msg = ' THE DIRECTION {0:10.6} {1:10.6} {2:10.6} is parallel to THE PLANE NORMAL {3:10.6} {4:10.6} {5:10.6}'.format(U, V, W, H, K, L)
                                
        #print(msg)

        JYES = yes_or_no(' Continue Analysis?')
        if JYES == False:
            break;

    return


#
# Contains rotation martices defining the symmetry operations of a cubic lattice
#
def MAP_CRYSTAL_ORIENT():

    BEC = np.zeros((24, 9))

    BEC[1,1]=1.0
    BEC[1,5]=1.0
    BEC[1,6]=1.0
    BEC[2,2]=1.0
    BEC[2,3]=1.0
    BEC[2,7]=1.0
    BEC[3,0]=-1.0
    BEC[3,5]=-1.0
    BEC[3,7]=-1.0
    BEC[4,1]=-1.0
    BEC[4,3]=-1.0
    BEC[4,8]=-1.0
    BEC[5,2]=-1.0
    BEC[5,4]=-1.0
    BEC[5,6]=-1.0
    BEC[6,0]=1.0
    BEC[6,5]=1.0
    BEC[6,7]=-1.0
    BEC[7,2]=1.0
    BEC[7,4]=1.0
    BEC[7,6]=-1.0
    BEC[8,1]=1.0
    BEC[8,3]=1.0
    BEC[8,8]=-1.0
    BEC[9,2]=-1.0
    BEC[9,3]=-1.0
    BEC[9,7]=1.0
    BEC[10,1]=-1.0
    BEC[10,5]=-1.0
    BEC[10,6]=1.0
    BEC[11,0]=-1.0
    BEC[11,4]=-1.0
    BEC[11,8]=1.0
    BEC[12,0]=1.0
    BEC[12,5]=-1.0
    BEC[12,7]=1.0
    BEC[13,2]=1.0
    BEC[13,4]=-1.0
    BEC[13,6]=1.0
    BEC[14,1]=1.0
    BEC[14,3]=-1.0
    BEC[14,8]=1.0
    BEC[15,1]=-1.0
    BEC[15,5]=1.0
    BEC[15,6]=-1.0
    BEC[16,2]=-1.0
    BEC[16,3]=1.0
    BEC[16,7]=-1.0
    BEC[17,0]=-1.0
    BEC[17,4]=1.0
    BEC[17,8]=-1.0
    BEC[18,2]=-1.0
    BEC[18,4]=1.0
    BEC[18,6]=1.0
    BEC[19,0]=-1.0
    BEC[19,5]=1.0
    BEC[19,7]=1.0
    BEC[20,1]=-1.0
    BEC[20,3]=1.0
    BEC[20,8]=1.0
    BEC[21,1]=1.0
    BEC[21,5]=-1.0
    BEC[21,6]=-1.0
    BEC[22,2]=1.0
    BEC[22,3]=-1.0
    BEC[22,7]=-1.0
    BEC[23,0]=1.0
    BEC[23,4]=-1.0
    BEC[23,8]=-1.0

    return BEC


#
#
#
def PAIR():

    AA = np.zeros(9)
    AINVR = np.zeros(9)
    B = np.zeros(9)
    R = np.zeros(9)
    BEC = np.zeros((24,9))
    
    while True:
        LOGO(8, 1.0)

        print(''' 2B       2A       1B
\\     |     /
 \\    |    /
  \\   |   /
   \\  |  /
    \\ | /
     \\  ------------ 1A
 Arrange your data according to this diagram
 VECTORS 1A & 2A are from CRYSTAL A
 VECTORS 1B & 2B are from CRYSTAL B
 You also need the acute angle between 2A & 2B''')
              
        input(' Press Return when ready ...')

        H, K, L = list_of_numbers(' Components of VECTOR 1A ?', 3)
        H1, K1, L1 = list_of_numbers(' Components of VECTOR 2A ?', 3)
        U1, V1, W1 = list_of_numbers(' Components of VECTOR 1B ?', 3)
        U, V, W = list_of_numbers(' Components of VECTOR 2B ?', 3)
        PHI = REED(' Acute angle (degrees) between VECTORS 2A and 2B ?')
        PHI = PHI * 2.0 * 3.14159 / 360.0
        H, K, L = map.NORM(H, K, L)
        H1, K1, L1 = map.NORM(H1, K1, L1)
        U1, V1, W1 =  map.NORM(U1, V1, W1)
        U, V, W = map.NORM(U, V, W)
        H2 = K*L1-L*K1
        K2 = L*H1-H*L1
        L2 = H*K1-K*H1
        #print(H2, K2, L2)
        H2, K2, L2 = map.NORM(H2, K2, L2)
        Q = H*H1+K*K1+L*L1
        Q = math.acos(Q) if math.isnan(Q) else float("nan")
        M = U*U1+V*V1+W*W1 
        M = math.acos(M) if math.isnan(M) else float("nan")
        P = M-PHI
        A1=H2
        A2=K2
        A3=L2
        A4=H1
        A5=K1
        A6=L1
        A7=H
        A8=K
        A9=L
        C2=math.cos(PHI)
        C3=math.cos(PHI+Q)
        C4=math.cos(P)
        C5=math.cos(Q-P)
        B[0] = A5*A9-A6*A8
        B[1] = A6*A7-A4*A9
        B[2] = A4*A8-A7*A5
        B[3] = A8*A3-A9*A2
        B[4] = A9*A1-A7*A3
        B[5] = A7*A2-A8*A1
        B[6] = A2*A6-A3*A5
        B[7] = A3*A4-A1*A6
        B[8] = A1*A5-A2*A4
        C=A1*B[0]+A2*B[1]+A3*B[2]

        for I in range(0,9):
            AA[I]=B[I]/C

        H3=AA[3]*C2+AA[6]*C3
        K3=AA[4]*C2+AA[7]*C3
        L3=AA[5]*C2+AA[8]*C3
        H4=AA[3]*C4+AA[6]*C5
        K4=AA[4]*C4+AA[7]*C5
        L4=AA[5]*C4+AA[8]*C5
        H3, K3, L3 = map.NORM(H3,K3,L3)
        H4, K4, L4 =  map.NORM(H4,K4,L4)
        X=(K4-V1)*(L3-W)-(L4-W1)*(K3-V)
        Y=(L4-W1)*(H3-U)-(H4-U1)*(L3-W)
        Z=(H4-U1)*(K3-V)-(K4-V1)*(H3-U)
        C10=(H4+U1)*(H3-U)+(K4+V1)*(K3-V)+(L3-W)*(L4+W1)
        X=X/C10
        Y=Y/C10
        Z=Z/C10
        M=math.sqrt(X*X+Y*Y+Z*Z)
        N=2.0*math.atan(M)
        N1=N*360.0/(2.0* math.pi)
        X, Y, Z = map.NORM(X,Y,Z)
        #print(' ******   Rotation Axis = [ {0} {1} {2} ]  ****** Right-handed rotation angle = {3} Degrees'.format(X, Y, Z, N1))
        R = map.MAP_UTIL_ROT(X, Y, Z, N, R)

        XXX=X
        YYY=Y
        ZZZ=Z
        NNNN1=N1
        for I20 in range(0,24):
            for I21 in range(1,9):
                BEC[I20,I21]=0.0

        BEC = MAP_CRYSTAL_ORIENT()

        #print(' Rotation Matrix')
        #print(R[0],R[1],R[2])
        #print(R[3],R[4],R[5])
        #print(R[6],R[7],R[8])
        
        #print(' Inverse of Rotation Matrix')
        AINVR = map.MAP_UTIL_INVERS(R)
        #print(AINVR[0], AINVR[1], AINVR[2])
        #print(AINVR[3], AINVR[4], AINVR[5])
        #print(AINVR[6], AINVR[7], AINVR[8])
        BEC = map.ROTAT(XXX, YYY, ZZZ, NNNN1, BEC)

        JYES = yes_or_no(' Continue Analysis?')
        if JYES == False:
            break
    return

#
# 
#
def DOC():

    print('''** PRACTICAL CRYSTALLOGRAPHY **'/
**               (Documentation)      **
** 1.  Analysis of Electron Diffraction
**
** 2.  Calculation of Metric Tensor  
**
** 3.  Interplanar Spacings
**
** 4.  Angles between Vectors
**
** 5.  Convert between Real & Reciprocal
**
** 6.  Axis-Angle Pairs
**
** 7.  Orientation Relationships
**
** 8.  Four index notation              ',
**''')
    J = REEDI(' Choose document number ... ', 1, 8)
    output = ""
    if J == 1:
        output = '''[1] Analysis of electron diffraction patterns - the camera constant may be known, in which case the program needs an input of two d-spacings from a pair of reciprocal lattice  vectors and the angle between the two vectors concerned.
 The crystal system and lattice type have to be specified.  It is often the case that the camera constant is not  known. The program then operates without d-spacings, needing instead an input of the ratio of the lengths of two reciprocal  lattice vectors, as measured from the diffraction pattern.  If the lattice type is not known then it can be assumed to be  primitive in the first instance. Crystallographically  equivalent solutions can be avoided, but this may increase the computing time.  The program also asks for a "measurement error", which is typically 3%. The program is able to handle any crystal system. The trigonal system is represented by its hexagonal cell.'''

    if J == 2:
        output = '''[2] Computes the Metric Tensor which is in effect a coordinate transformation relating the real and reciprocal axes.'''
          
    if J == 3:
        output='''[3] Calculates the spacing of planes in crystals of arbitrary system, taking account of absences due to lattice types.'''
          
    if J == 4:
        output = '''[4] This program calculates the angles between two vectors in an arbitrary crystal system, whether the vectors are plane normals (i.e., reciprocal vectors) or directions (i.e., real vectors) - any combination of these two kinds of vectors serves as the input.'''
          
    if J == 5:
        ouput = '''[5] Converts the components of a vector from real space to reciprocal space or vice versa, for any arbitrary crystal system. The program can therefore be used to identify a plane normal which is parallel to a specified direction, or vice versa.'''
          
    if J == 6:
        output = '''[6] A powerful program for the calculation of axis-angle pairs relating crystals belonging to the Cubic system. Also calculates rotation matrices and symmetry related descriptions.'''
          
    if J == 7:
        output = '''[7] A very powerful program which calculates the orientation relationship between two crystals, each of which may belong to any arbitrary crystal system, from a knowledge of a pair of vectors (any combination of real or reciprocal vectors will do) from each crystal e.g., (011)A || (111)B and [100]A || [-1 0 1]B. Calculates coordinate transformation matrices, generates lists of planes (and directions) which are parallel in the two crystals, and allows chosen vectors to be examined in detail.'''
          
    if J ==8:
        output = '''[9] The conversion of four index vector notation to three index notation, and vice versa. Miller-Bravais notation is used for plane normals and Weber notation for directions. This is conventional since it allows the correct use of  the Weiss zone law.'''

    #print(textwrap.fill(output))
    #print(textwrap.fill('''Details of the crystallographic methods used can be found in "The Geometry of Crystals", 1987,  by H. K. D. H. Bhadeshia, and in Chapter 1 of "Microstructural Characterisation of High Temperature Materials, both published by the Institute of Metals, 1 Carlton House Terrace, London SW1Y 5DB'''))

    return


#
# Converts from three or four index notation to four or three index notation
# for a hexagonal lattice.
#
def NOTAT():
    
    LOGO(9, 1.0)
    
    while True:
        
        JR1 = REEDI(' Four Index to Three Index Notation (= 0) Three Index to Four Index Notation (= 1) ?', 0, 1)

        if JR1 != 1:
            U, V, T, W = list_of_numbers(' Components of Vector ?', 4)
            JR2 = REEDI(' REAL (= 0) or RECIPROCAL (= 1) ?', 0, 1)
            if JR2 != 1:
                U1=U-T
                U2=V-T
                U3=W
            else:
                U1=U
                U2=V
                U3=W
            #print(' Vector in four index notation  = {0} {1} {2} {3}  Vector in three index notation = {4} {5} {6}'.format(U,V,T,W,U1,U2,U3))
        else:
            U, V, W = list_of_numbers(' Components of Vector ?', 3)
            JR2 = REEDI(' REAL (= 0) or RECIPROCAL (= 1) ?', 0, 1)
            if JR2 == 1:
                U1=U
                U2=V
                U3=W
                U4=-1.0*(U+V)
            else:

                U1=(2.0*U-V)/3.0
                U2=(2.0*V-U)/3.0
                U3=W
                U4=-1.0*(U+V)/3.0
            #print(' Vector in four index notation  = {0} {1} {2} {3}  Vector in three index notation = {4} {5} {6}'.format(U1,U2,U4,U3,U,V,W))

        JYES = yes_or_no(' Continue Analysis ?')

        if JYES == 0:
            break

    return


# To calculate the magnitude of a vector defined in real space
#
# Reference
#  Worked Examples in the Geometry of Crystals
#  by H. K. D. H. Bhadeshia
#  Published by the Institute of Materials, London, 1987.
#  pages 5,6 and equation 4a-c.
# The metric tensor G converts the components from real to reciprocal space.
#
def MAP_CRYSTAL_VECMAG(G,U,V,W,M):

    H, K, L = MAP_UTIL_TRANS(G,U,V,W)
    M=math.sqrt(U*H+V*K+W*L)

    return M


#
#
#
def VECMAG():

    LOGO(10,1.0)
    J, JTYPE = QUES(1,1)
    G = MET(J)

    while True:

        U,V,W = list_of_numbers(' Components of Real Space Vector ?', 3)
        H, K, L = MAP_UTIL_TRANS(G, U, V, W)
        M=math.sqrt(U*H + V*K + W*L)

        #print('    Magnitude of Vector [ {0} {1} {2} ]  = {3} Angstroms'.format(U, V, W, M))

        JYES = yes_or_no(' Continue Analysis ?')
        if JYES == False:
            break

    return


if __name__ == '__main__':
    I = 0
    while I != 11:
        AVER = 1.7
        I = logo3(AVER)
        
        if I == 1:
            electron_diffraction()
        elif I == 2:
            calculate_metric_tensor()
        elif I == 3:
            DSPACE()
        elif I == 4:
            AANG()
        elif I == 5:
            CONV()
        elif I == 6:
            PAIR()
        elif I == 7:
            CORD()
        elif I == 8:
            DOC()
        elif I == 9:
            NOTAT()
        elif I == 10:
            VECMAG()
        elif I == 11:
            exit()

def call_crystal(Crystal_Systemp, Crystal_System_Typep, APp, BPp, CPp, ALPHAp, BETAp, GAMMAp, Camera_Constant_Knownp, Ratio_of_two_g_vectorps, Acute_Angle_Between_Vectorsp, Accuracyp, Max_hkl_indexp, Avoid_Crystallographically_Equivalent_Solutionsp):    # Added this to call crystal from another file
    global Crystal_System, Crystal_System_Type, AP, BP, CP, ALPHA, BETA, GAMMA, Camera_Constant_Known, Ratio_of_two_g_vectors, Acute_Angle_Between_Vectors, Accuracy, Max_hkl_index, Avoid_Crystallographically_Equivalent_Solutions

    Crystal_System=Crystal_Systemp
    Crystal_System_Type = Crystal_System_Typep
    AP = APp
    BP = BPp
    CP = CPp
    ALPHA = ALPHAp
    BETA = BETAp
    GAMMA = GAMMAp
    Camera_Constant_Known=Camera_Constant_Knownp
    Ratio_of_two_g_vectors= Ratio_of_two_g_vectorps
    Acute_Angle_Between_Vectors=Acute_Angle_Between_Vectorsp
    Accuracy= Accuracyp
    Max_hkl_index = Max_hkl_indexp
    Avoid_Crystallographically_Equivalent_Solutions = Avoid_Crystallographically_Equivalent_Solutionsp

    I = 0
    while I != 11:
        AVER = 1.7
        I = logo3(AVER)
        
        if I == 1:
            answers = electron_diffraction()
            return answers
            I=11  # This way it only runs electron diffraction analysis once
        elif I == 2:
            calculate_metric_tensor()
        elif I == 3:
            DSPACE()
        elif I == 4:
            AANG()
        elif I == 5:
            CONV()
        elif I == 6:
            PAIR()
        elif I == 7:
            CORD()
        elif I == 8:
            DOC()
        elif I == 9:
            NOTAT()
        elif I == 10:
            VECMAG()
        elif I == 11:
            exit()
    #return answer
