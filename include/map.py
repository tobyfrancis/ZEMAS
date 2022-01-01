import re
import cmath
#import keyboard
import math
import numpy as np
import sys
import textwrap

#
#
#
def TEST(ACCU, F, GV1, GV2, RR, AA, H, K, L, HH, KK, LL, JTYPE, J40, J50, J60, AEQ):
    Answer = []
    AZ = np.zeros(9)
    
    if JTYPE != 1:
        JLOG = TYPE(H, K, L, JTYPE)
        JLOG2 = TYPE(HH, KK, LL, JTYPE)
        if JLOG == 1 or JLOG2 == 1:
            return [J60, AEQ, Answer]

    H1, K1, L1 = TRANS(F, H, K, L)
    A = MAG(H, K, L, H1, K1, L1)
    A = math.sqrt(A)
    AD = 1.0 / A
    HH1, KK1, LL1 = TRANS(F, HH, KK, LL)
    B = MAG(HH, KK, LL, HH1, KK1, LL1)
    B = math.sqrt(B)
    BD = 1.0/B
    if J40 == True:
        
        DUM1 = abs( (BD/GV1) - 1.0)
        DUM2 = abs( (AD/GV1) - 1.0)
        DUM3 = abs( (BD/GV2) - 1.0)
        DUM4 = abs( (AD/GV2) - 1.0)

        if DUM3 >= ACCU or DUM4 >= ACCU:
            return [J60, AEQ, Answer]

        if DUM1 < ACCU or DUM2 < ACCU:
            return [J60, AEQ, Answer]

    R=A/B
    if R < 1.0:
        R = B/A
        
    if abs( (RR/R) - 1.0 ) >= ACCU:
        return [J60, AEQ, Answer]
    
    C = MAG(H, K, L, HH1, KK1, LL1)
    ANGLE=C/(A*B)
    AA2=math.cos(AA)
    if abs(abs(ANGLE) - AA2) >= ACCU:
        return [J60, AEQ, Answer]

    U1 = K * LL - L * KK
    V1 = L * HH - H * LL
    W1 = H * KK - K * HH
    
    U1, V1, W1 = NORM(U1, V1, W1)
    if J50 == True:

        AZ[0]=abs(U1)
        AZ[1]=abs(V1)
        AZ[2]=abs(W1)
        AZ = PIKSRT(AZ)
        for J70 in range(0, J60):
            if J60 == 1:
                break
            if AEQ[J70, 0] == AZ[0] and AEQ[J70, 0] == AZ[0] and AEQ[J70, 1] == AZ[1]:
                return [J60, AEQ, Answer]

        AEQ[J60, 0] = AZ[0]
        AEQ[J60, 1] = AZ[1]
        AEQ[J60, 2] = AZ[2]

    if J60 != 5:
        J60 = J60 + 1

    AAA = math.acos(ANGLE) * 360.0 / (2.0 * math.pi)
    #print('{0:2.0f}. {1:2.0f}. {2:2.0f}. {3:1.4f}     {4:2.0f}. {5:2.0f}. {6:2.0f}. {7:1.4f}     {8: } {9: } {10: } {11:2.1f}'.format(H, K, L, AD, HH, KK, LL, BD, U1, V1, W1, AAA))
    Answer = [H, K, L, AD, HH, KK, LL, BD, U1, V1, W1, AAA]    
    return [J60, AEQ, Answer]

#
# Sorts an array arr into ascending numerical order, by straight insertion.
# n is  input;
# arr is replaced on output by its sorted rearrangement
#
def PIKSRT(ARR):

    for i in range(0, len(ARR)):
        save = ARR[i]
        j = i
        while j > 0 and ARR[j - 1] > save:
            ARR[j] = ARR[j - 1]
            j -= 1
        ARR[j] = save    
    return ARR


#
# 
#
def MAG(H, K, L, H1, K1, L1):

    val = H*H1+K*K1+L*L1
    return val


#
#
#
def TYPE(H, K, L, JTYPE):

    calc = True
    cond = False
    cond_val = False
    JLOG = -1
    if JTYPE == 1:
        calc = False
        JLOG = 0
    elif JTYPE == 2:
        A = H + K + L
    elif JTYPE == 3:
        cond = True
        J1 = ODD(H)
        J2 = ODD(K)
        J3 = ODD(L)
        cond_val = (J1 == 0  and  J2 == 0  and  J3 == 0) or (J1 == 1  and  J2 == 1  and  J3 == 1)        
    elif JTYPE == 4:
        A=K+L
    elif JTYPE == 5:
        A=H+L
    if JTYPE == 6:
        A = H + K
        
    if cond == True:
        if cond_val:
            JLOG = 0
        else:
            JLOG = 1
    else:
        JODD = ODD(A)
        if JODD == 1:
            JLOG = 1
        else:
            JLOG = 0
        
    return JLOG


#
# 
#
def ODD(A):
    JODD = 1
    if int(A/2.0) == (A/2.0):
        JODD=0
    return JODD


#
# 
#

def analyze(J21, J22, J20, ACCU, F, GV1, GV2, RR, AA, JTYPE, J40, J50, J60, AEQ):

    #print('------------------------------------------------------------------')
    #print('    Vector 1       Vector 2          Zone Axis     ')
    #print('  (h  k  l)      d   (h  k  l)     d    [U      V      W]    Angle')
    answers = []
    for I1 in range(1, J22+1):
        H = I1 - 1
        for I2 in range(1, J20+1):
            K = I2 + J21
            for I3 in range(1, J20+1):
                L = I3 + J21
                for J1 in range(1, J22+1):
                    HH = J1 - 1
                    for J2 in range(1, J20+1):
                        KK = J2 + J21
                        for J3 in range(1, J20+1):
                            LL = J3 + J21
                            first = (H == 0 and K == 0 and L == 0)
                            second = (HH == 0 and KK == 0 and LL == 0)
                            third = (HH == H and KK == K and LL == L)
                            fourth = (-1*HH == H and -1*KK == K and -1*LL == L)
                            cond =  first or second or third or fourth
                            if cond == True:
                                #print( first, second, third, fourth, 'so TEST method is not run and we get no result')
                                continue
                            else:
                                J60, AEQ, Answer = TEST(ACCU, F, GV1, GV2, RR, AA, H, K, L, HH, KK, LL, JTYPE, J40, J50, J60, AEQ)
    
                                if len(Answer)>0:
                                    answers.append(Answer)

    return answers


#
#
#
def TRANS(LR, X1, X2, X3):

    X4 = LR[0] * X1 + LR[1] * X2 + LR[2] * X3
    X5 = LR[3] * X1 + LR[4] * X2 + LR[5] * X3
    X6 = LR[6] * X1 + LR[7] * X2 + LR[8] * X3
    return [X4, X5, X6]


#
#
#
def NORM(X, Y, Z):

    C = math.sqrt(X * X + Y * Y + Z * Z)
    if C == 0:
        X = float("nan")
        Y = float("nan")
        Z = float("nan")
    else:
        X = X/C
        Y = Y/C
        Z = Z/C
    
    return [X, Y, Z]


#
# Invert a 3x3 matrix  is the inverse of B
#
def MAP_UTIL_INVERS(B):

    A = np.zeros(9)
    
    C1=(B[4]*B[8]-B[5]*B[7])*B[0]+(B[5]*B[6]-B[3]*B[8])*B[1]+(B[3]*B[7]-B[4]*B[6])*B[2]
    A[0]=(B[4]*B[8]-B[5]*B[7])/C1
    A[3]=(B[5]*B[6]-B[3]*B[8])/C1
    A[6]=(B[3]*B[7]-B[4]*B[6])/C1
    A[1]=(B[7]*B[2]-B[8]*B[1])/C1
    A[4]=(B[8]*B[0]-B[6]*B[2])/C1
    A[7]=(B[6]*B[1]-B[7]*B[0])/C1
    A[2]=(B[1]*B[5]-B[2]*B[4])/C1
    A[5]=(B[2]*B[3]-B[0]*B[5])/C1
    A[8]=(B[0]*B[4]-B[1]*B[3])/C1

    return A


#
#
#
def ROTAT(P1, P2, P3, THETA, BEC):

    SYM = np.zeros(9)
    ANS = np.zeros(9)
    R = np.zeros(9)
    
    P1, P2, P3 =  NORM(P1,P2,P3)
    THETA=THETA*0.017453292
    A5=1.0
    R = MAP_UTIL_ROT(P1,P2,P3,THETA,R)
    print(' The 23 Equivalent Axis Angle Pairs  No.  Axis Angle')
    for I1 in range(1, 24):
        for I2 in range(0, 9):
            SYM[I2] = BEC[I1, I2]
        
        ANS = MAP_UTIL_PROD(SYM, R)
        ANS[0], ANS[1], ANS[2] = NORM(ANS[0], ANS[1], ANS[2])
        ANS[3], ANS[4], ANS[5] = NORM(ANS[3], ANS[4], ANS[5])
        ANS[6], ANS[7], ANS[8] = NORM(ANS[6], ANS[7], ANS[8])
        THETA=ANS[0] + ANS[4] + ANS[8]
        if THETA >= -0.999999:#)  GOTO 999
            P1=ANS[7]-ANS[5]
            P2=ANS[2]-ANS[6]
            P3=ANS[3]-ANS[1]
            return BEC
        
# ******* NEXT FEW LINES FOR SYMMETRIC ROT MAT

        THETA=-1.0
        P1=math.sqrt((ANS[0] + 1.0)/2.0)
        P2=math.sqrt((ANS[4] + 1.0)/2.0)
        P3=math.sqrt((ANS[8] + 1.0)/2.0)

        P1, P2, P3 = NORM(P1,P2,P3)
        THETA=(THETA-1.0)/2.0
        THETA=math.acos(THETA)*360.0/(2.0*math.pi)

        print(I1,P1,P2,P3,THETA)

    return BEC

  
#
#
#
def MAP_UTIL_ROT(P1, P2, P3, THETA, R):

    A5=1.0
    R[0]=P1*P1*(A5-math.cos(THETA))+math.cos(THETA)
    R[1]=P1*P2*(A5-math.cos(THETA))-P3*math.sin(THETA)
    R[2]=P1*P3*(A5-math.cos(THETA))+P2*math.sin(THETA)
    R[3]=P2*P1*(A5-math.cos(THETA))+P3*math.sin(THETA)
    R[4]=P2*P2*(A5-math.cos(THETA))+math.cos(THETA)
    R[5]=P2*P3*(A5-math.cos(THETA))-P1*math.sin(THETA)
    R[6]=P3*P1*(A5-math.cos(THETA))-P2*math.sin(THETA)
    R[7]=P3*P2*(A5-math.cos(THETA))+P1*math.sin(THETA)
    R[8]=P3*P3*(A5-math.cos(THETA))+math.cos(THETA)

    return R

#
# Evaluate the product A of two 3x3 matrices B and C
#
def MAP_UTIL_PROD(B,C):

    A = np.zeros(9)
    
    A[0] = B[0]*C[0]+B[1]*C[3]+B[2]*C[6]
    A[1] = B[0]*C[1]+B[1]*C[4]+B[2]*C[7]
    A[2] = B[0]*C[2]+B[1]*C[5]+B[2]*C[8]
    A[3] = B[3]*C[0]+B[4]*C[3]+B[5]*C[6]
    A[4] = B[3]*C[1]+B[4]*C[4]+B[5]*C[7]
    A[5] = B[3]*C[2]+B[4]*C[5]+B[5]*C[8]
    A[6] = B[6]*C[0]+B[7]*C[3]+B[8]*C[6]
    A[7] = B[6]*C[1]+B[7]*C[4]+B[8]*C[7]
    A[8] = B[6]*C[2]+B[7]*C[5]+B[8]*C[8]

    return A
