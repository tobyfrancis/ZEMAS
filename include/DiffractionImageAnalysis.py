import cv2
import numpy as np
import os
import matplotlib.pyplot as plt

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib import patches
import math

import pandas as pd



filepath = r'June10_2021_STO\June10_2021_STO\SAED'
directory = os.listdir(filepath)


def FindCenterOfDiffractionPattern(filepath,filename):
    temp=filepath+'\\' + filename
    im = cv2.imread(temp, cv2.IMREAD_GRAYSCALE)
    imr = im.ravel()
    th = np.mean(imr)
    #print (th , ' is the mean value of the image')
    for ii in range(10):
    # an auotomated threshold
        foreg = imr >= th
        backg = imr < th

        thf = np.mean(imr[foreg])
        thg = np.mean(imr[backg])

        th2 = (thf + thg) / 2.0

        err = np.abs(th - th2)
            
        if err < 0.1:
            break
        else:
           th = th2
           #print (th)


    im_th = im.copy()
    im_th[im >= th2] = 1
    im_th[im < th2] = 0

    # Find the center of mass
    mn = im_th.sum()

    XX, YY = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    cx = np.sum(XX * im_th)/mn
    cy = np.sum(YY * im_th)/mn
 

    x=cx
    y=cy
    

    return x,y

def compare(center, keypoint):   #return a quadrant number: 1, 2, 3, or 4, 0 means unidentified
    quadrant = 0
    if center[0] > keypoint.pt[0]:    # 
        if center[1] < keypoint.pt[1]:   # BOTTOM LEFT
            quadrant = 1
        if center[1] > keypoint.pt[1]: # TOP LEFT
            quadrant = 2
    if center[0] < keypoint.pt[0]:    # 
        if center[1] < keypoint.pt[1]:   # BOTTOM RIGHT
            quadrant = 3
        if center[1] > keypoint.pt[1]:  # TOP RIGHT
            quadrant = 4   
    return quadrant


def distance(kpt1, kpt2):
    #create numpy array with keypoint positions
    arr = np.array([kpt1.pt, kpt2.pt])  # in pixels, not a real unit yet
    #print(arr)
    #scale pixels to array unit
    #arr = arr*unit/pixel_size
    #return distance, calculted by pythagoras
    return np.sqrt(np.sum((arr[0]-arr[1])**2))
    

def distance_from_a_point(point, kpt2):
    #create numpy array with keypoint positions
    arr = np.array([point, kpt2.pt])  # in pixels, not a real unit yet
    #scale pixels to array unit
    #arr = arr*unit/pixel_size
    #return distance, calculted by pythagoras
    return np.sqrt(np.sum((arr[0]-arr[1])**2))

def angle_from_a_keypoint(kpt1, kpt2):
    #create numpy array with keypoint positions
    arr = np.array([kpt1.pt, kpt2.pt])  # in pixels, not a real unit yet
    #scale pixels to array unit
    #arr = arr*unit/pixel_size
    #return distance, calculted by pythagoras
    angle_rad = math.atan2(arr[0][1]-arr[1][1], arr[0][0]-arr[1][0])
    angle_deg = angle_rad *180/np.pi
    return (angle_deg)

def angle_from_a_point(point, kpt2):
    arr = np.array([point, kpt2.pt])  # in pixels, not a real unit yet
    angle_rad = math.atan2(arr[0][1]-arr[1][1], arr[0][0]-arr[1][0])
    angle_deg = angle_rad *180/np.pi
    return (angle_deg)
    
def all_values(obj):
  for attr in dir(obj):
    #print("obj.%s = %r" % (attr, getattr(obj, attr)))
    temp=1

def GetDiffractionSpotProperties(filepath,filename):
    temp=filepath+'\\' + filename
    # Read image
    im = cv2.imread(temp, cv2.IMREAD_GRAYSCALE)

    
    im=cv2.bitwise_not(im)
    params = cv2.SimpleBlobDetector_Params()
    params.filterByArea = True
    #params.minArea = 8
    #print(all_values(params))
    detector = cv2.SimpleBlobDetector_create(params)
    keypoints = detector.detect(im)
    #print(type(keypoints)) # <class 'cv2.KeyPoint'>
    im=cv2.bitwise_not(im)
    # Draw detected blobs as red circles.
    # cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS ensures the size of the circle corresponds to the size of blob
    im_with_keypoints = cv2.drawKeypoints(im, keypoints, np.array([]), (0,0,255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

    # Show keypoints
    #cv2.imshow("Keypoints", im_with_keypoints)
    #cv2.waitKey(0)

    # get distances between blobs
    distances=[]
    j=0
   
    
    
    
    com_x, com_y = FindCenterOfDiffractionPattern(filepath,filename)
    center = (int(com_x), int(com_y))
    #print('Center' , center)
    #im_with_keypoints =cv2.circle(im_with_keypoints, center, 10, (255, 0, 0), 1)  #image, center, radius, color, thickness with -1 for fillin  #FOR COM Center-of-mass
    number_of_keypoints = len(keypoints[0:])
    for ii in range(number_of_keypoints-1):
        #print(number_of_keypoints , ' is the number of key points')
        for i,keypoint in enumerate(keypoints[j+1:]):
            
            d=distance(keypoints[j], keypoint)
            #print(keypoints[0]) # is the same everytime
            #print(keypoint)
            #print("Distance: {0:6.3f} units".format(d))
            distances.append(d)
        j+=1


    # if no beam stopper, the brightest spot is the center
    #find the index of the key point with the largest size
    s=0
    #directbeam_xpos=0
    for i in keypoints:
        temp=i.size
        if s<temp:
            s=temp
            directbeam_xpos=i.pt[0]
            directbeam_ypos=i.pt[1]
            directbeam_keypoint_index=i
    center=(directbeam_xpos,directbeam_ypos)

    distances_from_center=[]
    spot_intensity = []
    min_d = 1000000
    
    angles_from_center=[]
    
    min_angle = 1000000

    for i,keypoint in enumerate(keypoints[0:]):
            
            d=distance_from_a_point(center, keypoint)
            #print("Distance: {0:6.3f} units".format(d))
            distances_from_center.append(d)
            spot_intensity.append(keypoint.size)
            if np.absolute(d) < np.absolute(min_d):
                min_d=d
                min_d_keypoint=i

            theta=angle_from_a_point(center, keypoint)

            if theta <0:
                theta = 360+theta

            #print( i, theta)
            angles_from_center.append(theta)
            
            if np.absolute(theta) < np.absolute(min_angle):
                min_angle=theta
                min_angle_keypoint=i


    # quadrants

    quad1 = [] # indices_of_keypoints_in_quad  # NO REASON TO MAKE FOUR LISTS - CAN CHANGE THIS
    quad2 = [] # indices_of_keypoints_in_quad
    quad3 = [] # indices_of_keypoints_in_quad
    quad4 = [] # indices_of_keypoints_in_quad
    for i,keypoint in enumerate(keypoints[0:]):
        if i != directbeam_keypoint_index:
            
            quadrant=compare(center, keypoint)
            if quadrant == 1:
                quad1.append(i)
            if quadrant == 2:
                quad2.append(i)
            if quadrant == 3:
                quad3.append(i)
            if quadrant == 4:
                quad4.append(i)
        
    #get the size of the spots in quad k
    

    intensities=[]
    for i in quad1:
        intensities.append(keypoints[i].size)
        #im_with_keypoints =cv2.circle(im_with_keypoints, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 10, (0, 255, 255), 1)
    quad1intensity=sum(intensities)

    intensities=[]
    for i in quad2:        
        intensities.append(keypoints[i].size)
        #im_with_keypoints =cv2.circle(im_with_keypoints, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 10, (0, 255, 0), 1)
    quad2intensity=sum(intensities)

    intensities=[]
    for i in quad3:        
        intensities.append(keypoints[i].size)
        #im_with_keypoints =cv2.circle(im_with_keypoints, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 10, (255, 0, 0), 1)
    quad3intensity=sum(intensities)

    intensities=[]
    for i in quad4:
        #print('indices' , i)
        #print(keypoints[i].size, 'INTENSITY OF DIFFRACTED BEAM')
        #im_with_keypoints =cv2.circle(im_with_keypoints, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 10, (255, 255, 0), 1)
        intensities.append(keypoints[i].size)
    quad4intensity=sum(intensities)

    #where are the alpha and beta axis, that depends on the microscope reading , also need some calibration from intensity to degrees
    tilt_top_right = - (quad2intensity - quad3intensity)
    tilt_top_left =  - (quad1intensity - quad4intensity)


    

    angles_from_a_key_point = []
    j=0
    for ii in range(number_of_keypoints-1):
        for i, keypoint in enumerate(keypoints[0:]):
            theta = angle_from_a_keypoint(keypoints[j], keypoint)  # the first argument's index shouldn't be 1 necessarily, it should be min_d_keypoint ideally to get a spot from the first ring
            angles_from_a_key_point.append(theta)
        j+=1



    #print('min_d_keypoint', min_d_keypoint)
    x= keypoints[min_d_keypoint].pt[0]
    y=keypoints[min_d_keypoint].pt[1]
    radius_spot=(keypoints[min_d_keypoint].size)/2
    im_with_keypoints = cv2.circle(im_with_keypoints, (int(x), int(y)), int(radius_spot), (255, 0, 0), 1)  # green circle is the diffraction spot closest to center of mass

    #print('Number of distances from center', len(distances_from_center))   
    count_distances=len(distances)
    #print (count_distances , ' is the number of measured distances')

    if len(distances) > 0:
    #print(type(distances))
        plt.subplot(2, 3, 1)
        plt.imshow(im_with_keypoints)
        plt.xlabel(filename.replace('.png',''))
        plt.subplot(2, 3, 2)
        #plt.scatter(range(count_distances), distances)
        #plt.gca().set_ylim(bottom=0)
        plt.scatter(y=angles_from_center, x=distances_from_center)
        plt.xlabel('Distances from center')
        plt.ylabel('Angles from center')
        plt.gca().set_ylim(bottom=0)  #-180 
        plt.gca().set_ylim(top=360)  # 180 
        plt.subplot(2, 3, 3)
        plt.hist(distances, density=False, bins=max(2,int(max(distances))))  # density=False would make counts
        plt.ylabel('Count')
        plt.xlabel('Distance')
        plt.subplot(2, 3, 3) # if you change to plt.subplot(2, 2, 3), the subplots will overlay
        plt.hist(distances_from_center, density=False, bins=max(2,int(max(distances))))  # density=False would make counts
        plt.ylabel('Count')
        plt.xlabel('Distance betw. points / from center')
        plt.subplot(2, 3, 4) # if you change to plt.subplot(2, 2, 3), the subplots will overlay
        plt.hist(angles_from_a_key_point, density=False, bins=max(2,int(max(angles_from_a_key_point))))  # density=False would make counts
        plt.ylabel('Count')
        plt.xlabel('Angles from one point to all others (degrees)')
        plt.xlabel('Angle (deg)')
        plt.subplot(2, 3, 5)
        plt.scatter(y=spot_intensity,x=distances_from_center)
        plt.xlabel('distance from center' )
        plt.ylabel('diffracted beam intensity')

        
        #plt.scatter(range(len(angles_from_a_key_point)), angles_from_a_key_point)
        #plt.gca().set_ylim(bottom=0)
        #plt.ylabel('Angle')
        #plt.xlabel('Index (no meaning)')
        plt.subplot(2, 3, 6)
        #plt.bar([1,2], [tilt_top_left,tilt_top_right])
        #plt.xlabel('Alpha is towards top left, Beta is towards top right')
        plt.pie([quad1intensity, quad2intensity,quad3intensity,quad4intensity], explode=None, labels=['bottom left', 'top left',' bottom right', 'top right'])
        #plt.show()
        plt.subplot(1, 1, 1)
        plt.imshow(im_with_keypoints)
        plt.xlabel(filename.replace('.png',''))
        #plt.show()

        df=pd.DataFrame()
        df['Angle']=angles_from_center
        df['Distance']=distances_from_center
        df=df[df['Distance'].gt(0)]  #  Ineffecient
        #print(df)

        return df, im_with_keypoints

def GetAngleANdDHKLDistances(df):
        
        if df is not None:
            #print(df.shape[0])
            temp=df.shape[0]
            
            df2=pd.DataFrame()
            k=0
            Spot1 =[]
            Spot2 =[]
            A=[]
            B=[]
            C=[]
            for i in range(temp):
                
                for j in range(i+1,temp):
                    k=+1
                    angle_betw_two_vectors = df['Angle'].values[j]-df['Angle'].values[i]
                    angle_betw_two_vectors=abs(angle_betw_two_vectors)
                    if angle_betw_two_vectors < 90:
                        Spot1.append(i)
                        Spot2.append(j)
                        A.append(angle_betw_two_vectors)
                        B.append(df['Distance'].values[i])
                        C.append(df['Distance'].values[j])

            df2['Spot 1'] = Spot1
            df2['Spot 2'] = Spot2
            df2['Angle'] = A
            df2['Vector 1 1/dhkl'] = B
            df2['Vector 2 1/dhkl'] = C
            #print(df2)
            return df2
        else:
            print('No acute angles')



directory=directory[25:27]
#print(directory)
for filename in directory:
    
    if ' ' in filename:
        continue
    if '.png' in filename:
        
        df=GetDiffractionSpotProperties(filepath,filename)
        df2 = GetAngleANdDHKLDistances(df)
        df2=df2[df2['Angle']>10]
        



        