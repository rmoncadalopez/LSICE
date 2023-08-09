import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
#from grain import Grain
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')

#Use code for SimInput
from SimInput import genParam

import math
import datetime
import csv
import cv2
import numpy as np
import glob
import shapely.geometry
import random

#GSD Error Plotting and Outputting Img versus Sim

def GSD_Comparison_Error(FileImg, FileData, DayNumber):

    #Import Image GSD File Data
    dataImg = np.loadtxt(FileImg, delimiter = " ")

    #Import data as vectors
    img_gsdD = []
    img_gsdP = []
    img_gsdC = []
    img_gsdConc = []

    for i in range(len(dataImg)):
        img_gsdD.append(dataImg[i][0])
        img_gsdP.append(dataImg[i][1])
        img_gsdC.append(dataImg[i][2])
        img_gsdConc.append(dataImg[i][3])

    #Import simulation GSD Data
    dataSim = np.loadtxt(FileData, delimiter = " ")

    #Import data as vectors
    gsdD = []
    gsdP = []
    gsdC = []
    gsdConc = []

    for i in range(len(dataSim)):
        gsdD.append(dataSim[i][0])
        gsdP.append(dataSim[i][1])
        gsdC.append(dataSim[i][2])
        gsdConc.append(dataSim[i][3])

    #WARNING!!!
    #Only for initial conditions the fine concentration be made equal in the Simfile
    gsdConc[-1] = img_gsdConc[-1]
    
    #Export concentration
    Fine_Conc = img_gsdConc[-1]
    Coarse_Conc = 0
    for i in range(len(gsdD)-1):
        Coarse_Conc += gsdConc[i]
    
    conc_text = [(Fine_Conc, Coarse_Conc)]
    np.savetxt(outDir+'Initial_Concentration.dat', conc_text, fmt='%.5f %.5f')
    
    #Resave sim file to update this minor change
    vec_export_text = []
    for i in range(len(gsdD)):
        vec_export_text.append([gsdD[i], gsdP[i], gsdC[i], gsdConc[i]])

    print(vec_export_text)
    np.savetxt(outDir+'GSD_SimInput.dat', vec_export_text, fmt='%.5f %.5f %d %.5f')

    #Compare with SIM Init. (% Error wrt to img)

    error_gsdP = np.zeros(len(gsdD))
    error_gsdConc = np.zeros(len(gsdD))
    Perror_gsdP = np.zeros(len(gsdD))
    Perror_gsdConc = np.zeros(len(gsdD))
    #Count does not have to be the same but Percent pass and Concentration have to be similar

    #Positive is pending grains, negative is exceeding
    for i in range(len(gsdD)):
        print("img_gsdP: ", img_gsdP[i])
        print("gsdP: ", gsdP[i])
        if int(img_gsdP[i]) == 0 and int(gsdP[i]) > 0:
            Perror_gsdP[i] = (img_gsdP[i] - gsdP[i])*100/gsdP[i]
        elif int(img_gsdP[i]) == 0 and int(gsdP[i]) == 0:
            Perror_gsdP[i] = 0
        else:
            Perror_gsdP[i] = (img_gsdP[i] - gsdP[i])*100/img_gsdP[i]
            
        if int(img_gsdConc[i]) == 0 and int(gsdConc[i]) > 0:
            Perror_gsdConc[i] = -100
        elif int(img_gsdConc[i]) == 0 and int(gsdConc[i]) == 0:
            Perror_gsdConc[i] = 0
        else:
            Perror_gsdConc[i] = (img_gsdConc[i] - gsdConc[i])*100/img_gsdConc[i]
           
        error_gsdP[i] = (img_gsdP[i] - gsdP[i])
        error_gsdConc[i] = (img_gsdConc[i] - gsdConc[i])

    Mean_error_gsdP = np.mean(error_gsdP)
    Mean_error_gsdConc = np.mean(error_gsdConc)
    Mean_Perror_gsdP = np.mean(Perror_gsdP)
    Mean_Perror_gsdConc = np.mean(Perror_gsdConc)

    #Export Error as text for Reference
    vec_export_text = [[Mean_error_gsdP, Mean_error_gsdConc, Mean_Perror_gsdP, Mean_Perror_gsdConc]]
    for i in range(len(error_gsdP)):
        vec_export_text.append([error_gsdP[i], error_gsdConc[i], Perror_gsdP[i], Perror_gsdConc[i]])

    np.savetxt(outDir+'GSD_Error_Day_'+ str(DayNumber)+'.dat', vec_export_text, fmt='%.5f %.5f %5f %.5f')

    #Plot and export Error graph
    fig = plt.subplots(figsize=(10,10))
    plt.plot(gsdD, Perror_gsdP, 'r.', label = 'Percent Pass % Error')
    plt.plot(gsdD, Perror_gsdConc, 'b.', label = 'Concentration % Error')
    plt.axhline(y=0.0, color='k', linestyle='--')
    plt.title('Initial Grain Size Distribution % Error: Day ' + str(DayNumber))
    plt.savefig(outDir + 'Initial_GSD_Per_Error_Day_'+ str(DayNumber)+'.png', format='png')
    plt.close()

    fig = plt.subplots(figsize=(10,10))
    plt.plot(gsdD, error_gsdP, 'r.', label = 'Percent Pass Error')
    plt.axhline(y=0.0, color='k', linestyle='--')
    plt.title('Initial Grain Size Distribution Error (% Pass): Day ' + str(DayNumber))
    plt.savefig(outDir + 'Initial_GSD_Error_Pass_Day_'+ str(DayNumber)+'.png', format='png')
    plt.close()

    fig = plt.subplots(figsize=(10,10))
    plt.plot(gsdD, error_gsdConc, 'b.', label = 'Concentration Error')
    plt.axhline(y=0.0, color='k', linestyle='--')
    plt.title('Initial Grain Size Distribution Error (Concentration): Day ' + str(DayNumber))
    plt.savefig(outDir + 'Initial_GSD_Error_Conc_Day_'+ str(DayNumber)+'.png', format='png')
    plt.ylim(-5, 5)
    plt.close()
    
    #Plot FSD
    fig, ax1 = plt.subplots(figsize=(10,10))
       #plt.figure(figsize=(10,10))
    color = 'tab:red'
    ###plt.xscale('log') ###Do for GSD only, not histogram
    plt.xlim(0, MaxD*1.1)
    ax1.set_xlabel('Grain Size')
    ax1.set_ylabel('% Pass', color=color)
    ax1.plot(gsdD, gsdP, 'ro-')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.bar(gsdD, gsdConc, align='edge', color = 'red', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.25) ) # A bar chart
    plt.ylim(0,100)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('# Floes', color=color)  # we already handled the x-label with ax1
    ax2.plot(gsdD, gsdC, 'b*-')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.bar(gsdD, gsdC, align='edge', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.5) ) # A bar chart
    plt.ylim(0,350)

    #fig.tight_layout()  # otherwise the right y-label is slightly clipped

    #Revert X-Axis
    #axG = plt.gca()
    #axG.set_xlim(axG.get_xlim()[::-1])
    plt.title('Initial Grain Size Distribution and Count')
    plt.savefig(outDir + 'Initial_GSD.png', format='png')
    plt.close()
    
    return Mean_error_gsdP, Mean_error_gsdConc, Mean_Perror_gsdP, Mean_Perror_gsdConc


#Points area function for Initial GSD
def PointsArea(pointList):
    n = len(pointList)
    Area = 0
    for i in range(n-1):
        Area += ( float(pointList[i][0]) * float(pointList[i+1][1]) -  float(pointList[i][1]) * float(pointList[i+1][0]) )
    Area += ( float(pointList[n-1][0]) * float(pointList[0][1]) -  float(pointList[n-1][1]) * float(pointList[0][0]) )
    return 0.5*abs(Area)


#Mean Caliper Diameter for Initial GSD
def MeanCaliperD(Area):
    return 2 * np.sqrt(Area/math.pi)

#Intersection Function
def intersection_grains(p1,nGrains,pointsDat):
    #Check intersection with existing grains
    for j in range(nGrains):
        #p1 = shapely.geometry.Polygon(pointsDat[i])
        p2 = shapely.geometry.Polygon(pointsDat[j])
        if (p1.intersects(p2) == True):
            return True
#    #With Walls #NO WALLS HERE
#    for (j) in range(nWalls):
#        #p1 = shapely.geometry.Polygon(pointsDat[i])
#        p2 = shapely.geometry.Polygon(pointsWall[j])
#        if (p1.intersects(p2) == True):
#            return True
                
#    #With Additional
#    p2 = shapely.geometry.Polygon(polynew1)
#    if (p1.intersects(p2) == True):
#        return True
#
#    #With Additional 2
#    p2 = shapely.geometry.Polygon(polynew2)
#    if (p1.intersects(p2) == True):
#        return True
        
    return False
    

def AppendVec (VecText):
    Vecfull = []
    for i in range(len(VecText)):
        Vecfull.append(VecText[i])
    return Vecfull

#Code to Populate an area with random grains
#Start using some real info and fill in the blanks with additional grains
# Output contact and updated position file
outDir = "./Input/SeaIceWGJ2018/"

figDir             = "./Visualization/Initial/"
#grainDir            = "./Input/grains/"
grainDir2            = "./Input/grainsIceWGJ2018/"

inDir1 = "./Input/grainsIceWGJ2018/"
File1 = "initialpositions.dat"
posFile = np.loadtxt(inDir1 + File1, delimiter = " ")

#stats_file = np.loadtxt("grainspropsimage.dat", delimiter = " ")
#wallDir             = "./Input/grainsWallsF/"
#wallDir             = "./Input/grains 0WallsTS/"
#VelFile           =  "./Input/SeaIceOcean/initialvelocities.dat"
#ThickFile           =  "./Input/SeaIceOcean/initialthick.dat"
#TempFile           =  "./Input/SeaIceOcean/initialtemper.dat"
#posFile            =  "./Input/SeaIceOcean/initialpositions.dat"
#wallposFile            =  "./Input/SeaIce/initialpositionsWallFrev.dat"
#morphFile          =  "./Input/SeaIceOcean/morphologies.dat" #Pool of New Morphs
#morphFile2          =  "./Input/SeaIce/morphologiesAug.dat" #With Pool of New Morphs
#wallmorphFile          =  "./Input/SeaIce/morphologiesWallF.dat"
#morphFile          =  "./Input/mainCompression3/morphologies.dat"
#pointsFile         = "./Output/SeaIce/pointsout0.dat"
#posFile0           = "../input//positionsEquil.dat"

#VelV               =  np.loadtxt(VelFile) #For writing
#ThickV               =  np.loadtxt(ThickFile) #For writing
#TempV                =  np.loadtxt(TempFile) #For writing
#morphsVOr           =  np.loadtxt(morphFile2) #For choosing
#morphsV            =  np.loadtxt(morphFile) #For writing
#nGrains            = np.loadtxt(morphFile, dtype = int)[0]
##nWalls            = np.loadtxt(wallmorphFile, dtype = int)[0]
#posRot             = np.loadtxt(posFile) #For writing
##posWall            = np.loadtxt(wallposFile)


#morphsVOr = [103, 122, 8, 0, 183, 108, 136, 114, 44, 147, 192, 102, 124, 152, 145, 15, 127, 158, 55] #Diverse sizes
#morphsVOr = [102, 192, 124, 152, 15, 127, 145, 1, 158, 55, 186, 200, 25, 187, 164, 163] #Medium size
#morphsVOr = [114, 44, 147, 192, 102, 124, 152, 145, 15, 127, 158, 55]
#morphsVOr = [103, 8]

#morphsVTX = AppendVec(morphsV)
#ThickVTX = AppendVec(ThickV)
#TempVTX = AppendVec(TempV)
#VelVTX = AppendVec(VelV)
#posRotVTX = AppendVec(posRot)

morphsVTX = [0]
Thickness = random.random()+1 #For all floes
ThickVTX = []
Temperature = -1.8 #For all floes
TempVTX = []
VelVTX = []
posRotVTX = []


npoints = 100
#npointsWall = 700

#Number of grains to input
nNewGrains = 493#154 #34
morphsVOr = range(0,nNewGrains) #Choose from all these options

#nGrainsMorph = nGrains
nGrainsMorph = 0

multSc_im = (0.08425538203)**0.5  #Ref size for 2018 is 1943 km2 versus 23060.84137573685 u2 so the scale is 0.08425538203
#extent_conc = 259*339.66  #2020
extent_conc = 0.076011266724044 * 1474000 #For 2018 Example, see imageprocess2_label_v2018.py (1340 x, 11000 y in pixels)
x1 = 0
x2 = 1400*multSc_im #360 #Based on image values
y1 = 0
y2 = 1400*multSc_im #300 #Based on image values

#polynew1 = [[-10,5100],[-10,2400],[100,2500],[500,2500],[600,2650],[1000,3100],[1500,3500],[1800,4000],[2000,4300],[2100,4900],[1900, 5150]]
#
#polynew2 = [[6000,4500],[6000,4800],[4000,3800],[3000,3500],[2500,3000],[2000,2800],[2000,2800],[2000,2000],[2100,1500],[2000,1000],[6000, 1000]]

#Import all points
pointsDat = []
#for i in range(nGrains):
#    positionP = posRot[i]
#    print(positionP)
#    ng = np.loadtxt(morphFile, dtype = int)[i+1]
#    pointsFile = grainDir + "grainproperty" + str(ng) + ".dat"
#    pointsExtract = open(pointsFile,"r")
#    pointList = pointsExtract.readlines()[4]
#    pointList = pointList.split()
#    pointVec = []
#    jj = 0
#    for j in range(npoints):
#        floatP = [float(pointList[jj]), float(pointList[jj+1])]
#        Xg = floatP[0] + positionP[0]
#        Yg = floatP[1] + positionP[1]
#        vec = [Xg,Yg]
#        pointVec.append(vec)
#        jj = jj + 2
#    pointsDat.append(pointVec)

#pointsWall = []
#for i in range (nWalls):
#    positionP = posWall[i]
#    print(positionP)
#    ng = np.loadtxt(wallmorphFile, dtype = int)[i+1]
#    pointsFile = wallDir + "grainproperty" + str(ng) + ".dat"
#    pointsExtract = open(pointsFile,"r")
#    pointList = pointsExtract.readlines()[4]
#    pointList = pointList.split()
#    pointVec = []
#    jj = 0
#    for j in range(npointsWall):
#        floatP = [float(pointList[jj]), float(pointList[jj+1])]
#        Xg = floatP[0] + positionP[0]
#        Yg = floatP[1] + positionP[1]
#        vec = [Xg,Yg]
#        pointVec.append(vec)
#        jj = jj + 2
#    pointsWall.append(pointVec)
    
tol = 20

#Pending TODO LIST !!!!!
#0 Delete empty grains by using mass as filter (get mass) and diagnose those of Xpos = 0 or relocate (done)
#2 Insert extra grains to Match Sat and Sim GSD (Pending for Coarse and Fine) and print Initial GSD and error wrt to Sat. Data (done)
#3 Match units with WG picture or scale (done for Initial Geometry)  (Pending for Simulation critical mass and time stepping parameters and constraints, space seems to be fine)
#1 Transform pixel data of WG to Meter units and export GSD that way (done)
#1.5 Re-scale grains and re-run this (check area) (done)
#5 Fit parameters to evolve Big floe GSD similar
#4 Define single GSD time evolving indicators to compare satellite data img points and sim
#6 Consider finer ice
#7 Run side by side


#Loop through number of new grains
#invertG = range(nNewGrains-1,-1,-1) #Image processing require inverting Y coord and grain order
invertG = range(nNewGrains)
#invertG1 = np.array(range(153,-1,-1))
#invertG2 = np.array(range(187,153,-1))
#invertG = []
##print(invertG1)
##print(invertG2)
#for i in range(len(invertG1)+len(invertG2)):
#    if i < 154:
#        invertG.append(invertG1[i])
#    else:
#        invertG.append(invertG2[i-len(invertG1)])
#print(invertG1)
#print(invertG2)
#print(invertG)
#exit(1)

TagG = []
GoodG = []
GoodGPoints = []

#diamlist = []
#for i in range(nNewGrains):
#    pointsFileTemp = grainDir2 + "grainproperty" + str(i) + ".dat"
#    pointsExtract = open(pointsFileTemp, "r")
#    pointList = pointsExtract.readlines()[4]
#    pointList = pointList.split()
#    areaAreaArea = np.reshape(pointList, (int(0.5*len(pointList)),2))
#    RefArea = PointsArea(areaAreaArea)
#    DD = MeanCaliperD(RefArea)
#    diamlist.append([i,DD])
#    print([i,DD])
#exit(1)

val_replace = 0 #77
for i in range(nNewGrains):
    print("Step: ", i)
    positionPFull = posFile[i]
#
#    #Exchange 1 with Zero
#    if invertG[i] == 1 or invertG[i] == 8 or invertG[i] == 13 or invertG[i] == 50 or invertG[i] == 28 or invertG[i] == 182 or invertG[i] == 172 or invertG[i] == 170 or invertG[i] == 158 or invertG[i] == 122 or invertG[i] == 73 or invertG[i] == 38 or invertG[i] == 25: # or invertG[i] == 135 or invertG[i] == 40 or invertG[i] == 8 or invertG[i] == 150:
#             #Choose adequate sample grain to replace
#        pointsFileTemp = grainDir2 + "grainproperty" + str(val_replace) + ".dat"
#        GoodG.append(val_replace)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 113 or invertG[i] == 177:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(151) + ".dat"
#        GoodG.append(151)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 40:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(31) + ".dat"
#        GoodG.append(31)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 44:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(55) + ".dat"
#        GoodG.append(55)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 45:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(39) + ".dat"
#        GoodG.append(39)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 19:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(23) + ".dat"
#        GoodG.append(23)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 40:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(36) + ".dat"
#        GoodG.append(36)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 51:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(42) + ".dat"
#        GoodG.append(42)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 75:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(74) + ".dat"
#        GoodG.append(74)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 20:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(16) + ".dat"
#        GoodG.append(16)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 24:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(81) + ".dat"
#        GoodG.append(81)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 132:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(129) + ".dat"
#        GoodG.append(129)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 117:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(23) + ".dat"
#        GoodG.append(23)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 100:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(56) + ".dat"
#        GoodG.append(56)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 87:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(12) + ".dat"
#        GoodG.append(12)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 121:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(130) + ".dat"
#        GoodG.append(130)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 140:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(111) + ".dat"
#        GoodG.append(111)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 144:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(37) + ".dat"
#        GoodG.append(37)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 146:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(69) + ".dat"
#        GoodG.append(69)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 72:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(60) + ".dat"
#        GoodG.append(60)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 22:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(37) + ".dat"
#        GoodG.append(37)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 157:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(109) + ".dat"
#        GoodG.append(109)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 64:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(65) + ".dat"
#        GoodG.append(65)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 171:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(26) + ".dat"
#        GoodG.append(26)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 107:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(147) + ".dat"
#        GoodG.append(147)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 109:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(125) + ".dat"
#        GoodG.append(125)
#    elif invertG[i] == 49:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(84) + ".dat"
#        GoodG.append(84)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 74:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(75) + ".dat"
#        GoodG.append(75)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 86 or invertG[i] == 84:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(83) + ".dat"
#        GoodG.append(83)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 41 or invertG[i] == 11 or invertG[i] == 4:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(0) + ".dat"
#        GoodG.append(0)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 13 or invertG[i] == 155 or invertG[i] == 167 or invertG[i] == 181 or invertG[i] == 187 or invertG[i] == 48 or invertG[i] == 186 or invertG[i] == 77 or invertG[i] == 183 or invertG[i] == 142 or invertG[i] == 92 or invertG[i] == 127 or invertG[i] == 74:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(55) + ".dat"
#        GoodG.append(55)
#        #GoodG.append(invertG[i])
#    elif invertG[i] == 166 or invertG[i] == 169:
#        pointsFileTemp = grainDir2 + "grainproperty" + str(183) + ".dat"
#        GoodG.append(183)
#        #GoodG.append(invertG[i])
   # else:
    pointsFileTemp = grainDir2 + "grainproperty" + str(invertG[i]) + ".dat"
    GoodG.append(invertG[i])
    
#    pointsFileTemp = grainDir2 + "grainproperty" + str(invertG[i]) + ".dat"
#    GoodG.append(invertG[i])
    TagG.append(invertG[i])
    pointsExtract = open(pointsFileTemp, "r")
    pointList0 = pointsExtract.readlines()
    #print(pointList0)
    pointList = pointList0[4]
    pointList = pointList.split()
    pointsExtract.close()
    #print(pointList)
    #print(np.reshape(pointList, (int(0.5*len(pointList)),2)))
    GoodGPoints.append(np.reshape(pointList, (int(0.5*len(pointList)),2)))

#    #Ref grain to scale with approx. 1039 km2 versus u2, then scale to km both grains and positions, then adjust only grains to match sat grain size range
#    if invertG[i] == 37:
#            areaAreaArea = np.reshape(pointList, (int(0.5*len(pointList)),2))
#            RefArea = PointsArea(areaAreaArea)
#            print("Ref Area: ", RefArea)
#
#    #Correct too small grains with a sample grain
#    if ((pointList[0] == pointList[2]) and (pointList[1] == pointList[3]) and (pointList[0]== pointList[4]) and (pointList[1] == pointList[5]) and (pointList[0] == pointList[6]) and (pointList[1] == pointList[7]) and (pointList[0]== pointList[8]) and (pointList[1] == pointList[9])):
#        print("Point Replace")
#
#        pointsFileTemp = grainDir2 + "grainproperty" + str(val_replace) + ".dat"
#        print(GoodG)
#        print(i)
#        GoodG[-1] = invertG[i]
#        #GoodG[-1] = val_replace  #Replace latest element
#        pointsExtract = open(pointsFileTemp, "r")
#        pointList = pointsExtract.readlines()[4]
#        pointList = pointList.split()
#        GoodGPoints[-1] = np.reshape(pointList, (int(0.5*len(pointList)),2))
    
    pointTest = []
    jj = 0
    Rot = positionPFull[2]
    for j in range(npoints):
        floatP0 = [float(pointList[jj]), float(pointList[jj+1])]
        #Angle rotation
        ppoint1 = floatP0[0]*math.cos(Rot) + floatP0[1]*math.sin(Rot)
        ppoint2 = -floatP0[0]*math.sin(Rot) + floatP0[1]*math.cos(Rot)
        floatP = [ppoint1, ppoint2]
        #Displace to final position for checking intersection
        Xg = floatP[0] + positionPFull[0]
        Yg = floatP[1] + positionPFull[1]
        vec = [Xg,Yg]
        pointTest.append(vec)
        jj = jj + 2

    #append Grain to PointsDat list
    pointsDat.append(pointTest)

    #Modify Morphology Vector
    morphsVTX.append(GoodG[-1])
      
    #Modify Velocity Vector
    VelVTX.append([0.0000,0.0000,0.0000])

    #Modify Thickness
    min_adm = 0.05 #Avoid zero thickness
    mu_th, sigma_th = 0.5, 0.1 # mean and standard deviation
    srand = max(np.random.normal(mu_th, sigma_th), min_adm)
    ThickVTX.append(srand) ##Around 0.5 m thickness
    #ThickVTX.append(random.random()+1) ##Around 1m thickness

    #Modify Temperature
    TempVTX.append(Temperature)
    
    #Increase positions
    posRotVTX.append(positionPFull)
    
    nGrainsMorph = nGrainsMorph + 1
        
print(len(posRotVTX))
print(len(VelVTX))
print(len(morphsVTX))
print("GoodG= ", GoodG)

#New grains to insert (grain Pool is zeroX)
zeroX = []
print("ZeroX= ", zeroX)
for i in range(len(zeroX)):
    inter = True #Assume Not Good
    #Choose Random Morphology, rechoose if it doesnt work
    #rand_N = int(np.random.choice(morphsVOr[1:]))
    #rand_N = int(np.random.choice(morphsVOr))
    #print("Try: ", rand_N)
    #niterbas = 0
    #nmax = 30 #Final number for giving up
    
    #Assume we can get it well given low populated space (WARNING!!!)
    while (inter == True):
        niterlimit = 0
        #maxiter = 30 #How many time to try the same morphology in different positions
        
    #While Loop
    #while (inter == True and niterlimit < maxiter):
        #print("Attemp #: ", niterlimit)
        
        #Choose Random Position
        tol = int(min(x2,y2)*0.07) #20 To avoid going out
        
        if i == 0:
            #positionP1 = 34
            #positionP2 = 40
            #rand_N = 103
            positionP1 = np.random.choice(range(x1+tol, x2-tol))
            positionP2 = np.random.choice(range(y1+tol, y2-tol))
        elif i == 1:
            #positionP1 = 46
            #positionP2 = 40
            #rand_N = 8
            positionP1 = np.random.choice(range(x1+tol, x2-tol))
            positionP2 = np.random.choice(range(y1+tol, y2-tol))
        else:
            positionP1 = np.random.choice(range(x1+tol, x2-tol))
            positionP2 = np.random.choice(range(y1+tol, y2-tol))
        
        positionP = [positionP1, positionP2]
        
        #Random Rotation(later)
        #Rot = 3.141592653589793*0.5
        #Rot = random.randrange(0, 2*math.pi, 0.01)
        #Rot = random.uniform(0, 2*3.141592653589793)
        Rot = 0.000
        
        positionPFull = [positionP1, positionP2, Rot]
        
        #Make Set of Points
        pointsFileTemp = grainDir2 + "grainproperty" + str(zeroX[i]) + ".dat"
        pointsExtract = open(pointsFileTemp, "r")
        pointList = pointsExtract.readlines()[4]
        pointList = pointList.split()
        pointTest = []
        jj = 0
        for j in range(npoints):
            floatP0 = [float(pointList[jj]), float(pointList[jj+1])]
            #Angle rotation
            ppoint1 = floatP0[0]*math.cos(Rot) + floatP0[1]*math.sin(Rot)
            ppoint2 = -floatP0[0]*math.sin(Rot) + floatP0[1]*math.cos(Rot)
            floatP = [ppoint1, ppoint2]
            #Displace to final position for checking intersection
            Xg = floatP[0] + positionP[0]
            Yg = floatP[1] + positionP[1]
            vec = [Xg,Yg]
            pointTest.append(vec)
            jj = jj + 2

        #Check if there is intersection for all point Sets
        p1 = shapely.geometry.Polygon(pointTest)
        inter = intersection_grains(p1, len(pointsDat), pointsDat)
         
        if inter == False:
            #print("Insert more grains: ", rand_N)
            
            #append Grain to PointsDat list
            pointsDat.append(pointTest)

            #Modify Morphology Vector
            morphsVTX.append(zeroX[i])
          
            #Modify Velocity Vector
            VelVTX.append([0.0000,0.0000,0.0000])

            #Modify Thickness
            min_adm = 0.05 #Avoid zero thickness
            mu_th, sigma_th = 0.5, 0.1 # mean and standard deviation
            srand = max(np.random.normal(mu_th, sigma_th), min_adm)
            ThickVTX.append(srand) ##Around 0.5 m thickness
            #ThickVTX.append(random.random()+1) ##Around 1 m thickness

            #Modify Temperature
            TempVTX.append(Temperature)
            
            #Increase positions
            posRotVTX.append(positionPFull)
            
            nGrainsMorph = nGrainsMorph + 1

        else:
            niterlimit = niterlimit + 1
            print("Tried ", niterlimit, " tries.")
    
#    niterbas = niterbas + 1
#    if (niterbas == nmax):
#        rand_N = 64 #Very small grain that will work
            

#update Number morphs
morphsVTX[0] = nGrainsMorph
print("Ngrains: ", nGrainsMorph)

#Write new input files from final result
np.savetxt(outDir+'initialpositions.dat', posRotVTX, fmt='%.5f %.5f %.5f')
np.savetxt(outDir+'initialvelocities.dat', VelVTX, fmt='%.5f %.5f %.5f')
np.savetxt(outDir+'initialtemper.dat', TempVTX, fmt='%.5f')
np.savetxt(outDir+'initialthick.dat', ThickVTX, fmt='%.5f')
np.savetxt(outDir+'morphologies.dat', morphsVTX, fmt='%d')
    
    
#FOR PLOTTING
# Collect in patches
patches = []
#print(pointsDat)
#print("npts0: ", len(pointsDat))
#pointIndiv = np.split(pointsDat, nGrains)
#print("npts: ", len(pointIndiv))
#print("npts2: ", len(pointIndiv[0]))
for n in range(len(pointsDat)):
    #print(grains[n]._points)
    poly = Polygon(pointsDat[n], True)
    patches.append(poly)
    
#patcheswall = []
#for i in range(nWalls):
#    poly = Polygon(pointsWall[i], True)
#    patcheswall.append(poly)
    
#patcheswall.append(Polygon(polynew1,True))
#patcheswall.append(Polygon(polynew2,True))

#PLOT OF INITIAL STATE
fig, ax = plt.subplots(figsize = (15,15))
#ax.set_xlim(-lX,lX)
#ax.set_ylim(-lY,lY)
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
ax.autoscale_view()

pCol = PatchCollection(patches, facecolors='white',edgecolors='black', lw=0.1)
ax.set_facecolor('xkcd:blue')
ax.add_collection(pCol)

#xstats = np.zeros(len(stats_file))
#ystats = np.zeros(len(stats_file))
#tagstats = np.zeros(len(stats_file))
#
#for i range (len(stats_file)):
#    xstats[i] = stats_file[i][1]
#    ystats[i] = stats_file[i][2]
#    tagstats[i] = stats_file[i][0]

for i in range(len(GoodG)):
    plt.text(posRotVTX[i][0], posRotVTX[i][1], str(TagG[i]), size = 8)

#for i in range(len(stats_file)):
#    plt.text(xstats[i], ystats[i], str(tagstats[i]), size = 10)

#pCol2 = PatchCollection(patcheswall, facecolors='black',edgecolors='black', lw=0.1)
#ax.add_collection(pCol2)
       
#plt.axis('off')
plt.savefig(figDir + 'InitialTryPWGJ2018.png', format='png')
plt.close(fig)

#Also generated case scenarios
#SimInputV = genParam()
#SimInputV[i][2] *= 0.06 #ADJUST MEANWHILE
#np.savetxt(outDir+'SimInput.dat', SimInputV, fmt='%.5f %.5f %.5f %.5f %.5f')
#print("--SCENARIOS GENERATED--")
#for i in range(len(SimInputV)):
#    print(SimInputV[i][0], SimInputV[i][1], SimInputV[i][2], SimInputV[i][3], SimInputV[i][4])



#Obtain FSD from Entered Grains (for consistency same functions as in LS-DEM Sea Ice Frac)
MaxD = 0.0
MinD = 4000000000000000000.000
totArea = 0
#Get Max and Min Data
for i in range(len(GoodG)):
    areaList = GoodGPoints[i]
    Area = PointsArea(areaList)
    Diameter = MeanCaliperD(Area)
    print(Area)
    totArea += Area
    if Diameter > MaxD:
        MaxD = Diameter
        max_index = i
    if Diameter < MinD:
        MinD = Diameter
    
#Fix for convenience
#Data
MaxD = 42.84503#42.64503
MinD = 1.70271#1.90271

#IMG
#MaxD = 48.66452
#MinD = 7.46926
    
DistRange = MaxD - MinD
bin_size = 10
gsdD = np.zeros(bin_size+1)
gsdConc = np.zeros(bin_size+1)
gsdC = np.zeros(bin_size+1)
gsdP = np.zeros(bin_size+1)

#Define Diameters (NOTE: Descending order)
for i in range(bin_size+1):
    gsdD[i] = MaxD - (float(i)/float(bin_size))*DistRange

#Define Percent Passes
for i in range(bin_size+1):
    AreaSum = 0.0 #Temporal Area for Adding at %Passes
    for j in range(len(GoodG)):
        areaList = GoodGPoints[j]
        Area = PointsArea(areaList)
        Diameter = MeanCaliperD(Area)
        if ( Diameter <= gsdD[i]):
            AreaSum += Area
    gsdP[i] = AreaSum*100.0/totArea

#Define Count and Concentration
Domain_Area = extent_conc
print(gsdD)
for i in range(bin_size):
    count_size = 0
    AreaSumConc = 0.0 #Temporal Area for Adding at %Concentration in a bin
    for j in range(len(GoodG)):
        areaList = GoodGPoints[j]
        Area = PointsArea(areaList)
        Diameter = MeanCaliperD(Area)

        if ( Diameter <= gsdD[i] and  Diameter >= gsdD[i+1] ):
            count_size += 1
            AreaSumConc += Area
    gsdConc[i] = (AreaSumConc * 100) / Domain_Area
    gsdC[i] = count_size;

gsdConc[bin_size] = 0
gsdC[bin_size] = 0


#Export FSD as text for Reference
vec_export_text = []
for i in range(len(gsdD)):
    vec_export_text.append([gsdD[i], gsdP[i], gsdC[i], gsdConc[i]])

print(vec_export_text)
np.savetxt(outDir+'GSD_SimInput.dat', vec_export_text, fmt='%.5f %.5f %d %.5f')

##Plot FSD
#fig, ax1 = plt.subplots(figsize=(10,10))
#   #plt.figure(figsize=(10,10))
#color = 'tab:red'
####plt.xscale('log') ###Do for GSD only, not histogram
#plt.xlim(0, MaxD*1.1)
#ax1.set_xlabel('Grain Size')
#ax1.set_ylabel('% Pass', color=color)
#ax1.plot(gsdD, gsdP, 'ro-')
#ax1.tick_params(axis='y', labelcolor=color)
#ax1.bar(gsdD, gsdConc, align='edge', color = 'red', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.25) ) # A bar chart
#plt.ylim(0,100)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#color = 'tab:blue'
#ax2.set_ylabel('# Floes', color=color)  # we already handled the x-label with ax1
#ax2.plot(gsdD, gsdC, 'b*-')
#ax2.tick_params(axis='y', labelcolor=color)
#ax2.bar(gsdD, gsdC, align='edge', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.5) ) # A bar chart
#plt.ylim(0,120)
#
##fig.tight_layout()  # otherwise the right y-label is slightly clipped
#
##Revert X-Axis
##axG = plt.gca()
##axG.set_xlim(axG.get_xlim()[::-1])
#plt.title('Initial Grain Size Distribution and Count')
#plt.savefig(outDir + 'Initial_GSD.png', format='png')
#plt.close()


#Find Error for Initial Conditions  #UNCOMMENT!!!

#Import Image Based GSD
FileImg = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/Papers/Sea Ice Gifs/Multi_spec/2018ExResults/0_GSD_Data.dat"
FileData = outDir + "GSD_SimInput.dat"
DayNumber = 0

#Function to compare image and simulation at specific steps based on temporal fitting, here it's easy because it's step 0
[Mean_error_gsdP, Mean_error_gsdConc, Mean_Perror_gsdP, Mean_Perror_gsdConc] = GSD_Comparison_Error(FileImg, FileData, DayNumber)
print(Mean_error_gsdP, Mean_error_gsdConc, Mean_Perror_gsdP, Mean_Perror_gsdConc)

##Check intersection with each other
#for i in range(nGrains):
#    #print(i)
#    for (j) in range(i+1,nGrains,1):
#        #print(j)
#        p1 = shapely.geometry.Polygon(pointsDat[i])
#        p2 = shapely.geometry.Polygon(pointsDat[j])
#        if (p1.intersects(p2) == True):
#            print("We got intersection between: ", i , " and: ", j)
#
##With Walls
#for i in range(nGrains):
#    #print(i)
#    for (j) in range(nWalls):
#        #print(j)
#        p1 = shapely.geometry.Polygon(pointsDat[i])
#        p2 = shapely.geometry.Polygon(pointsWall[j])
#        if (p1.intersects(p2) == True):
#            print("We got WALL intersection between: ", i , " and: ", j)
#
##With Additional
#for i in range(nGrains):
#        p1 = shapely.geometry.Polygon(pointsDat[i])
#        p2 = shapely.geometry.Polygon(polynew1)
#        if (p1.intersects(p2) == True):
#            print("We got POLYNEW 1 intersection BY: ", i)
#
##With Additional 2
#for i in range(nGrains):
#        p1 = shapely.geometry.Polygon(pointsDat[i])
#        p2 = shapely.geometry.Polygon(polynew2)
#        if (p1.intersects(p2) == True):
#            print("We got POLYNEW 2 intersection BY: ", i)
#
