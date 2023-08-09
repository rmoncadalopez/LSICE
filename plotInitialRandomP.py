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

import math
import datetime
import csv
import cv2
import numpy as np
import glob
import shapely.geometry

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
outDir = "./Input/SeaIceRandomP/"

figDir             = "./Visualization/Initial/"
grainDir            = "./Input/grains/"
grainDir2            = "./Input/grainsAug/"
#wallDir             = "./Input/grainsWallsF/"
#wallDir             = "./Input/grains 0WallsTS/"
VelFile           =  "./Input/SeaIce/initialvelocities.dat"
ThickFile           =  "./Input/SeaIce/initialthick.dat"
TempFile           =  "./Input/SeaIce/initialtemper.dat"
posFile            =  "./Input/SeaIce/initialpositions.dat"
#wallposFile            =  "./Input/SeaIce/initialpositionsWallFrev.dat"
morphFile          =  "./Input/SeaIce/morphologies.dat" #Pool of New Morphs
morphFile2          =  "./Input/SeaIce/morphologiesAug.dat" #With Pool of New Morphs
#wallmorphFile          =  "./Input/SeaIce/morphologiesWallF.dat"
#morphFile          =  "./Input/mainCompression3/morphologies.dat"
#pointsFile         = "./Output/SeaIce/pointsout0.dat"
#posFile0           = "../input//positionsEquil.dat"

VelV               =  np.loadtxt(VelFile) #For writing
ThickV               =  np.loadtxt(ThickFile) #For writing
TempV                =  np.loadtxt(TempFile) #For writing
morphsVOr           =  np.loadtxt(morphFile2) #For choosing
morphsV            =  np.loadtxt(morphFile) #For writing
nGrains            = np.loadtxt(morphFile, dtype = int)[0]
#nWalls            = np.loadtxt(wallmorphFile, dtype = int)[0]
posRot             = np.loadtxt(posFile) #For writing
#posWall            = np.loadtxt(wallposFile)

morphsVTX = AppendVec(morphsV)
ThickVTX = AppendVec(ThickV)
TempVTX = AppendVec(TempV)
VelVTX = AppendVec(VelV)
posRotVTX = AppendVec(posRot)

npoints = 100
#npointsWall = 700

nNewGrains = 500

nGrainsMorph = nGrains

x1 = 0
x2 = 6000
y1 = 0
y2 = 6000

#polynew1 = [[-10,5100],[-10,2400],[100,2500],[500,2500],[600,2650],[1000,3100],[1500,3500],[1800,4000],[2000,4300],[2100,4900],[1900, 5150]]
#
#polynew2 = [[6000,4500],[6000,4800],[4000,3800],[3000,3500],[2500,3000],[2000,2800],[2000,2800],[2000,2000],[2100,1500],[2000,1000],[6000, 1000]]

#Import all points
pointsDat = []
for i in range(nGrains):
    positionP = posRot[i]
    print(positionP)
    ng = np.loadtxt(morphFile, dtype = int)[i+1]
    pointsFile = grainDir + "grainproperty" + str(ng) + ".dat"
    pointsExtract = open(pointsFile,"r")
    pointList = pointsExtract.readlines()[4]
    pointList = pointList.split()
    pointVec = []
    jj = 0
    for j in range(npoints):
        floatP = [float(pointList[jj]), float(pointList[jj+1])]
        Xg = floatP[0] + positionP[0]
        Yg = floatP[1] + positionP[1]
        vec = [Xg,Yg]
        pointVec.append(vec)
        jj = jj + 2
    pointsDat.append(pointVec)

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
    
#Loop through number of new grains
for i in range(nNewGrains):
    print("Step: ", i)
    inter = True #Assume Not Good
    #Choose Random Morphology, rechoose if it doesnt work
    rand_N = int(np.random.choice(morphsVOr[1:]))
    print("Try: ", rand_N)
    niterbas = 0
    nmax = 20
    while (inter == True and niterbas < nmax):
        niterlimit = 0
        maxiter = 20 #How many time to try the same morphology in different positions
        
        #While Loop
        while (inter == True and niterlimit < maxiter):
            print(niterlimit)
            #Choose Random Position
            tol = 300 # To avoid going out
            positionP1 = np.random.choice(range(x1+tol, x2-tol))
            positionP2 = np.random.choice(range(y1+tol, y2-tol))
            positionP = [positionP1, positionP2]
            
            #Random Rotation(later)
            Rot = 0.0000
            positionPFull = [positionP1, positionP2, Rot]
            
            #Make Set of Points
            pointsFileTemp = grainDir2 + "grainproperty" + str(rand_N) + ".dat"
            pointsExtract = open(pointsFileTemp, "r")
            pointList = pointsExtract.readlines()[4]
            pointList = pointList.split()
            pointTest = []
            jj = 0
            for j in range(npoints):
                floatP = [float(pointList[jj]), float(pointList[jj+1])]
                Xg = floatP[0] + positionP[0]
                Yg = floatP[1] + positionP[1]
                vec = [Xg,Yg]
                pointTest.append(vec)
                jj = jj + 2

            #Check if there is intersection for all point Sets
            p1 = shapely.geometry.Polygon(pointTest)
            inter = intersection_grains(p1, len(pointsDat), pointsDat)
             
            if inter == False:
                print("Insert more grains: ", rand_N)
                
                #append Grain to PointsDat list
                pointsDat.append(pointTest)

                #Modify Morphology Vector
                morphsVTX.append(rand_N)
              
                #Modify Velocity Vector
                VelVTX.append([0.0000,0.0000,0.0000])

                #Modify Thickness
                ThickVTX.append(ThickV[0])

                #Modify Temperature
                TempVTX.append(TempV[0])
                
                #Increase positions
                posRotVTX.append(positionPFull)
                
                nGrainsMorph = nGrainsMorph + 1


            else:
                niterlimit = niterlimit + 1
        
        niterbas = niterbas + 1
        if (niterbas == nmax):
            rand_N = 64 #Very small grain that will work
            

#update Number morphs
morphsVTX[0] = nGrainsMorph

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
fig, ax = plt.subplots(figsize = (5,5))
#ax.set_xlim(-lX,lX)
#ax.set_ylim(-lY,lY)
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
ax.autoscale_view()

pCol = PatchCollection(patches, facecolors='white',edgecolors='black', lw=0.1)
ax.set_facecolor('xkcd:blue')
ax.add_collection(pCol)

#pCol2 = PatchCollection(patcheswall, facecolors='black',edgecolors='black', lw=0.1)
#ax.add_collection(pCol2)
       
#plt.axis('off')
plt.savefig(figDir + 'InitialTryP.png', format='png')
plt.close(fig)




















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
