import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
#from grain import Grain
from multiprocessing import Pool
import matplotlib
#matplotlib.use('Agg')
from scipy.optimize import curve_fit
import sys
import shutil

from PIL import Image

#import shapely.geometry

import math
import datetime
import csv
import cv2
import numpy as np
import glob


#Accumulated Sum
def acc_sum(vec):
    out_vec = np.zeros(len(vec))
    for i in range(len(vec)):
        if i < len(vec) - 1:
            out_vec[i] = np.sum(vec[0:i+1])
        else:
            out_vec[i] = np.sum(vec[0:len(vec)])
    return out_vec
    
#Zeros filter for better Nfloes Count and FSD
def zeros_filter(gsd_dataD0, gsd_dataNF0):
    gsd_dataD = gsd_dataD0
    gsd_dataNF = gsd_dataNF0
    ctnz = 0
    #nonz_st = 0
    #nonz_def = False
    newD = [] #Eliminating all zeros
    newNF = []  #Eliminating all zeros
    #Count how many nonzero values in array
    for i in range(len(gsd_dataNF0)):
        if gsd_dataNF0[i]>0.0:
            ctnz += 1
            #if nonz_def == False:
            #    nonz_st = i
            #    nonz_def = True
            newD.append(gsd_dataD0[i])  #Eliminating all zeros
            newNF.append(gsd_dataNF0[i])  #Eliminating all zeros
    #Eliminate zeros and use new data if at least 2 points, else leave zeros as it is
    if ctnz > 1:
        #gsd_dataD = gsd_dataD0[nonz_st:]
        #gsd_dataNF = gsd_dataNF0[nonz_st:]
        gsd_dataD = newD
        gsd_dataNF = newNF

    return gsd_dataD, gsd_dataNF

#Zeros filter for better SLope Calc. (remove flat ends)
def zeros_filter_slope(gsd_dataD0, gsd_dataNF0):
    gsd_dataD = gsd_dataD0
    gsd_dataNF = gsd_dataNF0
    ctnz = 0
    #nonz_st = 0
    #nonz_def = False
    newD = [] #Eliminating all zeros
    newNF = []  #Eliminating all zeros
    #Count how many nonzero values in array
    for i in range(len(gsd_dataNF0)):
        if i < len(gsd_dataNF0) - 1 and i > 0:
            if gsd_dataNF0[i]!=gsd_dataNF0[i+1] or gsd_dataNF0[i]!=gsd_dataNF0[i-1]:
                ctnz += 1
                #if nonz_def == False:
                #    nonz_st = i
                #    nonz_def = True
                newD.append(gsd_dataD0[i])  #Eliminating all zeros
                newNF.append(gsd_dataNF0[i])  #Eliminating all zeros
        elif i == 0:
            if gsd_dataNF0[0]!=gsd_dataNF0[1]:
                ctnz += 1
                #if nonz_def == False:
                #    nonz_st = i
                #    nonz_def = True
                newD.append(gsd_dataD0[i])  #Eliminating all zeros
                newNF.append(gsd_dataNF0[i])  #Eliminating all zeros
        else:
            if gsd_dataNF0[-2]!=gsd_dataNF0[-1]:
                ctnz += 1
                #if nonz_def == False:
                #    nonz_st = i
                #    nonz_def = True
                newD.append(gsd_dataD0[i])  #Eliminating all zeros
                newNF.append(gsd_dataNF0[i])  #Eliminating all zeros
    #Eliminate zeros and use new data if at least 2 points, else leave zeros as it is
    if ctnz > 1:
        #gsd_dataD = gsd_dataD0[nonz_st:]
        #gsd_dataNF = gsd_dataNF0[nonz_st:]
        gsd_dataD = newD
        gsd_dataNF = newNF

    return gsd_dataD, gsd_dataNF



# objective function
def objective(x, a, b):
    return a * x + b


#GSD Slope Linear
def gsd_slope_linear(objective, x, y):
    # curve fit for all results
    poptMean, _ = curve_fit(objective, x, y)
    # summarize the parameter values
    aMean = poptMean[0]
    bMean = poptMean[1]
    return aMean

#GSD Slope Exponential
def gsd_slope_exp(x, y):
    # curve fit for all results
    fit = np.polyfit(x, np.log(y), 1)
    return np.exp(fit[0])



#GSD Error Plotting and Outputting Img versus Sim
def GSD_Comparison_Error(FileImg, FileData, DayNumber, outDir):

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

    print(FileData)
    #Import simulation GSD Data
    dataSim = np.loadtxt(FileData, delimiter = " ")
    print(dataSim)

    #Import data as vectors
    gsdD = []
    gsdP = []
    gsdC = []
    gsdConc = []

    for i in range(len(dataSim)):
        if math.isnan(dataSim[i][0]):
           gsdD.append(0)
        else:
            gsdD.append(dataSim[i][0])
        if math.isnan(dataSim[i][1]):
           gsdP.append(0)
        else:
            gsdP.append(dataSim[i][1])
        if math.isnan(dataSim[i][2]):
           gsdC.append(0)
        else:
            gsdC.append(dataSim[i][2])
        if math.isnan(dataSim[i][3]):
           gsdConc.append(0)
        else:
            gsdConc.append(dataSim[i][3])

    #Compare with SIM Init. (% Error wrt to img)

    error_gsdP = np.zeros(len(gsdD))
    error_gsdConc = np.zeros(len(gsdD))
    Perror_gsdP = np.zeros(len(gsdD))
    Perror_gsdConc = np.zeros(len(gsdD))
    #Count does not have to be the same but Percent pass and Concentration have to be similar

    #Positive is pending grains, negative is exceeding
    for i in range(len(gsdD)):
        if img_gsdP[i] == 0 and gsdP[i] != 0:
            Perror_gsdP[i] = (img_gsdP[i] - gsdP[i])*100/gsdP[i]
        elif img_gsdP[i] == 0 and gsdP[i] == 0:
            Perror_gsdP[i] = 0
        else:
            Perror_gsdP[i] = (img_gsdP[i] - gsdP[i])*100/img_gsdP[i]
            
        if img_gsdConc[i] == 0:
            Perror_gsdConc[i] = -100
        elif img_gsdConc[i] == 0 and gsdConc[i] == 0:
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
    plt.close()
    
    return Mean_error_gsdP, Mean_error_gsdConc, Mean_Perror_gsdP, Mean_Perror_gsdConc



#def plotSim(figDir,morphFile,posFile,morphDir,vidName,makePics,video):

#Folders for output flexibility

#Main directory for visualization

#Choose year for Analysis
#year = 2020;
#year = 2018;
#year = int(sys.argv[2])

create_folder = True
year = 2018

if create_folder:
    folder_title = sys.argv[1]
    mainVfolder = "./ProcessGSD_" + str(folder_title) + "/"
    mainOutputfolder = "./GSD_textfiles/" #2018    
    #Check if folder already exists
    if os.path.exists(mainVfolder) == False:
        #Create folders
        os.mkdir(mainVfolder)
else:
    mainVfolder = "./ProcessGSD_" + str(folder_title) + "/" 
    mainOutputfolder = "./GSD_textfiles/" #2018  


# if year == 2018:
#     #Field Data from Input (2018)
#     #inputFolder0 = "./Input/SeaIceWGJ2018_MoreUniform1/"
#     inputFolder0 = "./Input/SeaIceWGJ2018/"
#     inputFolder = "./Input/SeaIceWGJ2018/"
#     inputFolder2 = "./Input/SeaIceWGJ2018/2018Ex_Rev3_Results/" #2 older
#     dataGSDFile = inputFolder + "GSD_FinalResultsData_2018v2.dat"  #Source: MODIS Data #v2 modified with different days
# elif year == 2020:
#     #Field Data from Input (2020)
#     inputFolder = "./Input/SeaIceWGJ2020/"
#     inputFolder2 = "./Input/SeaIceWGJ2020/2020Ex_Rev1_Results/"
#     dataGSDFile = inputFolder + "GSD_FinalResults2020.dat" #Source: MODIS Data #v1 modified with different days 2020

# dataTempFile = inputFolder + "Temp_data.dat"  #Source: Copernicus Data
# dataThickFile = inputFolder + "Thick_data.dat"  #Source: Copernicus Data
# dataConcV            = np.loadtxt(dataGSDFile)
# dataTempV            = np.loadtxt(dataTempFile)
# dataThickV            = np.loadtxt(dataThickFile)

# Output contact and updated position file
figDir              = mainVfolder
figDir4             = mainOutputfolder
figDir5             = mainVfolder+"GSD/"
# figConcDir          = mainVfolder+"FluxOnly/"
# figDir6             = mainVfolder+"LevelSet0/"
# figDir7             = mainVfolder+"OceanGrid/"
figDir8             = mainVfolder+"Output_Comparison/"
figDir9             = mainVfolder+"Output_Comparison/GSD_Compare/"


# posFile            =  mainOutputfolder+"positions0.dat"
# thickFile          =  mainOutputfolder+"thickness0.dat"

# LSFile             =  mainOutputfolder+"samplels0.dat"

# numbergFile        =  mainOutputfolder+"numberg0.dat"
# DiamFile           =  mainOutputfolder+"Diameters0.dat"
# npergFile          = mainOutputfolder+"npointsperg0.dat"
# morphFile          =  inputFolder0 + "morphologies.dat"
# #morphFile         =  "./Input/mainCompression3/morphologies.dat"
# pointsFile         = mainOutputfolder+"pointsout0.dat"
# pointsLSFile         = mainOutputfolder+"pointsoutLS0.dat"
# concFile         = mainOutputfolder+"normConc0.dat"
# oceanGridFile         = mainOutputfolder+"oceanTGrid0.dat"
# oceanGridDimFile         = mainOutputfolder+"oceanTGridDim0.dat"
# oceanGridCoordFile         = mainOutputfolder+"oceanTGridCoord0.dat"
# testParamsFile             = mainOutputfolder+"testParams0.dat"

# #Plot fine points
# posFinesFile             = mainOutputfolder+"posFines0.dat"
# numberFinesFile             = mainOutputfolder+"numberFines0.dat"

#File output for comparing GSD
#fileImgDir = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/Papers/Sea Ice Gifs/Multi_spec/WG_nasa-worldview-2020-06-01T00 00Z-to-2020-06-30T00 00Z_label_edit222Jv2/" #2020
#For computer
#fileImgDir = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/Papers/Sea Ice Gifs/Multi_spec/2018ExResults/" #2018
#For HPC
# if year == 2018:
#     fileImgDir = "./Input/SeaIceWGJ2018/2018Ex_Rev3_Results/" #2 older  #2018 #Use information from Input folder related to MODIS image analysis
# elif year == 2020:
#     fileImgDir = "./Input/SeaIceWGJ2020/2020Ex_Rev1_Results/"  #2020
# outErrorDir = mainVfolder+"GSD_Error/"

#Only if it is created we make these new folder (otherwise not create)
if create_folder:
    if os.path.exists(figDir5) == False:
        os.mkdir(figDir5)
    # if os.path.exists(figConcDir) == False:
    #     os.mkdir(figConcDir)
    # if os.path.exists(figDir6) == False:
    #     os.mkdir(figDir6)
    # if os.path.exists(outErrorDir) == False:
    #     os.mkdir(outErrorDir)
    # if os.path.exists(figDir7) == False:
    #     os.mkdir(figDir7)
    if os.path.exists(figDir8) == False:
        os.mkdir(figDir8)
    if os.path.exists(figDir9) == False:
        os.mkdir(figDir9)

# #Walls are immutable
# #wallDir             = "./Input/grainsWallsFPlot/"
# wallDir             = "./Input/grainsWallsF/"
# wallposFile            =  "./Input/SeaIce/initialpositionsWallFrev.dat"
# wallmorphFile          =  "./Input/SeaIce/morphologiesWallF.dat"
# nWalls            = np.loadtxt(wallmorphFile, dtype = int)[0]
# posWall            = np.loadtxt(wallposFile)

# concV            = np.loadtxt(concFile)
# thickV            = np.loadtxt(thickFile)

# numberg            = np.loadtxt(numbergFile)
# nperg             = np.loadtxt(npergFile)
# DiamV            = np.loadtxt(DiamFile)
# testParamsV      = np.loadtxt(testParamsFile)
# posFinesV        = np.loadtxt(posFinesFile)
# numberFinesV     = np.loadtxt(numberFinesFile)

###########################################################################
#Turn on or off certain modules
oceanGridPlot = False #Plot grid setup
oceanGridPlot_img = False #Plot grid images
ocean_temp_grid_sample = False #For result output only
Thick_Output = False #Plot thickness of each coarse floe
Fine_Output = False #Output fine positions
fine_grad = False #Provide details of thickness changes
LS_plotting = 0 #TURN OFF, FOR MELT DEBUGGING
floe_plot = False #Plot floe images as patches
centroid_ave_cal = False #Calculate floe velocities for calibration
video_gen = False  #Generate video from floe images
###########################################################################

if oceanGridPlot:
    oceanGridV            = np.loadtxt(oceanGridFile)
    oceanGridDimV            = np.loadtxt(oceanGridDimFile)
    oceanGridCoordV            = np.loadtxt(oceanGridCoordFile)

    #Get ocean grid dimensions
    nOx = int(oceanGridDimV[0])
    nOy = int(oceanGridDimV[1])
    print("Ocean grid dimensions: ", nOx, nOy)

    #Get ocean grid coordinates
    xOGrid = np.zeros((nOx*nOy))
    yOGrid = np.zeros((nOx*nOy))
    ct_grid = 0
    for i in range(nOy):
        for j in range(nOx):
            xOGrid[j+i*nOx] = oceanGridCoordV[ct_grid]
            ct_grid += 1
            yOGrid[j+i*nOx] = oceanGridCoordV[ct_grid]
            ct_grid += 1

#plt.scatter(xOGrid, yOGrid)
#plt.show()

#print(oceanGridV)
#print((oceanGridV[3]))
#print(len(oceanGridV[3]))
#print(len(oceanGridV[2]))
#print(len(oceanGridV))
#print(len(xOGrid))
#exit(1)

# nGrains            = np.loadtxt(morphFile, dtype = int)[0]
# posRot             = np.loadtxt(posFile)

if LS_plotting == 1:
    LSVec              = np.loadtxt(LSFile)

if floe_plot:
    pointsDatOr          = np.loadtxt(pointsFile)
# pointsLS          = np.loadtxt(pointsLSFile)
nSteps             = 49
print(nSteps)
#posRot             = np.split(posRot,nSteps)
#pointsDat          = np.split(pointsDat,nSteps)
#We need to split unevently just in case

#npoints = 100 IDEALLY
#npointsWall = 700 #THIS ONE CAN BE FIXED

if floe_plot:
    count = 0
    npergT = []
    numberFinesT = []
    for i in range(nSteps):
        npergStep = []
        #numberFinesStep = []
        for j in range(int(numberg[i])):
            npergStep.append(int(nperg[count]))
            #numberFinesStep.append(int(numberFinesV[count]))
            count = count + 1
        npergT.append(npergStep)
        #numberFinesT.append(numberFinesStep)
        numberFinesT.append(int(numberFinesV[i]))
    #print("Number for nperg: ", count)
    #print("A: ", nperg)
    #print("B: ", npergT[99][0])

if floe_plot:
    pointsDat = []
    thickDat = []
    positDat = []
    oceanGridDat = []
    dot = 0 #to go through all point file pointsDatOr
    dot2 = 0 #to go through all data
    for i in range(nSteps):
        pointsDatT = []
        thickDatT = []
        positDatT = []
        #print(numberg[i])
        for j in range(int(numberg[i])):
            pointsDatP = []
            for k in range(int(npergT[i][j])):
                vec = [float(pointsDatOr[dot][0]), float(pointsDatOr[dot][1])]
                dot += 1
                pointsDatP.append(vec)
            pointsDatT.append(pointsDatP)
            thickDatT.append(thickV[dot2])
            positDatT.append(posRot[dot2][:2])
            dot2 += 1
        pointsDat.append(pointsDatT)
        thickDat.append(thickDatT)
        positDat.append(positDatT)

    #Same but for fines
    dot2 = 0
    finepositDat = []
    for i in range(nSteps):
        finepositDatT = []
        for j in range(int(numberFinesV[i])):
            finepositDatT.append(posFinesV[dot2])
            dot2 += 1
        finepositDat.append(finepositDatT)

#print(pointsLS)
dot2 = 0 #to go through all data

if LS_plotting == 1:
    npLS = 100
    pointsLS_step = []
    ct_ptsLS = 0

    adj = -1
    
    szx = int(LSVec[0])
    szy = int(LSVec[1])
    LS_step = []
    jls = 2 #Start information at 2, 0 and 1 were dimension sizes
    for i in range(nSteps):
        used_grain = 242
        refC = positDat[i][used_grain]
        LS_now = np.zeros((szy,szx))
        pointsLS_now = np.zeros((npLS,2))
        for j in range(szy):
            for k in range(szx):
                LS_now[j][k] = LSVec[jls]  #DONT USE DURING BREAKAGE, IT'S FOR DEBUGGING MELTING
                jls = jls + 1
        for jj  in range(npLS):
            pointsLS_now[jj] = [pointsLS[ct_ptsLS][0]-refC[0]+0.5*szx+adj, pointsLS[ct_ptsLS][1]-refC[1]+0.5*szy+adj] #Equation to convert points to LS grid coordinates and viceversa
            ct_ptsLS += 1
        LS_step.append(LS_now)
        pointsLS_step.append(pointsLS_now)

#print('See')
#print((pointsDatT))
#print('See2')
#print((pointsDat[0]))

#print('Size')
#print(dot)
#print(len(pointsDatOr))
#posRot0            = np.loadtxt(posFile)



    #if not os.path.exists(figDir):
    #    os.mkdir(figDir)

    # Read grain morphology data
#morphID = np.loadtxt(morphFile, dtype=int, skiprows=1)

    # Instantiate grains
#grains = np.empty(nGrains, dtype=object)
#for n in range(nGrains):
#    propFile =  "../input/grainproperty" + str(morphID[n])+ ".dat"
#    #print(propFile)
#    grains[n] = Grain(n, morphID[n], propFile, posRot0[n])

    # Plot range
#offset_val = 4000
#lX = offset_val #80000
#lY = offset_val #80000
#
#x1 = 0
#x2 = offset_val #80000
#y1 = 0
#y2 = offset_val #80000

if year == 2020:
    offset_val = 380 #2020
elif year == 2018:
    offset_val = 400 #2018

lX = offset_val #80000
lY = offset_val #80000

x1 = 0
x2 = offset_val #80000
y1 = 0
y2 = offset_val #80000

# #if not (makePics):
# #nSteps = 0

# polynew1 = [[-10,5100],[-10,2400],[100,2500],[500,2500],[600,2650],[1000,3100],[1500,3500],[1800,4000],[2000,4300],[2100,4900],[1900, 5150]]

# polynew2 = [[6000,4500],[6000,4800],[4000,3800],[3000,3500],[2500,3000],[2000,2800],[2000,2800],[2000,2000],[2100,1500],[2000,1000],[6000, 1000]]

def find_pos_ave(positDat):
    ave_posx = 0
    ave_posy = 0
    for i in range(len(positDat)):
        ave_posx += positDat[i][0]
        ave_posy += positDat[i][1]
    return ave_posx/len(positDat), ave_posy/len(positDat)

# for step in range(nSteps):
def write_img(step):
    if floe_plot:
        # Set up the figure
        #print(step)
        fig, ax = plt.subplots(figsize = (10,10))
        #ax.set_xlim(-lX,lX)
        #ax.set_ylim(-lY,lY)
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
        ax.autoscale_view()
                    
        # Update grain positions
        #posRotStep = posRot[step]
        #print(posRotStep)
    #    for n in range(nGrains):
    #        grains[n].updateXY(posRotStep[n])
            
        # Collect in patches
        patches = []
        pointsDatstep = pointsDat[step]
        print(len(pointsDatstep[0]))
        #print("PtsBIG: ",(pointsDatstep))
        #print("Pts0: ", len(pointsDatstep))
        #print("Npts0: ", int(numberg[step]))
        #pointIndiv = np.split(pointsDatstep, int(numberg[step]))
        
        #pointIndiv = []
        #jj = 0
        #for i in range(int(numberg[step])):
            #pointIndivSmall = []
            #print(npergT[step][i])
            #print(len(pointsDatstep[i]))
        #    pointsDatP = pointsDatstep[i]
        #    print(pointsDatP)
        #    for j in range(npergT[step][i]):
        #        #print(pointsDatstep[j])
        #        pointIndivSmall.append(pointsDatstep[j])
        #        jj += 1
        #    pointIndiv.append(pointIndivSmall)
        
        #print("PtsBIG: ",(pointIndiv[0]))
        print(step)
        #print("npts00: ", len(pointsDat[step]))
        #print("npts0: ", len(pointIndiv))
        print(int(numberg[step]))
        #print("npts2: ", len(pointIndiv[0]))
        
    #    #Check intersection with each other
    #    for i in range(int(numberg[step])-1):
    #        j = i+1
    #        for (j) in range(int(numberg[step])):
    #            p1 = shapely.geometry.Polygon(pointsDatstep[i])
    #            p2 = shapely.geometry.Polygon(pointsDatstep[i+1])
    #            if (p1.intersects(p2) == True):
    #                print("We got intersection")
    #            else:
    #                print("No intersection")
    
        patches_sel = []
        patches_sel2 = []
        for n in range(int(numberg[step])):
            #print(grains[n]._points)
            #poly = Polygon(pointIndiv[n], True)
            poly = Polygon(pointsDatstep[n], True)
            patches.append(poly)
            #if (1 < 0):  #Not really used with this
            if (n == 0):
            #if (n == 40 or n == 41 or n == 42 or n == 43 or n == 44 or n == 45):
                poly = Polygon(pointsDatstep[n], True)
                patches_sel.append(poly)
            if (n == 1):
                poly = Polygon(pointsDatstep[n], True)
                patches_sel2.append(poly)
    
            
        # pointsWall = []
        # for i in range (nWalls):
        #     positionP = posWall[i]
        #     #print(positionP)
        #     ng = np.loadtxt(wallmorphFile, dtype = int)[i+1]
        #     pointsFile = wallDir + "grainproperty" + str(ng) + ".dat"
        #     pointsExtract = open(pointsFile,"r")
        #     pointList = pointsExtract.readlines()[4]
        #     pointList = pointList.split()
        #     pointVec = []
        #     jj = 0
        #     for j in range(npointsWall):
        #         floatP = [float(pointList[jj]), float(pointList[jj+1])]
        #         Xg = floatP[0] + positionP[0]
        #         Yg = floatP[1] + positionP[1]
        #         vec = [Xg,Yg]
        #         pointVec.append(vec)
        #         jj = jj + 2
        #     pointsWall.append(pointVec)
    
        # patcheswall = []
        # for i in range(nWalls):
        #     poly = Polygon(pointsWall[i], True)
        #     #print(pointsWall[i])
        #     patcheswall.append(poly)
            
        pCol = PatchCollection(patches, facecolors='white',edgecolors='black', lw=0.1)
        
        if fine_grad == False:
            ax.set_facecolor('xkcd:blue')
        else:
            if year == 2018:
                ax.set_facecolor(color = (0.1, 0.2, 0.99*step/49, 0.99*step/49))
            elif year == 2020:
                ax.set_facecolor(color = (0.1, 0.2, 0.99*step/25, 0.99*step/25))
    
        
        ax.add_collection(pCol)
        # pCol2 = PatchCollection(patcheswall, facecolors='black',edgecolors='black', lw=0.1)
        # ax.add_collection(pCol2)
        pCol_sel = PatchCollection(patches_sel, facecolors='red',edgecolors='black', lw=0.1)
        ax.add_collection(pCol_sel)
        pCol_sel2 = PatchCollection(patches_sel2, facecolors='pink',edgecolors='black', lw=0.1)
        ax.add_collection(pCol_sel2)
        
        # patcheswall2 = []
        # polynewP = Polygon(polynew1, True)
        # patcheswall2.append(polynewP)
        # pCol3 = PatchCollection(patcheswall2, facecolors='black',edgecolors='black', lw=0.1)
        # ax.add_collection(pCol3)
        
        # patcheswall3 = []
        # polynewPP = Polygon(polynew2, True)
        # patcheswall3.append(polynewPP)
        # pCol4 = PatchCollection(patcheswall3, facecolors='black',edgecolors='black', lw=0.1)
        # ax.add_collection(pCol4)
        
        # Setup up patchCollection
    #        if (n==3):
    #            print("ICE")
    #            pCol = PatchCollection(patches, facecolors='white',edgecolors='black', lw=0.1)
    #        else:
    #            print("LAND")
    #            pCol2 = PatchCollection(patches, facecolors='green',edgecolors='black', lw=0.1)
    
    
        #ax.add_collection(pCol2)
        
        #Add thickness text
        if Thick_Output:
            for ith in range(int(numberg[step])):
                if step == 0:
                    plt.text(positDat[step][ith][0], positDat[step][ith][1], str(round(thickDat[step][ith]*1,2)), size = 10)
                else:
                    plt.text(positDat[step][ith][0], positDat[step][ith][1], str(round(thickDat[step][ith]*1,2)), size = 10)
               
        #Plot fines as points
        if Fine_Output:
            fpList = finepositDat[step] #Customized for each step based on number of fines
            xs = [x[0] for x in fpList]
            ys = [x[1] for x in fpList]
            plt.scatter(xs, ys) #plt.plot??
    
        #plt.axis('off')
        plt.savefig(figDir + 'step_' + ('%d' % step) + '.png', format='png')
        plt.close(fig)
    
    if LS_plotting == 1:
    #    #Saving LS Image
    #    fig3 = plt.figure(figsize=(10,10))
    #    axx = fig3.add_subplot(1, 1, 1)
    #    #matplotlib.image.imsave(figDir6 + '/step_' + str(step) + '.png', LS_step[step])
    #    axx.imshow(LS_step[step], origin='lower')
    #    axx.contour(LS_step[step], levels=[0], colors='white', alpha=0.5)
    #    #axx.contour(LS_step[step], levels=np.logspace(-5, 0, 10), colors='white', alpha=0.5)
    #    plt.savefig(figDir6 + '/step_' + ('%d' % step) + '.png', format='png')
    #    plt.close(fig3)
    #    ##im = Image.fromarray(LS_step[step])
    #    ##print(im)
    #    ##im.save(figDir6 + '/step_' + str(step) + '.jpg')

        #Saving LS Image
        vmin = -20
        vmax = 30
        fig3 = plt.figure(figsize=(10,10))
        fig3.add_subplot(1, 1, 1)
        #matplotlib.image.imsave(figDir6 + '/step_' + str(step) + '.png', LS_step[step])
        cb = plt.imshow(LS_step[step], origin='lower')
        plt.contour(LS_step[step], levels=[0], colors='white', alpha=0.5)
        plt.colorbar(cb)
        #plt.clim(vmin,vmax)
        #axx.contour(LS_step[step], levels=np.logspace(-5, 0, 10), colors='white', alpha=0.5)
        
        datapLS = pointsLS_step[step]
        x, y = datapLS.T
        plt.scatter(x,y)
        
        plt.savefig(figDir6 + '/step_' + ('%d' % step) + '.png', format='png')
        plt.close(fig3)
        ##im = Image.fromarray(LS_step[step])
        ##print(im)
        ##im.save(figDir6 + '/step_' + str(step) + '.jpg')
        
    #Plot Ocean Grid
    if oceanGridPlot_img:
        figGrid = plt.figure(figsize=(10,10))
        plt.scatter(xOGrid, yOGrid, c = oceanGridV[step])
        #print(xOGrid, len(xOGrid))
        #print(yOGrid)
        #print(oceanGridV[step])
        plt.title('Ocean Temperature at step: ' + str(step))
        plt.colorbar()
        plt.clim(-2.0,4.0)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig(figDir7 + '/step_' + ('%d' % step) + '.png', format='png')
        plt.close(figGrid)
        
    # #Temporary plot for explaining method
    if ocean_temp_grid_sample:        
        z_x0 = 58
        z_x1 = 136
        z_y0 = 111
        z_y1 = 153
        ct = 0
        z_xOGrid = []
        z_yOGrid = []
        z_oceanGridV = []
        for ii in range(len(xOGrid)):
            if xOGrid[ii] >= z_x0 and xOGrid[ii] <= z_x1 and yOGrid[ii] >= z_y0 and yOGrid[ii] <= z_y1:
                z_xOGrid.append(xOGrid[ii])
                z_yOGrid.append(yOGrid[ii])
                z_oceanGridV.append(oceanGridV[step][ii]) 
                ct += 1
        
        #print((ct))
        #print(z_y1-z_y0+1)
        #print(z_x1-z_x0+1)
        ct2 = 0
        TempLS_Matrix = np.zeros((z_y1-z_y0+1, z_x1-z_x0+1))
        for i in range(z_y1-z_y0+1): #Y
             for j in range(z_x1-z_x0+1): #X
                TempLS_Matrix[i][j] = z_oceanGridV[ct2]
                ct2 += 1
        
        #z_oceanGridV = np.zeros((offset_val,offset_val))
        
        # for i in range(offset_val): #Y
        #     for j in range(offset_val): #X
        #         z_oceanGridV[offset_val-1-i][offset_val] = oceanGridV[step][ct]
        #         ct += 1
                
        # TempLS_Matrix = np.zeros((z_y1-z_y0+1, z_x1-z_x0+1))
        # for i in range(z_y1-z_y0+1):  #Y
        #     for j in range(z_x1-z_x0+1):  #X
        #         TempLS_Matrix[i][j] = z_oceanGridV[z_y0 + i][z_x0 + j]
        # figGrid = plt.figure(figsize=(10,5))
        # plt.scatter(z_xOGrid, z_yOGrid, c = z_oceanGridV)
        # #plt.imshow(z_oceanGridV)
        # clb=plt.colorbar()
        # plt.clim(-2.0,4.0)
        # clb.ax.tick_params(labelsize=10) 
        # clb.ax.set_title('Temperature (Celsius)',fontsize=12)
        # plt.title("Ocean Temperature Grid", fontsize=18)
        # plt.xlabel("Grid x units", fontsize=15)
        # plt.ylabel("Grid y units", fontsize=15)
        # plt.savefig(figDir7 + 'Temp_Matrix0_step_'+str(step)+'.png', format='png')
        # plt.close()
        
        figGrid = plt.figure(figsize=(10,5))
        plt.imshow(TempLS_Matrix)
        clb=plt.colorbar()
        plt.clim(-2.0,2.0)

        clb.ax.tick_params(labelsize=15) 
        clb.ax.set_title('Temperature (\N{DEGREE SIGN}C)',fontsize=17)
        plt.title("Ocean Temperature Grid", fontsize=20)
        plt.xlabel("Distance in x (km)", fontsize=18)
        plt.ylabel("Distance in y (km)", fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(figDir7 + 'Temp_Matrix_step_'+str(step)+'.png', format='png')
        plt.close()


#Parallel
#a_pool = Pool()
#a_pool.map(write_img, range(nSteps))

#In order for GSD
gsdFileMaxS = figDir4 + "GrainSizeDist_iter0_0.dat"
gsdFile0 = figDir4 + "GrainSizeDist_iter0_0.dat"
gsdFileMinS = figDir4 + "GrainSizeDist_iter" + str(nSteps-1) + "_0.dat"
gsd0Max = np.loadtxt(gsdFileMaxS)
gsd0Init = np.loadtxt(gsdFile0)
gsd0Min = np.loadtxt(gsdFileMinS)
MaxR = gsd0Max[0][0]*1.05
MinR = gsd0Min[len(gsd0Min)-1][0]

MaxR = 55 #42.84503 
MinR = 1.70271

#Just ranges for plotting for the limits of the x axis
print("Max Size: ", MaxR)
print("Min Size: ", MinR)
if MinR == 0:
    MinR = 0


gsd0Init = np.loadtxt(gsdFile0)

sizeV = []
countV = []
for i in range(len(gsd0Init)-1):
    sizeV.append(gsd0Init[i][0])
    countV.append(gsd0Init[i][2])


#Evolve bin functions

gsdC_Num = np.zeros((nSteps, len(sizeV)))
gsdCumC_Num = np.zeros((nSteps, len(sizeV)))
bin_size = sizeV[0] - sizeV[1]


#Function of Probabililty of getting sizeS from sizeL, sizeL > sizeS
def Bprob(sizeL, sizeS):
    probBh = 0
    probSh = 0
    #Linearly the floe splits into a big half equal to sizeL/2 + z * sizeL/2 and small half equal to sizeL/2 -z * sizeL/2
    # z = 0 ~ 0.8 * sizeL/2
    #If one or both floes are within the range of sizeS, then it is fine. 
    #Hence the probablity is the sum of both halves such that the length of the segment is between sizeS and sizeS - bin_size

    #Fine pass #Smaller Half only
    if sizeS == sizeL:
        probSh = ( (0.5*sizeL*0.8) - (0.5*sizeL - (sizeL - bin_size)) ) / (0.5*sizeL*0.8)
        return probBh + probSh

    #Bigger Half
    if sizeS > 0.5*sizeL:
        probBh = bin_size / (0.5*sizeL*0.8)
    #Smaller Half
    else:
        probSh = bin_size / (0.5*sizeL*0.8)
    
    return probBh + probSh

def Mprob(sizeL, DLostacc):
    if DLostacc > sizeL:
        return 1
    elif DLostacc < sizeV[-1] - bin_size:
        return 0
    else:
        return max(DLostacc, bin_size)*0.50 / bin_size

#Initialize
gsdC_Num[0] = countV

#Define parameters
Tadjust = 0.7
Lthermal = 2000 #Vert, Lateral rate ##km of side lost given km of thickness lost
qv = 25 #W/m2 Celsius
rho_ice = 970 #kg/m3
Lf = 330000 #J/kg
kappa = 20 #m2/s
Qf = 310 #W/m2 
convert_f = 86400/1000 #To convert m/2 to km/day
deltaT = 2 #Dif. between ocean and ice temperature for melt
thick_loss_rate = (qv / (rho_ice*Lf)) * deltaT * convert_f  #How many km of ice are lost per day
lat_loss_rate = Lthermal * thick_loss_rate  #How many km of ice are lost given thickness loss
DLostacc = 0

Brate = 86400 / 1000
Badjust = 0.0035  #0.0035 alone

for day in range(1,nSteps,1):
    #gsdC_Num[day] = gsdC_Num[0] #For test only
    GainM = 0
    LossM = 0
    GainB = 0
    LossB = 0

    
    for j in range(len(sizeV)):
        # #Melt (immediate bin interaction)
        # if j == 0:
        #     GainM = 0
        # else:
        #     GainM = (2/sizeV[j-1]**2) * Lthermal * thick_loss_rate * gsdC_Num[day-1][j-1] * Tadjust
        # LossM = -(2/sizeV[j]**2) * Lthermal * thick_loss_rate * gsdC_Num[day-1][j] * Tadjust

        #Melt (immediate bin interaction)
        if j == 0:
            GainM = 0
        else:
            #GainM = (1/sizeV[j-1]) * lat_loss_rate * gsdC_Num[day-1][j-1] * Mprob(sizeV[j-1], DLostacc) * Tadjust
        #LossM = -(1/sizeV[j]) * lat_loss_rate * gsdC_Num[day-1][j] * Mprob(sizeV[j], DLostacc) * Tadjust
            #GainM = (2/sizeV[j-1])  * gsdC_Num[day-1][j-1] * Mprob(sizeV[j-1], DLostacc) * Tadjust
        #LossM = -(2/sizeV[j])  * gsdC_Num[day-1][j] * Mprob(sizeV[j], DLostacc) * Tadjust
            GainM = (2/1)  * gsdC_Num[day-1][j-1] * Mprob(sizeV[j-1], DLostacc) * Tadjust
        LossM = -(2/1)  * gsdC_Num[day-1][j] * Mprob(sizeV[j], DLostacc) * Tadjust
        
        #Break (all versus all)
        if j == 0:   #Largest bin (no gain, only loss)
            GainB = 0
            LossB = 0
            for k in range(j+1, len(sizeV), 1):
                LossB += - Brate * Bprob(sizeV[j], sizeV[k]) * gsdC_Num[day-1][j] * Badjust
        elif j < len(sizeV)-1: #(gain and loss)
            GainB = 0
            LossB = 0
            for k in range(0, j, 1):
                GainB += Brate * Bprob(sizeV[j], sizeV[k]) * gsdC_Num[day-1][k] * Badjust
            for k in range(j+1, len(sizeV),1):
                LossB += - Brate * Bprob(sizeV[j], sizeV[k]) * gsdC_Num[day-1][j] * Badjust
        else: #(last bin, lost to fines, gains from the rest)
            GainB = 0
            LossB = 0
            for k in range(len(sizeV)-1):
                GainB +=  Brate * Bprob(sizeV[j], sizeV[k]) * gsdC_Num[day-1][k] * Badjust
            LossB += - Brate * Bprob(sizeV[j], sizeV[j]) * gsdC_Num[day-1][j] * Badjust #Bprob(sizeV[j], sizeV[k]) ~= 1 FOR MIN SIZE
            print("Rate of smallest size gainM: ", GainM)
            print("Rate of smallest size lossM: ", LossM)
            print("Rate of smallest size gainB: ", GainB)
            print("Rate of smallest size lossB: ", LossB)
            print("NET: ", GainM + LossM + GainB + LossB)

        #Exact    
        #gsdC_Num[day][j] = max( (gsdC_Num[day-1][j] + GainM + LossM + GainB + LossB) , 0 )
        #Rounded
        #gsdC_Num[day][j] = max( (gsdC_Num[day-1][j] + round(GainM,0) + round(LossM,0) + round(GainB,0) + round(LossB,0)) , 0 )
        gsdC_Num[day][j] = max(  round( ( gsdC_Num[day-1][j] + (GainM + LossM + GainB + LossB) ), 0) , 0 )
        print("NET: ", GainM + LossM + GainB + LossB)

    # #Break
    # for j in range(len(sizeV)):
    #     if j == 0:
    #         Gain = 0
    #     else:
    #         for k in range(len(size)-j):
    #             Gain = break_rate * gsdC_Num[day-1][j-k+1] 
    #             Loss  = -break_rate * gsdC_Num[day-1][j+k] 
    #     gsdC_Num[day][j] = gsdC_Num[day-1][j] + Gain + Loss


    #Accumulate perimeter loss
    DLostacc += 2*lat_loss_rate
    print("Lateral Loss km:", DLostacc)
    #Make cumulative floe count for later
    ncumct = 0
    for j in range(len(sizeV)):
        ncumct += gsdC_Num[day][j]
        gsdCumC_Num[day][j] = ncumct
    



#CORE OF EVERYTHING


# #List of parameters
# nT = testParamsV[0]
# nTemp = testParamsV[1]
# nBreak = testParamsV[2]
# qvert = testParamsV[3]
# Khor = testParamsV[4]
# Qatm = testParamsV[5]
# Aatm = testParamsV[6]
# Batm = testParamsV[7]
# alpha_ice = testParamsV[8]
# alpha_ocean = testParamsV[9]
# limitMaxDiam = testParamsV[10]
# limitMinDiam = testParamsV[11]


# #For concentration variation
# varC = np.zeros((len(gsd0Max),nSteps))
# #For Total mass variation
# gsdMelt = np.zeros((len(gsd0Max),nSteps))
# gsdBreak = np.zeros((len(gsd0Max),nSteps))

# ##TEMPORARY FIX
# TotIceConcTemp = []
# for i in range(len(concV)):
#     TotIceConcTemp.append(concV[i][2])

#Vector for slope
t_sim = []
slope_gsd_sim1 = []
slope_gsd_sim2 = []

for i in range(nSteps):
    #write_img(i)
    
    #Generate GSD per Step
    gsdFile = figDir4 + "GrainSizeDist_iter" + str(i) + "_0.dat"
    gsd0 = np.loadtxt(gsdFile)
    gsdD = [] #Floe Size
    gsdP = [] #Cumulative Percent Pass
    gsdC = [] #Count
    gsdConc = [] #Concentration
    for j in range(len(gsd0)-1):
        #varC[j][i] = gsd0[j][3]/(offset_val*offset_val) #Must be 3
        gsdD.append(gsd0[j][0])
        gsdP.append(gsd0[j][1])
        gsdC.append(gsd0[j][2]) #Count
        gsdConc.append(gsd0[j][3])
        #gsdConc.append(gsd0[j][3]/(offset_val*offset_val))
        # gsdMelt[j][i] = gsd0[j][4]
        # gsdBreak[j][i] = gsd0[j][5]
    gsdCumC = []
    ncum = 0
    for j in range(len(gsd0)-1):
        ncum += gsdC[j]
        gsdCumC.append(ncum)

#        if j < len(gsd0) - 1:  #TEMP
#            gsdConc.append( (gsd0[j][1] - gsd0[j+1][1]) * TotIceConcTemp[i]  )
#            varC[j][i] =  (gsd0[j][1] - gsd0[j+1][1]) * TotIceConcTemp[i]
#        else:   #TEMP
#            gsdConc.append(  gsd0[j][1] * TotIceConcTemp[i] )
#            varC[j][i] =  gsd0[j][1] * TotIceConcTemp[i]
        
#    figG = plt.figure(figsize=(10,10))
#    plt.plot(gsdD,gsdP,'ro-')
#    plt.title('Evolution Grain Size Distribution')
#    plt.xscale('log')
#    #plt.yscale('log') #For log-log plot
#    plt.xlim(MinR,MaxR) #Must be log limit so avoid zeros
#    plt.ylim(0.0,100)  #Must be log limit so avoid zeros
#    plt.xlabel('Grain Size')
#    plt.ylabel('% Pass')
#    ax2 = plt.twinx()
#    ax2.plot(gsdD,gsdC,'b*-')
#    ax2.xlim(MinR,MaxR) #Must be log limit so avoid zeros
#    ax2.ylim(0.0,100)  #Must be log limit so avoid zeros
#    ax2.xlabel('Grain Size')
#    ax2.ylabel('# Floes')
#    ax2.ylim(0.0, 30)  #Regulate using floe Number
      
    #For easier plotting
    t_sim.append(i)


    # #Function to get slope of GSD Sim
    # slope_data1 = gsd_slope_linear(objective,gsdD,gsdP)  
    # slope_data2 = gsd_slope_exp(gsdD,gsdP)    
    # slope_gsd_sim1.append(slope_data1) 
    # slope_gsd_sim2.append(slope_data2)    


    fig, ax1 = plt.subplots(figsize=(10,10))
    #plt.figure(figsize=(10,10))

    color = 'tab:red'
    ###plt.xscale('log') ###Do for GSD only, not histogram
    plt.xlim(0, MaxR)
    ax1.set_xlabel('Grain Size')
    ax1.set_ylabel('Cumul. # Floes', color=color)
    ax1.plot(gsdD, gsdCumC, 'ro-')
    ax1.tick_params(axis='y', labelcolor=color)
    plt.ylim(0,500)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('# Floes', color=color)  # we already handled the x-label with ax1
    ax2.plot(gsdD, gsdC, 'b*-')
    ax2.bar(gsdD, gsdC, align='edge', color = 'blue', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.25) ) # A bar chart
    ax2.tick_params(axis='y', labelcolor=color)
    #ax2.bar(gsdD, gsdC, align='edge', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.5) ) # A bar chart
    plt.ylim(0,500)
 
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    #Revert X-Axis
    #axG = plt.gca()
    #axG.set_xlim(axG.get_xlim()[::-1])
    plt.title('Evolution Grain Size Distribution and Count')
    plt.savefig(figDir5 + '/step_' + ('%d' % i) + '.png', format='png')
    plt.close()


    fig, ax1 = plt.subplots(figsize=(10,10))
    #plt.figure(figsize=(10,10))

    color = 'tab:red'
    ###plt.xscale('log') ###Do for GSD only, not histogram
    plt.xlim(0, MaxR)
    ax1.set_xlabel('Grain Size')
    ax1.set_ylabel('Cumul. # Floes', color=color)
    ax1.plot(gsdD, gsdC_Num[i], 'ro-')
    ax1.tick_params(axis='y', labelcolor=color)
    plt.ylim(0,500)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('# Floes', color=color)  # we already handled the x-label with ax1
    ax2.plot(gsdD, gsdC, 'b*-')
    #ax2.bar(gsdD, gsdC, align='edge', color = 'blue', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.25) ) # A bar chart
    ax2.tick_params(axis='y', labelcolor=color)
    #ax2.bar(gsdD, gsdC, align='edge', width = -(gsdD[0]-gsdD[1])+2, fc=(0, 0, 1, 0.5) ) # A bar chart
    plt.ylim(0,500)
 
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    #Revert X-Axis
    #axG = plt.gca()
    #axG.set_xlim(axG.get_xlim()[::-1])
    plt.title('Evolution Count versus Numerical')
    plt.savefig(figDir5 + '/GSD_step_' + ('%d' % i) + '.png', format='png')
    plt.close()




#Calculate average velocities to control ocean drag
if centroid_ave_cal:
    velV = []
    meanx = np.zeros(nSteps)
    meany = np.zeros(nSteps)
    for i in range(nSteps):
        [meanx[i], meany[i]] = find_pos_ave(positDat[i])
    for i in range(nSteps-1):
        difx = meanx[i+1] - meanx[i]
        dify = meany[i+1] - meany[i]
        velV.append([i, ((difx**2 + dify**2)**0.5)/1, difx, dify, math.atan(dify/difx)]) #Time(day), Speed km/d, movex km, movey km, angle
    velV.append([nSteps-1, ((difx**2 + dify**2)**0.5)/1, difx, dify, math.atan(dify/difx)])
    np.savetxt(figDir+'velocity_ave.dat',velV, fmt='%d %.8f %.8f %.8f %.8f')


#exit(1)

#Plot Sim GSD vs. Field Data GSD
t_gsd=[] #Vector of days with field data
#Import GSD from field

if year == 2020:
    len_range = 25 #Consider Incomplete Days or Gaps #2020
elif year == 2018:
    len_range = 49 #Consider Incomplete Days or Gaps #2018

# for i in range(len_range+1):
#     gds_data_file = inputFolder2 + str(i) + "_GSD_Data.dat"
#     if (os.path.isfile(gds_data_file) != 0):
#         t_gsd.append(i)

# #Vector for slope
# slope_gsd_data1 = []
# slope_gsd_data2 = []

#Vector count slope
#slope_NF_data = []
slope_NF_sim = []
t_num = []

slope_NF_num = []

for i in range(nSteps): #Consider Incomplete Days or Gaps
    # if i > 470: 
    # #if i == 47:  #Should be needed
    #     gds_data_file = inputFolder2 + "48_GSD_Data.dat" #Before LS-DEM step modification
    # else:
    #     gds_data_file = inputFolder2 + str(i) + "_GSD_Data.dat"
    
    #All sim data
    gsdFileF = figDir4 + "GrainSizeDist_iter" + str(i) + "_0.dat"
    gsd0F = np.loadtxt(gsdFileF)
    log_gsdDF = []
    log_gsdNFF = []
    gsdNF_or = []

    #Do the same for numerical model for gsdC_Num
    log_gsdDFNum = []
    log_gsdNFFNum = []
    gsdNF_or_Num = []

    for j in range(len(gsd0F)-1):
        if gsd0F[j][0] <= 0:
            log_gsdDF.append(0.0)
        else:
            log_gsdDF.append(math.log(gsd0F[j][0]))
        if gsd0F[j][2] <= 0:
            gsdNF_or.append(0.0)
        else:
            gsdNF_or.append(gsd0F[j][2])
    
    for j in range(len(gsdC_Num[i])):
        if gsd0F[j][0] <= 0:
            log_gsdDFNum.append(0.0)
        else:
            log_gsdDFNum.append(math.log(gsd0F[j][0]))
        if gsdC_Num[i][j] <= 0:
            gsdNF_or_Num.append(0.0)
        else:
            gsdNF_or_Num.append(gsdC_Num[i][j])


    #Convert to cumulative sum
    #log_gsdNFFp = log_gsdNFF
    #log_gsdNFF = acc_sum(log_gsdNFFp)
    gsdNF_orp = gsdNF_or
    gsdNF_or = acc_sum(gsdNF_orp)

    gsdNF_orp_Num = gsdNF_or_Num
    gsdNF_or_Num = acc_sum(gsdNF_orp_Num)


    for idx in range(len(gsdNF_or)):
        if gsdNF_or[idx] <=0:
            log_gsdNFF.append(0.0)
        else:
            log_gsdNFF.append(math.log(gsdNF_or[idx]))
    
    for idx in range(len(gsdNF_or_Num)):
        if gsdNF_or_Num[idx] <=0:
            log_gsdNFFNum.append(0.0)
        else:
            log_gsdNFFNum.append(math.log(gsdNF_or_Num[idx]))
    
    #Filter zero results if occurring
    log_gsdDF0 = log_gsdDF
    log_gsdNFF0 = log_gsdNFF
    #[log_gsdDF, log_gsdNFF] = zeros_filter(log_gsdDF0, log_gsdNFF0)
    [log_gsdDF, log_gsdNFF] = zeros_filter_slope(log_gsdDF0, log_gsdNFF0)
    t_num.append(i)
    slope_simNO  = gsd_slope_linear(objective, log_gsdDF, log_gsdNFF) 
    slope_NF_sim.append(abs(slope_simNO))   

    #FOr numerical PDE results
    log_gsdDF0Num = log_gsdDFNum
    log_gsdNFF0Num = log_gsdNFFNum
    [log_gsdDFNum, log_gsdNFFNum] = zeros_filter_slope(log_gsdDF0Num, log_gsdNFF0Num)
    slope_numNO  = gsd_slope_linear(objective, log_gsdDFNum, log_gsdNFFNum) 
    slope_NF_num.append(abs(slope_numNO))   
    
    #For comparison
    #if (os.path.isfile(gds_data_file) != 0):
    if 1 == 1:

        # #Field Data
        # gsd_dataV = np.loadtxt(gds_data_file)
        # gsd_dataD = [] #Floe Size
        # gsd_dataP = [] #Cumulative Percent Pass
        # gsd_dataNF = [] # of Floes Percent Pass
        # #gsdConc = [] #Concentration
        
        # log_gsd_dataD = []
        # log_gsd_dataNF = []
        # for j in range(len(gsd_dataV)):
        #     gsd_dataD.append(gsd_dataV[j][0])
        #     gsd_dataP.append(gsd_dataV[j][1])
        #     gsd_dataNF.append(gsd_dataV[j][2])
        #     if gsd_dataV[j][0] <= 0:
        #         log_gsd_dataD.append(0.0)
        #     else:
        #         log_gsd_dataD.append(math.log(gsd_dataV[j][0]))
        #     #if gsd_dataV[j][2] <= 0:
        #     #    log_gsd_dataNF.append(0.0)
        #     #else:
        #     #    log_gsd_dataNF.append(math.log(gsd_dataV[j][2]))

        # #Convert to cumulative sum
        # #log_gsd_dataNFp = log_gsd_dataNF
        # #log_gsd_dataNF = acc_sum(log_gsd_dataNFp)
        # gsd_dataNFp = gsd_dataNF
        # gsd_dataNF = acc_sum(gsd_dataNFp)
        
        # for j in range(len(gsd_dataNF)):
        #     if gsd_dataNF[j] <=0:
        #         log_gsd_dataNF.append(0.0)
        #     else:
        #         log_gsd_dataNF.append(math.log(gsd_dataNF[j]))
    
        # #Function to get slope of GSD Data
        # slope_data1 = gsd_slope_linear(objective, gsd_dataD,gsd_dataP) 
        # slope_data2 = gsd_slope_exp(gsd_dataD,gsd_dataP)    
        # slope_gsd_data1.append(slope_data1)
        # slope_gsd_data2.append(slope_data2)
        
        #Sim Data
        gsdFile = figDir4 + "GrainSizeDist_iter" + str(i) + "_0.dat"
        gsd0 = np.loadtxt(gsdFile)
        gsdD = [] #Floe Size
        gsdP = [] #Cumulative Percent Pass
        gsdNF = [] #Count
        #gsdConc = [] #Concentration
        
        log_gsdD = []
        log_gsdNF = []
        for j in range(len(gsd0)-1):
            gsdD.append(gsd0[j][0])
            gsdP.append(gsd0[j][1])
            gsdNF.append(gsd0[j][2])
            if gsd0[j][0] <= 0:
                log_gsdD.append(0.0)
            else:
                log_gsdD.append(math.log(gsd0[j][0]))
            # if gsd0[j][2] <= 0:
            #     log_gsdNF.append(0.0)
            # else:
            #     log_gsdNF.append(math.log(gsd0[j][2]))


        #Convert to cumulative sum
        #log_gsdNFp =  log_gsdNF
        #log_gsdNF = acc_sum(log_gsdNFp)
        gsdNFp = gsdNF
        gsdNF = acc_sum(gsdNFp)

        for j in range(len(gsdNF)):
            if gsdNF[j] <=0:
                log_gsdNF.append(0.0)
            else:
                log_gsdNF.append(math.log(gsdNF[j]))

        #print("gsdD: ", gsdD)
        #print("gsdNF: ", gsdNF)
        #print("LOG gsdD: ", log_gsdD)
        #print("LOG gsdNF: ", log_gsdNF)
        print(i)
        print("gsdD: ", gsdD)
        print("gsdNF: ", gsdC_Num[i])
        # print("LOG gsdD: ", log_gsdD)
        # print("LOG gsdNF: ", log_gsdNF)
        
        #Filter original data to remove zero bins (remove flat)
        # gsd_dataD0 = gsd_dataD
        # gsd_dataNF0 = gsd_dataNF
        # log_gsd_dataD0 = log_gsd_dataD
        # log_gsd_dataNF0 = log_gsd_dataNF
        gsdD0 = gsdD
        gsdNF0 = gsdNF
        log_gsdD0 = log_gsdD
        log_gsdNF0 = log_gsdNF
        
        # [gsd_dataD, gsd_dataNF] = zeros_filter(gsd_dataD0, gsd_dataNF0)
        # [log_gsd_dataD, log_gsd_dataNF] = zeros_filter(log_gsd_dataD0, log_gsd_dataNF0)
        [gsdD, gsdNF] = zeros_filter(gsdD0, gsdNF0)
        [log_gsdD, log_gsdNF] = zeros_filter(log_gsdD0, log_gsdNF0)
        
        
        # [log_gsd_dataDslope, log_gsd_dataNFslope] = zeros_filter_slope(log_gsd_dataD0, log_gsd_dataNF0)
        #Function to get slope of Log Curve of Count Data
        # slope_dataNO = gsd_slope_linear(objective, log_gsd_dataDslope, log_gsd_dataNFslope) 
        #slope_simNO  = gsd_slope_linear(objective, log_gsdD, log_gsdNF) 
        # slope_NF_data.append(abs(slope_dataNO))
        #slope_NF_sim.append(abs(slope_simNO))   


        #Plot figure
        figg = plt.figure(figsize=(10,10))
        plt.xlabel('Diameter')
        plt.ylabel('Percent Pass (%)')

        # #Calculate RSME (assume all data has same length)
        # Pass_RMSE = 0
        # for ii in range(len(gsd_dataD0)):  #Use original value here
        #     Pass_RMSE += (gsd_dataP[ii]-gsdP[ii])**2
        # Pass_RMSE /= len(gsd_dataD)
        # Pass_RMSE = Pass_RMSE**0.5

        # # plt.plot(gsd_dataD0, gsd_dataP, 'r', marker='s', fillstyle='none', label='Sea Ice FSD: Data') #Same Size Data #Use all data points
        # plt.plot(gsdD0, gsdP, 'b*', label='Sea Ice FSD: Simulation RMSE = '+str(round(Pass_RMSE,4)) )
        # plt.legend()
        # plt.title('Floe Size Distribution Simulation and MODIS Data, Step: ' + str(i))
        # plt.savefig(figDir9 + "SeaIce_Concentration_Compare_step_"+str(i)+".png", format='png')
        # plt.close()
        
        # plt.plot(gsd_dataD, gsd_dataNF, 'r', marker='s', fillstyle='none', label='Sea Ice Floes Count Distrib.: Data') #Same Size Data
        plt.plot(gsdD, gsdNF, 'b*-', label='Sea Ice Floe Count Distrib.: Simulation')
        plt.legend()
        plt.title('Floe Count Distribution Simulation and MODIS Data, Step: ' + str(i))
        plt.savefig(figDir9 + "SeaIce_NoFloes_Compare_step_"+str(i)+".png", format='png')
        plt.close()
        
        # plt.plot(gsd_dataD, gsd_dataNF, 'r', marker='s', fillstyle='none', label='Sea Ice Floes Count Distrib.: Data') #Same Size Data
        plt.plot(gsdD, gsdNF, 'b*-', label='Sea Ice Floe Count Distrib.: Simulation')
        plt.legend()
        plt.yscale('log')
        plt.xscale('log')
        plt.title('Floe Count Distribution Simulation and MODIS Data LOG Axis, Step: ' + str(i))
        plt.savefig(figDir9 + "SeaIce_NoFloesLOGAxis_Compare_step_"+str(i)+".png", format='png')
        plt.close()
        
        # plt.plot(log_gsd_dataDslope, log_gsd_dataNFslope, 'r', marker='s', fillstyle='none', label='Sea Ice Floes Count Distrib.: Data') #Same Size Data
        plt.plot(log_gsdD, log_gsdNF, 'b*--', label='Sea Ice Floe Count Distrib.: Simulation')
        plt.plot(log_gsdDFNum, log_gsdNFFNum, 'rs-', label='Sea Ice Floe Count Distrib.: Analytical')
        plt.legend()
        plt.title('Floe Count Distribution Simulation and MODIS Data LOG, Step: ' + str(i))
        plt.savefig(figDir9 + "SeaIce_NoFloesLOG_Compare_step_"+str(i)+".png", format='png')
        plt.close()
        
        
#Compare Count Slope

figCount = plt.figure(figsize=(10,10))
plt.xlabel('Time (Days)')
plt.ylabel('FSD Slope')
#plt.plot(t_gsd, slope_NF_data, 'r', marker='s', fillstyle='none', label='Sea Ice F. Count Slope: Data') #Same Size Data
plt.plot(t_num, slope_NF_sim, 'b*-', label='Sea Ice F. Count Slope: Sim')
plt.plot(t_num, slope_NF_num, 'r*-', label='Sea Ice F. Count Slope: Num')
plt.legend()
plt.title('Evolution of Sea Ice Count Slope')
plt.savefig(figDir9 + 'SeaIce_NFloesCount_Slope_Evol.png', format='png')
plt.close()

#Plot Slope Vector
vec_text_out_data = []
vec_text_out_sim = []
# for i in range(len(t_gsd)):
#     vec_text_out_data.append((t_gsd[i],slope_NF_data[i]))
for i in range(len(t_num)):
    vec_text_out_sim.append((t_num[i],slope_NF_sim[i]))

# np.savetxt(figDir9+'Slope_data.dat', vec_text_out_data, fmt='%.8f %.8f')
np.savetxt(figDir9+'Slope_sim.dat', vec_text_out_sim, fmt='%.8f %.8f')


#Compare slopes Plot (Linear)
#Plot figure
figGSDS = plt.figure(figsize=(10,10))
plt.xlabel('Time (Days)')
plt.ylabel('FSD Slope')

# #Calculate RSME (assume all data has same length)
# Slope_RMSE = 0
# for i in range(len(slope_gsd_data1)):
#     tstepD = t_gsd[i]
#     for j in range(len(slope_gsd_sim1)):
#         if (tstepD == t_sim[j]):
#             indexstepN = j
#             break
#     Slope_RMSE += (slope_gsd_data1[i]-slope_gsd_sim1[indexstepN])**2
# Slope_RMSE /= len(slope_gsd_data1)
# Slope_RMSE = Slope_RMSE**0.5

# plt.plot(t_gsd, slope_gsd_data1, 'r', marker='s', fillstyle='none', label='Sea Ice FSD Linear Slope: Data') #Same Size Data
# plt.plot(t_sim, slope_gsd_sim1, 'r.', label='Sea Ice FSD Linea Slope: Simulation RMSE = '+str(round(Slope_RMSE,4)) )
# plt.legend()
# plt.title('Evolution of Sea Ice FSD Slope (Linear)')
# plt.savefig(figDir9 + 'SeaIce_FSD_Slope_Compare_Linear.png', format='png')
# plt.close()

# #Compare slopes Plot (Exponential)
# #Plot figure
# figGSDS2 = plt.figure(figsize=(10,10))
# plt.xlabel('Time (Days)')
# plt.ylabel('FSD Slope')

# #Calculate RSME (assume all data has same length)
# Slope_RMSE = 0
# for i in range(len(slope_gsd_data2)):
#     tstepD = t_gsd[i]
#     for j in range(len(slope_gsd_sim2)):
#         if (tstepD == t_sim[j]):
#             indexstepN = j
#             break
#     Slope_RMSE += (slope_gsd_data2[i]-slope_gsd_sim2[indexstepN])**2
# Slope_RMSE /= len(slope_gsd_data2)
# Slope_RMSE = Slope_RMSE**0.5

# plt.plot(t_gsd, slope_gsd_data2, 'r', marker='s', fillstyle='none', label='Sea Ice FSD Exponent: Data') #Same Size Data
# plt.plot(t_sim, slope_gsd_sim2, 'r.', label='Sea Ice FSD Exponent: Simulation RMSE = '+str(round(Slope_RMSE,4)) )
# plt.legend()
# plt.title('Evolution of Sea Ice FSD Exponent')
# plt.savefig(figDir9 + 'SeaIce_FSD_Slope_Compare_Expf.png', format='png')
# plt.close()

#Initial FSD
print("Initial FSD:")
print("Diameters: ", sizeV)
print("Number of floes", countV)

exit(1)

#Save Conc. and Normalized Area
timeC = []
nArea = []
Conc = []
ConcB = []
ConcF = []
surfArea = []
aveThick = []
globalOceanTemp = []
aveThickC = []
aveThickF = []
for i in range(len(concV)):
    timeC.append(concV[i][0])
    nArea.append(concV[i][1])
    Conc.append(concV[i][2]*100)
    ConcB.append(concV[i][3]*100)
    ConcF.append(concV[i][4]*100)
    surfArea.append(concV[i][5])
    aveThick.append(concV[i][6])
    globalOceanTemp.append(concV[i][7])
    aveThickC.append(concV[i][8])
    aveThickF.append(concV[i][9])
    
figG2 = plt.figure(figsize=(10,10))
plt.plot(timeC,nArea,'ro-')

plt.title('Evolution Normalized Area')

#plt.xlim(MinR,MaxR) #Must be log limit so avoid zeros
#plt.ylim(0.0,100)  #Must be log limit so avoid zeros

plt.xlabel('Time')
plt.ylabel('Normalized Area')
plt.ylim(0.0,1.0)  #Must be log limit so avoid zeros
plt.savefig(figConcDir + 'Normalized_area.png', format='png')
plt.close()


#Concentration results
#Big y = -1.33662 * x + 33.15290
#Small y = -1.26328 * x + 37.81386
#Total y = -2.59990 * x + 70.96676 #Slope = -2.5999 % per day or 0.025999 per day or 3.0081 e-7
#aR = -2.59990
#bR = 70.96676

#Plot Surface Area over time
figGS = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Surface Area (Perimeter)')
poptSf, _ = curve_fit(objective, timeC, surfArea)
# summarize the parameter values
aSf = poptSf[0]
bSf = poptSf[1]

plt.plot(timeC,surfArea,'r.', label='Surface Area, slope: ' + str(round(aSf,3)))

x_lineSf = np.arange(timeC[0], timeC[-1], 1)
y_lineSf = objective(x_lineSf, aSf, bSf)
plt.plot(x_lineSf, y_lineSf, '--', color='red')
plt.legend()
plt.title('Evolution of Surface Area')
plt.savefig(figConcDir + 'SeaIce_SurfArea.png', format='png')
plt.close()


########################################################################################################
#Plot Average Thickness over time
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Thickness(meters)')
poptAT, _ = curve_fit(objective, timeC, aveThick)
# summarize the parameter values
aAT = poptAT[0]
bAT = poptAT[1]

plt.plot(timeC,aveThick,'r.', label='Average Thickness, slope: ' + str(round(aAT,3)))

x_lineAT = np.arange(timeC[0], timeC[-1], 1)
y_lineAT = objective(x_lineAT, aAT, bAT)
plt.plot(x_lineAT, y_lineAT, '--', color='red')
plt.legend()
plt.title('Evolution of Average Thickness')
plt.savefig(figConcDir + 'SeaIce_AveThick.png', format='png')
plt.close()


#Plot Sim Thick vs. Field Data Thick
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Thickness(meters)')
dataThick = np.zeros(len(timeC))
for i in range(len(timeC)):
    dataThick[i] = dataThickV[i][1]
plt.plot(timeC, dataThick,'b', marker='s', label='Average Thickness: Data') #Same Size Data
plt.plot(timeC,aveThick,'k*', label='Average Thickness: Simulation')
plt.plot(timeC,aveThickC,'r*', label='Average Thickness Coarse: Simulation')
plt.plot(timeC,aveThickF,'g*', label='Average Thickness Fine: Simulation')
plt.legend()
plt.title('Evolution of Average Thickness')
plt.savefig(figDir8 + 'SeaIce_AveThick_Compare.png', format='png')
plt.close()

#Plot Just Thickness Data Series
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Thickness(meters)')
plt.plot(timeC, dataThick,'b', marker='s', label='Average Thickness: Data') #Same Size Data
plt.title('Evolution of Average Thickness')
plt.savefig(figDir8 + 'SeaIce_AveThick_Data.png', format='png')
plt.close()

#Plot simpler Image
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Thickness(meters)')
dataThick = np.zeros(len(timeC))
for i in range(len(timeC)):
    dataThick[i] = dataThickV[i][1]
plt.plot(timeC, dataThick,'b', marker='s', label='Average Thickness: Data') #Same Size Data
plt.plot(timeC,aveThick,'k*', label='Average Thickness: Simulation')
plt.legend()
plt.title('Evolution of Average Thickness')
plt.savefig(figDir8 + 'SeaIce_AveThick_CompareTot.png', format='png')
plt.close()

#############################################################################################################


#Plot Global Grid Temperature
figGTemp = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Temperature (Celsius)')
poptOT, _ = curve_fit(objective, timeC, globalOceanTemp)
# summarize the parameter values
aOT = poptOT[0]
bOT = poptOT[1]
print('globalOceanTemp: ', globalOceanTemp)

plt.plot(timeC,globalOceanTemp,'r.', label='Average Ocean Temp., slope: ' + str(round(aOT,3)))

x_lineOT = np.arange(timeC[0], timeC[-1], 1)
y_lineOT = objective(x_lineOT, aOT, bOT)
plt.plot(x_lineOT, y_lineOT, '--', color='red')
plt.legend()
plt.title('Evolution of Average Ocean Temp.')
plt.savefig(figConcDir + 'SeaIce_AveGlobalOceanTemp.png', format='png')
plt.close()

#Get data vector for RMSE and plot
dataTemp = np.zeros(len(timeC))
for i in range(len(timeC)):
    dataTemp[i] = dataTempV[i][1]

#Find temp RMSE
Temp_RMSE = 0
for i in range(len(dataTemp)):
    tstepD = timeC[i]
    for j in range(len(globalOceanTemp)):
        if (tstepD == timeC[j]):
            indexstepN = j
            break
    Temp_RMSE += (dataTemp[i]-globalOceanTemp[indexstepN])**2
Temp_RMSE /= len(dataTemp)
Temp_RMSE = Temp_RMSE**0.5

#Plot Sim Temp vs. Field Data Temp
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Temperature (Celsius)')
plt.plot(timeC, dataTemp,'b', marker='s', label='Average Ocean Temp.: Data') #Same Size Data
plt.plot(timeC,globalOceanTemp,'r*', label='Average Ocean Temp.: Simulation RMSE =' +str(round(Temp_RMSE,4)))
plt.legend()
plt.title('Evolution of Average Ocean Temperature')
plt.savefig(figDir8 + 'SeaIce_AveOceanTemp_Compare.png', format='png')
plt.close()

#Plot Just Thickness Data Series
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Temperature (Celsius)')
plt.plot(timeC, dataTemp,'b', marker='s', label='Average Ocean Temp.: Data') #Same Size Data
plt.title('Evolution of Average Ocean Temperature')
plt.savefig(figDir8 + 'SeaIce_AveTemp_Data.png', format='png')
plt.close()


figG3 = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Sea Ice Concentration')
#plt.ylim(0.0,1.0)  #Must be log limit so avoid zeros


poptT, _ = curve_fit(objective, timeC, Conc)
poptTB, _ = curve_fit(objective, timeC, ConcB)
poptTF, _ = curve_fit(objective, timeC, ConcF)
# summarize the parameter values
aT = poptT[0]
bT = poptT[1]
aB = poptTB[0]
bB = poptTB[1]
aF = poptTF[0]
bF = poptTF[1]
plt.plot(timeC,Conc,'k.', label='Total floe Concentration, slope: ' + str(round(aT,3)))
plt.plot(timeC,ConcB,'r.', label='Big floe Concentration, slope: ' + str(round(aB,3)))
plt.plot(timeC,ConcF,'g.', label='Small floe Concentration, slope: ' + str(round(aF,3)))
print('Concentration Linear Fit Total y = %.5f * x + %.5f' % (aT, bT))
print('Concentration Linear Fit Total BIG y = %.5f * x + %.5f' % (aB, bB))
print('Concentration Linear Fit Total FINE y = %.5f * x + %.5f' % (aF, bF))
x_line = np.arange(timeC[0], timeC[-1], 1)
y_line = objective(x_line, aT, bT)
y_lineB = objective(x_line, aB, bB)
y_lineF = objective(x_line, aF, bF)
plt.plot(x_line, y_line, '--', color='black')
plt.plot(x_line, y_lineB, '--', color='red')
plt.plot(x_line, y_lineF, '--', color='green')
plt.legend()

#Image based line (needs to be rescaled)
#x_line = np.arange(0, 30, 1)
#y_line = objective(x_line, aR, bR)
#plt.plot(x_line, y_line, '--', color='green')

plt.title('Evolution of Concentration with Slope Full = ' + str(round(aT,5)) + ' Big: ' + str(round(aB,5)) + ' Fine: ' + str(round(aF,5)))
plt.savefig(figConcDir + 'SeaIce_Concentration.png', format='png')
plt.close()

#Plot Sim Concentration vs. Field Data Concentration
figGAT = plt.figure(figsize=(10,10))
plt.xlabel('Time')
plt.ylabel('Concentration (%)')
timeCData = np.zeros(len(dataConcV))
dataConcC  = np.zeros(len(dataConcV))
dataConcF  = np.zeros(len(dataConcV))
dataConcT  = np.zeros(len(dataConcV))

for i in range(len(timeCData)):
    timeCData[i]    = dataConcV[i][0]
    dataConcC[i]    = dataConcV[i][5]
    dataConcF[i]    = dataConcV[i][6]
    dataConcT[i]    = dataConcV[i][4]

#Calculate RSME (assume all data has same length)
C_RMSE = 0
F_RMSE = 0
T_RMSE = 0
for i in range(len(dataConcC)):
    tstepD = timeCData[i]
    for j in range(len(ConcB)):
        if (tstepD == timeC[j]):
            indexstepN = j
            break
    C_RMSE += (dataConcC[i]-ConcB[indexstepN])**2
    F_RMSE += (dataConcF[i]-ConcF[indexstepN])**2
    T_RMSE += (dataConcT[i]-Conc[indexstepN])**2
C_RMSE /= len(dataConcC)
F_RMSE /= len(dataConcC)
T_RMSE /= len(dataConcC)
C_RMSE = C_RMSE**0.5
F_RMSE = F_RMSE**0.5
T_RMSE = T_RMSE**0.5


#Add parameter data text box
textstr = '\n'.join((
    r'$\mathrm{Total Steps}=%d$' % (nT, ),
    r'$\mathrm{nTemp}=%d$' % (nTemp, ),
    r'$\mathrm{nBreak}=%d$' % (nBreak, ),
    r'$\mathrm{qvert}=%.2f$' % (qvert, ),
    r'$\mathrm{Khor}=%.2f$' % (Khor, ),
    r'$\mathrm{Qatm}=%.2f$' % (Qatm, ),
    r'$\mathrm{Aatm}=%.2f$' % (Aatm, ),
    r'$\mathrm{Batm}=%.2f$' % (Batm, ),
    r'$\mathrm{alpha_ice}=%.2f$' % (alpha_ice, ),
    r'$\mathrm{alpha_ocean}=%.2f$' % (alpha_ocean, ),
    r'$\mathrm{MaxD}=%d$' % (limitMaxDiam, ),
    r'$\mathrm{MinD}=%.4f$' % (limitMinDiam, )
    ))


#Save parameters to our results folder
shutil.copy(testParamsFile, figDir8)


#Save parameter and Errors for Conc. and Temperature
vec_export_text = []
vec_export_text.append([ nTemp, nBreak, Qatm, qvert, C_RMSE, F_RMSE, T_RMSE, Temp_RMSE ])
print(vec_export_text)
np.savetxt(figDir8+'Param_Errors_'+str(caseNo)+'.dat', vec_export_text, fmt='%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f')


# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# place a text box in upper left in axes coords
plt.text(38, 67, textstr, fontsize=10,
        verticalalignment='top', bbox=props)

#Could add 1 for meanD, 2 for maxD, 3 for minD
plt.plot(timeCData, dataConcC, 'r', marker='s', fillstyle='none', label='Sea Ice Concentration Coarse: Data') #Same Size Data
plt.plot(timeC,ConcB, 'r.', label='Sea Ice Concentration Coarse: Simulation RMSE = '+str(round(C_RMSE,4)) )
plt.plot(timeCData, dataConcF, 'g', marker='s', fillstyle='none', label='Sea Ice Concentration Fine: Data') #Same Size Data
plt.plot(timeC,ConcF, 'g.', label='Sea Ice Concentration Fine: Simulation RMSE = '+str(round(F_RMSE,4)))
plt.plot(timeCData, dataConcT, 'k', marker='s', fillstyle='none', label='Sea Ice Concentration Total: Data') #Same Size Data
plt.plot(timeC,Conc, 'k.', label='Sea Ice Concentration Total: Simulation RMSE = '+str(round(T_RMSE,4)))
plt.legend()
plt.title('Evolution of Sea Ice Concentration')
plt.savefig(figDir8 + 'SeaIce_Concentration_Compare.png', format='png')
plt.close()

#Plot Only Coarse Concentration Performance
figG3 = plt.figure(figsize=(10,10))
plt.plot(timeC,ConcB, 'r.', label='Sea Ice Concentration Coarse: Simulation RMSE = '+str(round(C_RMSE,4)) )
plt.plot(timeCData, dataConcC, 'r', marker='s', fillstyle='none', label='Sea Ice Concentration Coarse: Data') #Same Size Data
plt.legend()
plt.xlabel('Time (days)')
plt.ylabel('Sea Ice Concentration (%)')
plt.title('Evolution of Sea Ice Concentration: Coarse')
plt.savefig(figDir8 + 'SeaIce_Concentration_CompareC.png', format='png')
plt.close()

#Plot Only Fines Concentration Performance
figG3 = plt.figure(figsize=(10,10))
plt.plot(timeC,ConcF, 'g.', label='Sea Ice Concentration Fine: Simulation RMSE = '+str(round(F_RMSE,4)))
plt.plot(timeCData, dataConcF, 'g', marker='s', fillstyle='none', label='Sea Ice Concentration Fine: Data') #Same Size Data
plt.legend()
plt.xlabel('Time (days)')
plt.ylabel('Sea Ice Concentration (%)')
plt.title('Evolution of Sea Ice Concentration: Fines')
plt.savefig(figDir8 + 'SeaIce_Concentration_CompareF.png', format='png')
plt.close()

#Plot Diameters (Max, Mean, Min, nFloes)
nFloes = numberg
Diam_Max = []
Diam_Mean = []
Diam_Min = []
for i in range(len(DiamV)):
    Diam_Max.append(DiamV[i][0])
    Diam_Mean.append(DiamV[i][1])
    Diam_Min.append(DiamV[i][2])

# curve fit for all results
poptMean, _ = curve_fit(objective, timeC, Diam_Mean)
# summarize the parameter values
aMean = poptMean[0]
bMean = poptMean[1]
print('MeanD y = %.5f * x + %.5f' % (aMean, bMean))

poptMax, _ = curve_fit(objective, timeC, Diam_Max)
# summarize the parameter values
aMax = poptMax[0]
bMax = poptMax[1]
print('MaxD y = %.5f * x + %.5f' % (aMax, bMax))

poptMin, _ = curve_fit(objective, timeC, Diam_Min)
# summarize the parameter values
aMin = poptMin[0]
bMin = poptMin[1]
print('MinD y = %.5f * x + %.5f' % (aMin, bMin))

print(len(timeC))
print((timeC))
print(len(nFloes))
print((nFloes))
poptN, _ = curve_fit(objective, timeC, nFloes)
# summarize the parameter values
aN = poptN[0]
bN = poptN[1]
print('Nfloes y = %.5f * x + %.5f' % (aN, bN))

## Plotting meanD and maxD-minD
fig, ax1 = plt.subplots(figsize=(10,10))
color = 'tab:red'
ax1.set_xlabel('Time')
ax1.set_ylabel('Diameter', color=color)
ax1.plot(timeC, Diam_Mean, 'r*', label='Mean Diameter, slope: ' + str(round(aMean,3)))
ax1.plot(timeC, Diam_Max, 'g*', label='MaxD, slope: ' + str(round(aMax,3)))
ax1.plot(timeC, Diam_Min, 'ko', label='MinD, slope: ' + str(round(aMin,3)))
ax1.tick_params(axis='y', labelcolor=color)

# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(timeC[0], timeC[-1], 1)
# calculate the output for the range
y_line = objective(x_line, aMean, bMean)
# create a line plot for the mapping function
ax1.plot(x_line, y_line, '--', color='red')


# define a sequence of inputs between the smallest and largest known inputs
x_line2 = np.arange(timeC[0], timeC[-1], 1)
# calculate the output for the range
y_line2 = objective(x_line2, aMax, bMax)
# create a line plot for the mapping function
ax1.plot(x_line2, y_line2, '--', color='green')


# define a sequence of inputs between the smallest and largest known inputs
x_line3 = np.arange(timeC[0], timeC[-1], 1)
# calculate the output for the range
y_line3 = objective(x_line3, aMin, bMin)
# create a line plot for the mapping function
ax1.plot(x_line3, y_line3, '--', color='black')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('# of Floes', color=color)  # we already handled the x-label with ax1
ax2.plot(timeC, nFloes, 'bo', label='No. of Floes, slope: ' + str(round(aN,3)))
ax2.tick_params(axis='y', labelcolor=color)

# define a sequence of inputs between the smallest and largest known inputs
x_line4 = np.arange(timeC[0], timeC[-1], 1)
# calculate the output for the range
y_line4 = objective(x_line4, aN, bN)
# create a line plot for the mapping function
ax2.plot(x_line4, y_line4, '--', color='blue')


#fig.tight_layout()  # otherwise the right y-label is slightly clipped

#Revert X-Axis
#axG = plt.gca()
#axG.set_xlim(axG.get_xlim()[::-1])
ax1.legend(loc = 'upper right')
ax2.legend(loc = 'lower right')

plt.title('Evolution of Diameters with Slope Full')
plt.savefig(figConcDir + 'SeaIce_Diameters.png', format='png')
plt.close()


#Compare Data and Simulation Diameters

#Get Field Data

#Plot Sim Diameter vs. Field Data Diameters
timeDData = np.zeros(len(dataConcV))
dataMeanD  = np.zeros(len(dataConcV))
dataMaxD = np.zeros(len(dataConcV))
dataMinD  = np.zeros(len(dataConcV))
dataN = np.zeros(len(dataConcV))

for i in range(len(timeDData)):
    timeDData[i]    = dataConcV[i][0]
    dataMeanD[i]    = dataConcV[i][1]
    dataMaxD[i]    = dataConcV[i][2]
    dataMinD[i]    = dataConcV[i][3]
    dataN[i]    = dataConcV[i][7]

#Calculate RSME (assume all data has same length)
RMSE_Mean = 0
RMSE_Max = 0
RMSE_Min = 0
RMSE_N = 0
for i in range(len(dataMeanD)):
    tstepD = timeDData[i]
    for j in range(len(Diam_Max)):
        if (tstepD == timeC[j]):
            indexstepN = j
            break
    RMSE_Mean += (dataMeanD[i] - Diam_Mean[indexstepN])**2
    RMSE_Max += (dataMaxD[i] - Diam_Max[indexstepN])**2
    RMSE_Min += (dataMinD[i] - Diam_Min[indexstepN])**2
    RMSE_N   += (dataN[i] - nFloes[indexstepN])**2
RMSE_Mean /= len(dataMeanD)
RMSE_Max /= len(dataMeanD)
RMSE_Min /= len(dataMeanD)
RMSE_N /= len(dataMeanD)
RMSE_Mean = RMSE_Mean**0.5
RMSE_Max = RMSE_Max**0.5
RMSE_Min = RMSE_Min**0.5
RMSE_N = RMSE_N**0.5


## Plotting meanD and maxD-minD
fig, ax1 = plt.subplots(figsize=(10,10))
color = 'tab:red'
ax1.set_xlabel('Time')
ax1.set_ylabel('Diameter', color=color)
ax1.plot(timeDData, dataMeanD, 'g', marker='s', fillstyle='none', label='MeanD = Data') #Same Size Data
ax1.plot(timeC, Diam_Mean, 'g*', label='MeanD = Simulation, RMSE: ' + str(round(RMSE_Mean,3)))
ax1.plot(timeDData, dataMaxD, 'r', marker='s', fillstyle='none', label='MaxD = Data') #Same Size Data
ax1.plot(timeC, Diam_Max, 'r*', label='MaxD = Simulation, RMSE: ' + str(round(RMSE_Max,3)))
ax1.plot(timeDData, dataMinD, 'k', marker='s', fillstyle='none', label='MinD = Data') #Same Size Data
ax1.plot(timeC, Diam_Min, 'k*', label='MinD = Simulation, RMSE:' + str(round(RMSE_Min,3)))
ax1.tick_params(axis='y', labelcolor=color)

color = 'tab:blue'
ax2.set_ylabel('# of Floes', color=color)  # we already handled the x-label with ax1
ax2.plot(timeDData, dataN, 'b', marker='s', fillstyle='none', label='No. of Floes = Data') #Same Size Data
ax2.plot(timeC, nFloes, 'b*', label='No. of Floes = Simulation, RMSE: ' + str(round(RMSE_N,3)))
ax2.tick_params(axis='y', labelcolor=color)

ax1.legend(loc = 'upper right')
ax2.legend(loc = 'lower right')

plt.title('Evolution of Diameters Simulation vrs. MODIS Data')
plt.savefig(figDir8 + 'SeaIce_Diameters_Compare.png', format='png')
plt.close()


## Plotting Nfloes Compare
fig, ax1 = plt.subplots(figsize=(10,10))
ax1.set_xlabel('Time')
ax1.set_ylabel('# of Floes')
ax1.plot(timeDData, dataN, 'b', marker='s', fillstyle='none', label='No. of Floes = Data') #Same Size Data
ax1.plot(timeC, nFloes, 'b*', label='No. of Floes = Simulation, RMSE: ' + str(round(RMSE_N,3)))
ax1.legend(loc = 'upper right')

plt.title('Evolution of Floe Number Simulation vrs. MODIS Data')
plt.savefig(figDir8 + 'SeaIce_FloeN_Compare.png', format='png')
plt.close()


#Plot variation of concentration for each size bin (constant bins for all times)
for j in range(len(gsd0Max)):
    
    figC2 = plt.figure(figsize=(10,10))
    plt.plot(timeC,varC[j],'r.-')  #varC is a vector of size time from a matrix len(gsD[0])*len(timeC)

    if ( j < len(gsd0Max)-1 ):
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: '  + str(round(gsd0Max[j+1][0],2)))
    else:
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: 0')
        
    #plt.xlim(MinR,MaxR) #Must be log limit so avoid zeros
    #plt.ylim(0.0,100)  #Must be log limit so avoid zeros

    plt.xlabel('Time')
    plt.ylabel('Concentration %')
    plt.ylim(0.0,20.0)  #Must be log limit so avoid zeros
    plt.savefig(figConcDir + 'Bin_Conc_'+str(j)+'.png', format='png')
    plt.close()


#Plot variation of Melt for each size bin (constant bins for all times)
for j in range(len(gsd0Max)):
    
    figC2 = plt.figure(figsize=(10,10))
    plt.plot(timeC,gsdMelt[j],'r.-')  #varC is a vector of size time from a matrix len(gsD[0])*len(timeC)

    if ( j < len(gsd0Max)-1 ):
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: '  + str(round(gsd0Max[j+1][0],2)))
    else:
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: 0')
        
    #plt.xlim(MinR,MaxR) #Must be log limit so avoid zeros
    #plt.ylim(0.0,100)  #Must be log limit so avoid zeros

    plt.xlabel('Time')
    plt.ylabel('Accumulated Mass Change')
    #plt.ylim(0.0,20.0)  #Must be log limit so avoid zeros
    plt.savefig(figConcDir + 'Bin_Melt_'+str(j)+'.png', format='png')
    plt.close()
    
#Plot variation of Break for each size bin (constant bins for all times)
for j in range(len(gsd0Max)):
    
    figC2 = plt.figure(figsize=(10,10))
    plt.plot(timeC,gsdBreak[j],'r.-')  #varC is a vector of size time from a matrix len(gsD[0])*len(timeC)

    if ( j < len(gsd0Max)-1 ):
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: '  + str(round(gsd0Max[j+1][0],2)))
    else:
        plt.title('Evolution of concentration from: ' + str(round(gsd0Max[j][0],2)) + ' to: 0')
        
    #plt.xlim(MinR,MaxR) #Must be log limit so avoid zeros
    #plt.ylim(0.0,100)  #Must be log limit so avoid zeros

    plt.xlabel('Time')
    plt.ylabel('Accumulated Mass Change')
    #plt.ylim(0.0,20.0)  #Must be log limit so avoid zeros
    plt.savefig(figConcDir + 'Bin_Break_'+str(j)+'.png', format='png')
    plt.close()


#GSD Error Evolution for Adjusting Simulation
fileImgEnd = "_GSD_Data.dat"
fileSimDir = figDir4 + "GrainSizeDist_iter"
fileSimEnd = "_0.dat"

if year == 2020:
    imgNo = [0, 1, 2, 9, 10, 14, 15, 16, 17, 22, 24] #2020
elif year == 2018:
    imgNo = [0, 5, 8, 9, 10, 14, 15, 17, 18, 22, 32, 35, 38, 42, 43, 44, 48] #2018
FileImg = []
FileData = []
DayNumber = imgNo

for i in range(len(imgNo)):
    FileImg.append(fileImgDir + str(imgNo[i]) + fileImgEnd)
    if (imgNo[i] - imgNo[0])/(imgNo[-1] - imgNo[0]) == 1:
        proportion = int(  ( (imgNo[i] - imgNo[0])/(imgNo[-1] - imgNo[0])  * nSteps) - 1 )
    else:
        proportion = int(  ( (imgNo[i] - imgNo[0])/(imgNo[-1] - imgNo[0])  * nSteps) )
    FileData.append(fileSimDir + str(proportion) + fileSimEnd)

print(FileImg)
print(FileData)
x = []
Mean_error_gsdP = np.zeros(len(imgNo))
Mean_error_gsdConc = np.zeros(len(imgNo))
Mean_Perror_gsdP = np.zeros(len(imgNo))
Mean_Perror_gsdConc = np.zeros(len(imgNo))
for i in range(len(imgNo)):
    x.append(imgNo[i])
    [Mean_error_gsdP[i], Mean_error_gsdConc[i], Mean_Perror_gsdP[i], Mean_Perror_gsdConc[i]] = GSD_Comparison_Error(FileImg[i], FileData[i], DayNumber[i], outErrorDir)

#Plot Error Evolution over time
figEr = plt.figure(figsize=(10,10))
plt.plot(x, Mean_error_gsdP ,'r*-')
plt.plot(x, Mean_error_gsdConc ,'b.-')
plt.title('Evolution of Error for Pass and Concentration')
plt.xlabel('Time')
plt.ylabel('Error')
plt.savefig(outErrorDir + 'PassandConc_Error.png', format='png')
plt.close()

figEr = plt.figure(figsize=(10,10))
plt.plot(x, Mean_error_gsdConc ,'b.-')
plt.title('Evolution of Error for Pass and Concentration')
plt.xlabel('Time')
plt.ylabel('Error')
plt.savefig(outErrorDir + 'OnlyConc_Error.png', format='png')
plt.close()


figEr2 = plt.figure(figsize=(10,10))
plt.plot(x, Mean_Perror_gsdP ,'r*-')
plt.plot(x, Mean_Perror_gsdConc ,'b.-')
plt.title('Evolution of Error % for Pass and Concentration')
plt.xlabel('Time')
plt.ylabel('% Error')
plt.savefig(outErrorDir + 'PassandConcPerc_Error.png', format='png')
plt.close()

#For paper Output comparison only (TURN FOR TEXT RESYLT COMPARISON)
# #Save Info Snapshot
# sample_days = [0, 8, 15, 22, 32]
# data_sample_days = [0, 2, 6, 9, 10]  #To avoid index failure  [0, 5, 8, 9, 10, 14, 15, 17, 18, 22, 32, 35, 38, 42, 43, 44, 48] 
# vec_export_text = []
# for i in range(len(sample_days)):
#     vec_export_text.append([ sample_days[i], dataN[data_sample_days[i]], nFloes[sample_days[i]], dataConcC[data_sample_days[i]], ConcB[sample_days[i]], dataConcF[data_sample_days[i]], ConcF[sample_days[i]], dataConcT[data_sample_days[i]], Conc[sample_days[i]], dataThick[data_sample_days[i]], aveThick[sample_days[i]], dataMeanD[data_sample_days[i]], Diam_Mean[sample_days[i]], dataTemp[data_sample_days[i]], globalOceanTemp[sample_days[i]] ])  ##  #Day  # no_floes_data  #no_floes_sim  #ConcC_data #ConcC_sim  #ConcF_data #ConcF_sim #ConcT_data #ConcT_sim #Avethick_data #Avethick_sim #Ave Diam Data # Ave Diam Sim  #aveTempData #AveTempSim
# np.savetxt(figDir8+'General_Data_'+str(caseNo)+'.dat', vec_export_text, fmt='%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f')

# vec_export_text2 = []
# #For 7 fields
# for i in range(len(sample_days)):
#     vec_export_text2.append([ sample_days[i] , dataN[data_sample_days[i]], nFloes[sample_days[i]]  ])
#     vec_export_text2.append([ sample_days[i] ,  dataConcC[data_sample_days[i]], ConcB[sample_days[i]] ])
#     vec_export_text2.append([ sample_days[i] , dataConcF[data_sample_days[i]], ConcF[sample_days[i]]  ])   
#     vec_export_text2.append([ sample_days[i] , dataConcT[data_sample_days[i]], Conc[sample_days[i]]  ])   
#     vec_export_text2.append([ sample_days[i] , dataThick[data_sample_days[i]], aveThick[sample_days[i]]   ])  
#     vec_export_text2.append([ sample_days[i] , dataMeanD[data_sample_days[i]], Diam_Mean[sample_days[i]]    ])  
#     vec_export_text2.append([ sample_days[i] , dataTemp[data_sample_days[i]], globalOceanTemp[sample_days[i]]    ])  
# np.savetxt(figDir8+'General_Data2_'+str(caseNo)+'.dat', vec_export_text2, fmt='%.5f %.5f %.5f')


if video_gen:
    ##VIDEOS
    #INPUT / OUTPUT LOCATIONS
    datadirin = mainVfolder
    datadirout = mainVfolder
    
    #DIRECTORY
    vid_add = datadirin
    #Number of frames
    frames = len(glob.glob1(vid_add,"*.png"))
    
    img_array = []
    
    for i in range(0,frames):
        vid_num = 'step_'+ str(i)
        filename = vid_add + vid_num + '.png'
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)
    
    
    out = cv2.VideoWriter(datadirout + 'VideoPoints.mp4',cv2.VideoWriter_fourcc(*'mp4v'), 5, size) #10 frame default
    #out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
    
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()
    
    
    ##INPUT / OUTPUT LOCATIONS
    datadirin3 = mainVfolder+"GSD/"
    datadirout3 = mainVfolder+"GSD/"
    
    #DIRECTORY
    vid_add = datadirin3
    #Number of frames
    frames = len(glob.glob1(vid_add,"*.png"))
    
    img_array = []
    
    for i in range(0,frames):
        vid_num = 'step_'+ str(i)
        filename = vid_add + vid_num + '.png'
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)
    
    
    out = cv2.VideoWriter(datadirout3 + 'VideoGSD.mp4',cv2.VideoWriter_fourcc(*'mp4v'), 10, size)
    #out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
    
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()
    
    
    ###INPUT / OUTPUT LS
    if LS_plotting == 1:
        datadirin4 = mainVfolder+"LevelSet0/"
        datadirout4 = mainVfolder+"LevelSet0/"
        #
        ##DIRECTORY
        vid_add = datadirin4
        #Number of frames
        frames = len(glob.glob1(vid_add,"*.png"))
    
        img_array = []
    
        for i in range(0,frames):
            vid_num = 'step_'+ str(i)
            filename = vid_add + vid_num + '.png'
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width,height)
            img_array.append(img)
    
    
        out = cv2.VideoWriter(datadirout4 + 'VideoLS.mp4',cv2.VideoWriter_fourcc(*'mp4v'), 10, size)
        #out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
    
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()












#GLOB ALTERNATIVE METHOD

# for filename in glob.glob(vid_add + '*.png'):
#     img = cv2.imread(filename)
#     height, width, layers = img.shape
#     size = (width,height)
#     img_array.append(img)
#     print(filename)


# To run code on VM make sure to enter this into commandline to work on CV environment:

#export WORKON_HOME=$HOME/.virtualenvs
#export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3
#source /usr/local/bin/virtualenvwrapper.sh
#mkvirtualenv cv -p python3
#workon cv

#export WORKON_HOME=$HOME/.virtualenvs
#export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3
#source /usr/local/bin/virtualenvwrapper.sh
#mkvirtualenv cv -p python3
#workon cv
#pip install opencv-python-headless
#pip3 install scikit-image

