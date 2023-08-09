import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.optimize import curve_fit
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
import sys
import shutil

# objective function
def objective(x, a, b):
    return a * x + b

#Accum. Sum
def acc_sum(vec):
    out_vec = np.zeros(len(vec))
    for i in range(len(vec)):
        if i < len(vec) - 1:
            out_vec[i] = np.sum(vec[0:i+1])
        else:
            out_vec[i] = np.sum(vec[0:len(vec)])
    return out_vec

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

#GSD Slope Linear
def gsd_slope_linear(objective, x, y):
    # curve fit for all results
    poptMean, _ = curve_fit(objective, x, y)
    # summarize the parameter values
    aMean = poptMean[0]
    bMean = poptMean[1]
    return aMean

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

##**************************************************************************##
##**************************************************************************##

#Based on Code to Populate an area with random grains but for Arbitrary GSD
#Start using some real info and fill in the blanks with additional grains
# Output contact and updated position file

rho_ice = 910e9

#Directories

caseNo = sys.argv[1]
year = int(sys.argv[2])

if year == 2018:
    inputNew = "./Input/"
    inputDir = "./Input/SeaIceWGJ2018/"
    figDir = "./Visualization/Initial_GSD_2018/"
    grainDir2  = "./Input/grainsIceWGJ2018/"
    ingDir1 = "./Input/grainsIceWGJ2018/"
    mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/" #2018
    #mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #2018 #PROB. CHANGE
elif year == 2020:
    inputNew = "./Input/"
    inputDir = "./Input/SeaIceWGJ2020/"
    figDir = "./Visualization/Initial_GSD_2020/"
    grainDir2  = "./Input/grainsIceWGJ/"
    ingDir1 = "./Input/grainsIceWGJ/"
    mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" #2020

if os.path.exists(figDir) == False:
    os.mkdir(figDir)

#Domain_Area in km2
if year == 2020:
    Domain_Area = 259.0 * 339.66 #For 2020
elif year == 2018:
    Domain_Area = 0.076011266724044 * 1474000; #For 2018 Example, see imageprocess2_label_v2018.py (1340 x, 11000 y in pixels)  #Is 112040.6071512409

#All grains (Area and Thickness)
TSDC_fileOr = mainOutputfolder + "Vertical_Thickness/" + "TSD_Coarse_step_0.dat"
TSDCV = np.loadtxt(TSDC_fileOr)

#Bins Concentration and Floe Number
GSDFileOr = mainOutputfolder + "GSD/" + "GrainSizeDist_iter0_0.dat"
gsd0F = np.loadtxt(GSDFileOr)

#print(TSDCV)
#print(GSD_OrV)

#Analyze Diameters
Dmax = 0
Dmin = 100000000000
total_mass = 0
ave_thick = 0
total_area = 0
for i in range(len(TSDCV)):
    total_mass += rho_ice * TSDCV[i][0] * TSDCV[i][1]
    ave_thick += TSDCV[i][1]
    total_area += TSDCV[i][0]
    if MeanCaliperD(TSDCV[i][0]) > Dmax:
        Dmax = MeanCaliperD(TSDCV[i][0])
    if MeanCaliperD(TSDCV[i][0]) < Dmin:
        Dmin = MeanCaliperD(TSDCV[i][0])
  
ave_thick /= len(TSDCV)    
Dmean = 0.5*(Dmax - Dmin)

print("Dmax/Dmin/Dmean = ", Dmax, Dmin, Dmean)

#For more resolution
# #Vector for Slope from all data
# diamV0 = []
# for i in range(len(TSDCV)):
#     diamV0.append(MeanCaliperD(TSDCV[i][0]))
# #Order list
# diamV0.sort(reverse=True)
# print(diamV0)

#Analyze Slope (GSD)
log_gsdDF = []
log_gsdNFF = []
gsdNF_or = []
coarse_conc = 0
fine_conc = 0
for j in range(len(gsd0F)):
    if j < len(gsd0F) - 1:
        coarse_conc += gsd0F[j][3]
    else:
        fine_conc += gsd0F[j][3]
    if gsd0F[j][0] <= 0:
        log_gsdDF.append(0.0)
    else:
        log_gsdDF.append(math.log(gsd0F[j][0]))
    if gsd0F[j][2] <= 0:
        gsdNF_or.append(0.0)
    else:
        gsdNF_or.append(abs(gsd0F[j][2]))
#Convert to cumulative sum
gsdNF_orp = gsdNF_or
print("Original Floe Count: ", gsdNF_or)
gsdNF_or = acc_sum(gsdNF_orp)
print("Original Cum Floe Count: ", gsdNF_or)

for idx in range(len(gsdNF_or)):
    if gsdNF_or[idx] <=0:
        log_gsdNFF.append(0.0)
    else:
        log_gsdNFF.append(math.log(gsdNF_or[idx]))

#Filter zero results if occurring
log_gsdDF0 = log_gsdDF
log_gsdNFF0 = log_gsdNFF
[log_gsdDF, log_gsdNFF] = zeros_filter_slope(log_gsdDF0, log_gsdNFF0)
slope_simNO  = abs(gsd_slope_linear(objective, log_gsdDF, log_gsdNFF))
D_max_slope = log_gsdDF[0]
D_min_slope = log_gsdDF[-1]
N_max_slope = log_gsdNFF[0]
N_min_slope = log_gsdNFF[-1]
D_mean_slope = 0.5*(math.exp(D_max_slope) - math.exp(D_min_slope))
Brute_Slope = (N_min_slope-N_max_slope) / (D_max_slope - D_min_slope)

print("Count Slope: ", slope_simNO)
print("Total Mass: ", total_mass)
print("Total Area: ", total_area)
print("Conc Area: ", total_area*100/Domain_Area)
print("Average Thickness: ", ave_thick)
print("Coarse Conc.: ", coarse_conc)
print("Fine Conc.: ", fine_conc)
print("LOG Dmax_slope: ", D_max_slope)
print("LOG Dmin_slope: ", D_min_slope)
print("LOG Nmax_slope: ", N_max_slope)
print("LOG Nmin_slope: ", N_min_slope)
print("Dmax_slope: ", math.exp(D_max_slope))
print("Dmin_slope: ", math.exp(D_min_slope))
print("Nmax_slope: ", math.exp(N_max_slope))
print("Nmin_slope: ", math.exp(N_min_slope))
print("Dmean_slope: ", D_mean_slope)
print("Brute Slope: ", Brute_Slope)

#exit(1)

#Specified Desire Outcome (FIXED DMean, Variable Slope or Viceversa)
var_slope = True
var_dmean = False
Dmean_obj = Dmean #Modify      #Dmean or D_mean_slope
slope_factor = 1.15  #1033 MoreUniform1 = 1.1 
larger = True #Larger alpha
slope_obj = slope_factor *  Brute_Slope   #OR slope_simNO #Modify #Brute_Slope ?
mass_obj = total_mass
Dmax_obj = Dmax
Dmin_obj = Dmin
thick_obj = ave_thick
dev_thick = 0.1
Area_obj = total_area
coarse_conc_obj = coarse_conc
fine_conc_obj = fine_conc
nbins = 12 #Includes Dmax and Dmin
bin_size = (Dmax_obj - Dmin_obj)/(nbins - 1)
nother_bins = int(((Dmax_obj - Dmin_obj)/bin_size) - 1)

#Consider both cases
if var_slope:
    new_log_dif = slope_obj * (D_max_slope - D_min_slope)

N_max_new = N_max_slope
N_min_new = N_max_new + new_log_dif

print("Nmax_slope New: ", math.exp(N_max_new))
print("Nmin_slope New: ", math.exp(N_min_new))
print("Test slope: ", (N_min_new - N_max_new)/(D_max_slope - D_min_slope)  )
print("Obj slope: ", slope_obj)

nMax_new = math.exp(N_max_new)
nMin_new = math.exp(N_min_new)
nGrains_new = int(round(nMin_new,0))

print("Print Bins for nGRains_new: ", nGrains_new, nother_bins)


#exit(1)

cbins = np.zeros(nbins)
cbins[0] = Dmax_obj
for i in range(nother_bins):
    cbins[i+1] = Dmax_obj - bin_size*(i+1)
cbins[-1] = Dmin_obj

print(cbins)

#Import Floe Statistics
morph_file = inputDir +"morphologies.dat"
morphV = np.loadtxt(morph_file)
nGrains = int(morphV[0])
npoints = 100 #Assuming all morphologies start with 100 points

#Generate a database of grains and size properties to then customize to objective metrics and then create input files from that list
AreaV = np.zeros(nGrains)
DiamV = np.zeros(nGrains)
pointsDat = []
total_areaG = 0
max_G_index = 0
max_area = 0
min_G_index = 0
min_area = 100000000
og_morphs = np.zeros(nGrains)
for i in range(nGrains):
    pointsFile = ingDir1 + "grainproperty" + str(int(morphV[i+1])) + ".dat"
    og_morphs[i] = morphV[i+1]
    pointsExtract = open(pointsFile,"r")
    pointList = pointsExtract.readlines()[4]
    pointList = pointList.split()
    pointVec = []
    jj = 0
    for j in range(npoints):
        floatP = [float(pointList[jj]), float(pointList[jj+1])]
        Xg = floatP[0]
        Yg = floatP[1]
        vec = [Xg,Yg]
        pointVec.append(vec)
        jj = jj + 2
    pointsDat.append(pointVec)
    AreaV[i] = PointsArea(pointVec)
    if AreaV[i] > max_area:
        max_area = AreaV[i]
        max_G_index = i
    if AreaV[i] < min_area:
        min_area = AreaV[i]
        min_G_index = i
    DiamV[i] = MeanCaliperD(AreaV[i])
    total_areaG += AreaV[i]
    #print(pointVec)
    #print(AreaV[i])
    #print(DiamV[i])
    #exit(1)


#Validation only
# print("Max Size Grain Index: ", max_G_index)
# print("Max Size Grain Diameter: ", DiamV[max_G_index])
print("Min Size Grain Index: ", min_G_index)
print("Min Size Grain Diameter: ", DiamV[min_G_index])
# print("Conc Area G: ", total_areaG*100/Domain_Area)

compliance = False #Our new distribution does not fit specifications (make sure it will)
niter_comp = 0

if larger:
    #For more uniform
    coeff_small = 0.85 #0.14 seems to good for 1.1 alpha (similar to 1.15)
else:
    #For more coarse
    coeff_small = 0.04 #0.08 works for 0.93 alpha #0.04 for 0.88
print("----Start iteration to get desired FSD----")
while compliance == False:
    #Create grain Database
    grainsArea = np.zeros(nGrains_new)
    grainsDiam = np.zeros(nGrains_new)
    grainsThick = np.zeros(nGrains_new)
    morph_new = np.zeros(nGrains_new)
    mass_new = 0
    area_new = 0
    #Smallest and Biggest grains
    grainsArea[0] = AreaV[max_G_index]
    grainsDiam[0] = DiamV[max_G_index]
    grainsThick[0] = thick_obj
    # grainsArea[1] = AreaV[min_G_index]
    # grainsDiam[1] = DiamV[min_G_index]
    # grainsThick[1] = thick_obj
    Ncum_curve = np.zeros(len(cbins))
    Ncum_curve[0] = 1
    morph_new[0] = max_G_index
    # Ncum_curve[-1] = 1
    
    
    area_tol = coeff_small * Area_obj
    
    if larger:
        reduce_small = 0.9 #0.8 - 1.0 #0.5
    else:
        reduce_small = 1.0
    
    for i in range(nGrains_new-1):
        #random grain from all list
        if area_new < Area_obj - area_tol:
            if larger:
                #Just choose random and the populate with small the remaining area
                rand_index = random.sample(range(nGrains), 1)
            else:
                #Just random but Dmean or smaller more or less, and then fill with small after reaching a tolerance
                good_size = False
                while good_size == False:
                    rand_index = random.sample(range(nGrains), 1)
                    if DiamV[rand_index] >= Dmean_obj*0.6:  #0.6 for 0.93 alpha, 
                        good_size = True
        elif area_new >= Area_obj - area_tol and area_new <= Area_obj:
            #Choose a small size but not necessarily the smallest to avoid sharp drops
            good_size = False      
            while good_size == False:
                rand_index = random.sample(range(nGrains), 1)
                if DiamV[rand_index] <= cbins[-2]*reduce_small:  #Small but not too small
                    good_size = True 
            #rand_index = min_G_index
        elif area_new > Area_obj:
            #print("#@#@##Area Exceedeed, increase area_tol##@##@#")
            rand_index = min_G_index
        
        grainsArea[1+i] =  AreaV[rand_index]
        grainsDiam[1+i] =  DiamV[rand_index]
        grainsThick[1+i] = thick_obj
        area_new += grainsArea[1+i]
        mass_new += rho_ice * grainsArea[1+i] * grainsThick[1+i]
        morph_new[1+i] = int(og_morphs[rand_index])
        for j in range(nother_bins+1):
            if grainsDiam[1+i] > cbins[0]:
                Ncum_curve[0] += 1
                break
            elif grainsDiam[1+i] <= cbins[-2]:
                Ncum_curve[-2] += 1
                break
            else:
                if grainsDiam[1+i] > cbins[j+1] and grainsDiam[1+i] <= cbins[j]:
                    Ncum_curve[j+1] += 1
                    break
    ave_thick_new = np.mean(grainsThick)
    
    # #Modify Thickness (Normal Thickness Distribution Tip)
    # min_adm = 0.05 #Avoid zero thickness
    # mu_th, sigma_th = 0.5, 0.1 # mean and standard deviation
    # srand = max(np.random.normal(mu_th, sigma_th), min_adm)
    # ThickVTX.append(srand) ##Around 0.5 m thickness
    # #ThickVTX.append(random.random()+1) ##Around 1 m thickness
    
    #Get Cumulative new for real
    log_cbins = []
    N_Ncum_curve = []
    for i in range(len(cbins)):
        if grainsDiam[i] <= 0:
            log_cbins.append(0.0)
        else:
            log_cbins.append(math.log(cbins[i]))
        if Ncum_curve[i] <= 0:
            N_Ncum_curve.append(0.0)
        else:
            N_Ncum_curve.append(abs(Ncum_curve[i]))
    #Convert to cumulative sum
    N_Ncum_curvep = N_Ncum_curve
    N_Ncum_curve = acc_sum(N_Ncum_curvep)
    
    log_Ncum_curve = []
    for idx in range(len(N_Ncum_curve)):
        if N_Ncum_curve[idx] <=0:
            log_Ncum_curve.append(0.0)
        else:
            log_Ncum_curve.append(math.log(N_Ncum_curve[idx]))
    
    #Filter zero results if occurring
    log_cbins0 = log_cbins
    log_Ncum_curve0 = log_Ncum_curve
    [log_cbins, log_Ncum_curve] = zeros_filter_slope(log_cbins0, log_Ncum_curve0)
    slope_new  = abs(gsd_slope_linear(objective, log_cbins, log_Ncum_curve))
    D_max_slope_new = log_cbins[0]
    D_min_slope_new = log_cbins[-1]
    N_max_slope_new = log_Ncum_curve[0]
    N_min_slope_new = log_Ncum_curve[-1]
    D_mean_slope_new = 0.5*(math.exp(D_max_slope_new) - math.exp(D_min_slope_new))
    Brute_Slope_new = (N_min_slope_new-N_max_slope_new) / (D_max_slope_new - D_min_slope_new)    

    #Verify compliance, otherwise redo until satisfy
    
    tolc_area = Area_obj * 0.01 #0.01 for 1.1, 0.025 for 1.15, 0.025 for 0.93, 0.025 for 0.88
    tolc_mass = mass_obj * 0.01
    tolc_slope = slope_obj * 0.025
    
    if abs(mass_new - mass_obj) < tolc_mass and abs(area_new - Area_obj) < tolc_area and abs(slope_new - slope_obj) < tolc_slope:
        compliance = True
        print("Compliance True after iteration # : ", niter_comp)
    else:
        niter_comp += 1
        if niter_comp % 200 == 0:
            print("Already tried to fit FSD this number of times: ", niter_comp)
            print(" OBJ / SIM / COND : ", abs(mass_new - mass_obj), " / ", tolc_mass, " / ",  abs(mass_new - mass_obj) < tolc_mass, abs(area_new - Area_obj), " / ", tolc_area, " / ",  abs(area_new - Area_obj) < tolc_area, abs(slope_new - slope_obj), " / ", tolc_slope, " / ", abs(slope_new - slope_obj) < tolc_slope)
        if niter_comp % 2000 == 0:
            coeff_small += 0.02
            print("Coeff Small increased to: ", coeff_small)
#Tranform to logarithm

#Get new slope

print("***Compare Outputs***")
print("Mass OBJ / SIM: ", mass_obj, mass_new)  #Might be related to a while loop
print("Area OBJ / SIM: ", Area_obj, area_new)  #Might be related to a while loop
print("Thickness OBJ / SIM: ", thick_obj, ave_thick_new) 
print("Slope OBJ / SIM: ", slope_obj, slope_new)   #Might be related to a while loop
print("Brute Slope OBJ / SIM: ", slope_obj, Brute_Slope_new) 
print("LOG Dmax_slope_new: ", D_max_slope_new)
print("LOG Dmin_slope_new: ", D_min_slope_new)
print("LOG Nmax_slope_new: ", N_max_slope_new)
print("LOG Nmin_slope_new: ", N_min_slope_new)
print("Dmax_slope_new: ", math.exp(D_max_slope_new))
print("Dmin_slope_new: ", math.exp(D_min_slope_new))
print("Nmax_slope_new: ", math.exp(N_max_slope_new))
print("Nmin_slope_new: ", math.exp(N_min_slope_new))
print("Dmean_slope_new: ", D_mean_slope_new)
print("Brute Slope_new: ", Brute_Slope_new)
print("Floe Bins: ", cbins)
print("Floe Count: ", Ncum_curve)
print("Cum Floe Count: ", N_Ncum_curve)
#print(morph_new)

#exit(1)



#Use AreaV and DiamV plus a new ThickV to fulfill objectives given and then make files
#With floe List let's now build our new initial conditions
morphV_save = np.zeros(nGrains_new+1)
morphV_save[0] = nGrains_new 
for i in range(nGrains_new):
    morphV_save[i+1] = morph_new[i] 

#Define other input files
ThickVTX = grainsThick
or_temp = -1.8
TempVTX = np.ones(nGrains_new) * or_temp
VelVTX = np.zeros((nGrains_new,3))

#Define extents for position
if year == 2020:
    offset_val = 380 #2020
elif year == 2018:
    offset_val = 400 #2018

x1 = 0
x2 = int(offset_val)
y1 = 0
y2 = int(offset_val)
    
#Using intersection assign positions to new list of Morphs avoiding overlaps (also save image to figDir with name)
pos_tol = float(Dmax * 0.5)
posRotVTX = np.zeros((nGrains_new,3))

#Reference configuration points for all morphs
pointsDat_new = []
for i in range(nGrains_new):
    pointsDat_new.append(pointsDat[int(morph_new[i])])

#Absolute configuration points for intersection and plotting
pointsControl = []
for i in range(nGrains_new):
    inter = True #Assume Not Good
    #Assume we can get it well given low populated space (WARNING!!!)
    if i % 100 == 0:
        print("Already on grain: ", str(i))
    
    while (inter == True):
        niterlimit = 0
        #Choose Random Position
        tol = int(pos_tol)
        
        positionP1 = np.random.choice(range(x1+tol, x2-tol))
        positionP2 = np.random.choice(range(y1+tol, y2-tol))

        positionP = [positionP1, positionP2]
        
        Rot = 0.000
        
        positionPFull = [positionP1, positionP2, Rot]
        
        pointTest = []
        pointList = pointsDat_new[i]
        for j in range(npoints):
            floatP0 = [float(pointList[j][0]), float(pointList[j][1])]
            #Angle rotation
            ppoint1 = floatP0[0]*math.cos(Rot) + floatP0[1]*math.sin(Rot)
            ppoint2 = -floatP0[0]*math.sin(Rot) + floatP0[1]*math.cos(Rot)
            floatP = [ppoint1, ppoint2]
            #Displace to final position for checking intersection
            Xg = floatP[0] + positionP[0]
            Yg = floatP[1] + positionP[1]
            vec = [Xg,Yg]
            pointTest.append(vec)

        #Check if there is intersection for all point Sets
        p1 = shapely.geometry.Polygon(pointTest)
        inter = intersection_grains(p1, len(pointsControl), pointsControl)
         
        if inter == False:
            #print("Insert more grains: ", rand_N)
            
            #append Grain to PointsDat list
            pointsControl.append(pointTest)
            
            #Increase positions
            posRotVTX[i][0] = positionPFull[0]
            posRotVTX[i][1] = positionPFull[1]
            posRotVTX[i][2] = positionPFull[2]
            
        else:
            niterlimit = niterlimit + 1
            #print("Tried ", niterlimit, " tries.")    
    
code_name = "MoreUniform2"
#code_name = "LessUniform1"
inputPath = inputNew + "SeaIceWGJ"+str(year)+"_"+code_name+"/"
if os.path.exists(inputPath) == False:
    os.mkdir(inputPath)

#Write new input files from final result
np.savetxt(inputPath+'initialpositions.dat', posRotVTX, fmt='%.5f %.5f %.5f')
np.savetxt(inputPath+'initialvelocities.dat', VelVTX, fmt='%.5f %.5f %.5f')
np.savetxt(inputPath+'initialtemper.dat', TempVTX, fmt='%.5f')
np.savetxt(inputPath+'initialthick.dat', ThickVTX, fmt='%.5f')
np.savetxt(inputPath+'morphologies.dat', morphV_save, fmt='%d')

#Generate other files (environmental data)
shutil.copyfile(inputDir + "Temp_data.dat", inputPath + "Temp_data.dat")
shutil.copyfile(inputDir + "WaveH_data.dat", inputPath + "WaveH_data.dat")
shutil.copyfile(inputDir + "WaveH_data_Day.dat", inputPath + "WaveH_data_Day.dat")
shutil.copyfile(inputDir + "WaveLength_data.dat", inputPath + "WaveLength_data.dat")
shutil.copyfile(inputDir + "WaveLength_data_Day.dat", inputPath + "WaveLength_data_Day.dat")
shutil.copyfile(inputDir + "Vel_Change.dat", inputPath + "Vel_Change.dat")
shutil.copyfile(inputDir + "SimInputParamsv2.dat", inputPath + "SimInputParamsv2.dat")

Fine_Conc = fine_conc
Coarse_Conc = area_new*100/Domain_Area
conc_text = [(Fine_Conc, Coarse_Conc)]
np.savetxt(inputPath+'Initial_Concentration.dat', conc_text, fmt='%.5f %.5f')

print("Plot Image") 
#FOR PLOTTING and Visualization
# Collect in patches
patches = []
#print(pointsDat)
#print("npts0: ", len(pointsDat))
#pointIndiv = np.split(pointsDat, nGrains)
#print("npts: ", len(pointIndiv))
#print("npts2: ", len(pointIndiv[0]))
pointsDat = pointsControl
for n in range(len(pointsDat)):
    #print(grains[n]._points)
    poly = Polygon(pointsDat[n], True)
    patches.append(poly)

#PLOT OF INITIAL STATE
fig, ax = plt.subplots(figsize = (15,15))
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
ax.autoscale_view()

pCol = PatchCollection(patches, facecolors='white',edgecolors='black', lw=0.1)
ax.set_facecolor('xkcd:blue')
ax.add_collection(pCol)

for i in range(nGrains_new):
    plt.text(posRotVTX[i][0], posRotVTX[i][1], str(morphV_save[i+1]), size = 8)

plt.savefig(figDir + 'InitialPWGJ'+str(year)+'_'+code_name+'.png', format='png')
plt.close(fig)

exit(1)


####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################






# File1 = "initialpositions.dat"
# posFile = np.loadtxt(inDir1 + File1, delimiter = " ")


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
