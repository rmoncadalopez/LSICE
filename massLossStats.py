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

#Folders for output flexibility
#Main directory for visualization

#Choose year for Analysis
#year = 2020;
#year = 2018;
year = int(sys.argv[3])


if year == 2020:
    main_dir = "./V_P_2020/" #2020
elif year == 2018:
    #main_dir = "./V_P_2018Prob/" #2018 
    main_dir = "./V_P_2018/" #2018 
    #main_dir = "./V_P/" #2018 
    #main_dir = "./V_P_2018h/" #2018 

create_folder = True

caseStart = int(sys.argv[1])
caseEnd = int(sys.argv[2])

if create_folder:
    folder_title = sys.argv[4]
    mainVfolder = main_dir + "MassLossAnalysis_" + str(folder_title) + "/"
    #mainVfolder = "./Visualizaton_" + str(folder_title) + "/"
    
    #Check if folder already exists
    if os.path.exists(mainVfolder) == False:
        #Create folders
        os.mkdir(mainVfolder)

#Main Ouput Vector
bkg_melt_ratioV = np.zeros((caseEnd-caseStart+1,7))
failed_case = []

for j in range(caseEnd-caseStart+1):
    if year == 2020:
        mainOutputfolder = "./Output/SeaIce2020_" + str(j+caseStart) + "/" #2020
        if os.path.exists(mainOutputfolder) == False:
            bkg_melt_ratioV[j][0] = j+caseStart
            failed_case.append(j+caseStart)
            print("No Folder")
            continue
    elif year == 2018:
        #mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #2018
        mainOutputfolder = "./Output/SeaIce2018_" + str(j+caseStart) + "/" #2018
        #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/" #2018
        #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/" #2018
        if os.path.exists(mainOutputfolder) == False:
            bkg_melt_ratioV[j][0] = j+caseStart
            failed_case.append(j+caseStart)
            print("No Folder")
            continue
    else:
        mainVfolder = main_dir    
        #mainVfolder = "./Visualization/"
        mainOutputfolder = "./Output/SeaIce/"
    
    
    #List of parameters
    testParamsFile = mainOutputfolder+"testParams0.dat"
    if os.path.exists(testParamsFile) == False:
        bkg_melt_ratioV[j][0] = j+caseStart
        failed_case.append(j+caseStart)
        print("No testParamsFile")
        continue
    testParamsV = np.loadtxt(testParamsFile)
    if len(testParamsV) == 0:
        bkg_melt_ratioV[j][0] = j+caseStart
        failed_case.append(j+caseStart)
        print("No testParamsFile correct length")
        continue
    nT = testParamsV[0]
    nTemp = testParamsV[1]
    nBreak = testParamsV[2]
    qvert = testParamsV[3]
    Khor = testParamsV[4]
    Qatm = testParamsV[5]
    Aatm = testParamsV[6]
    Batm = testParamsV[7]
    alpha_ice = testParamsV[8]
    alpha_ocean = testParamsV[9]
    limitMaxDiam = testParamsV[10]
    limitMinDiam = testParamsV[11]
    
    #Get number of steps
    numbergFile        =  mainOutputfolder+"numberg0.dat"
    if os.path.exists(numbergFile) == False:
        bkg_melt_ratioV[j][0] = j+caseStart
        failed_case.append(j+caseStart)
        print("No numbergFile")
        continue
    numberg            = np.loadtxt(numbergFile)
    nSteps             = int(len(numberg))
    print(j+caseStart)

    #Get Mass Loss File from Main_Periodic_Ocean
    massLossFile = mainOutputfolder + "massLoss0.dat"
    if os.path.exists(massLossFile) == False:
        bkg_melt_ratioV[j][0] = j+caseStart
        failed_case.append(j+caseStart)
        print("No massLossFile")
        continue
    massLossV = np.loadtxt(massLossFile)
    if len(massLossV) < nSteps:
        bkg_melt_ratioV[j][0] = j+caseStart
        failed_case.append(j+caseStart)
        print("No massLossFile good length")
        continue
    
    
    #Extract value from SIM Output file (Shape = time_plot, tmass, tmasscoarse, tmassfines, loss_mcl, loss_mcv, gain_fines, loss_fines, total_melt_coarse, tmassA, tmasscoarseA )
    time = np.zeros(nSteps) # Time (Days)
    tmass = np.zeros(nSteps) # Total Mass _mass (kg from here on)                     
    tmasscoarse = np.zeros(nSteps) # Total Mass Coarse _mass
    tmassfines = np.zeros(nSteps) # Tota Mass Fines
    loss_bkg = np.zeros(nSteps) # Loss Melt Coarse Lateral
    loss_mcv = np.zeros(nSteps) # Loss Melt Coarse Vertical (Thickness)
    loss_bc = np.zeros(nSteps) # Loss Breakage Coarse
    total_melt_coarse = np.zeros(nSteps) # Total Mass Loss Coarse due to Melt
    tmassA = np.zeros(nSteps) # Total Mass from Points
    tmasscoarseA = np.zeros(nSteps) # Total Mass Coarse from Points
    loss_fines = np.zeros(nSteps)
    gain_fines = np.zeros(nSteps)
    fine_net_loss = np.zeros(nSteps)
    dt_size = 86400 #1 day of snapshot or 86400 s
    
    for i in range(nSteps):
        time[i]              = massLossV[i][0]
        tmass[i]             = massLossV[i][1]
        tmasscoarse[i]       = massLossV[i][2]
        tmassfines[i]        = massLossV[i][3]
        loss_bkg[i]          = massLossV[i][4] 
        loss_mcv[i]          = massLossV[i][5] 
        gain_fines[i]        = massLossV[i][6] 
        loss_fines[i]        = massLossV[i][7] 
        total_melt_coarse[i] = massLossV[i][8] 
        tmassA[i]            = massLossV[i][9]
        tmasscoarseA[i]      = massLossV[i][10]
        fine_net_loss[i]     = massLossV[i][7] - massLossV[i][6] 
    
    
    bkg_melt_ratioV[j][0] = j+caseStart
    bkg_melt_ratioV[j][1] = nBreak
    bkg_melt_ratioV[j][2] = qvert
    bkg_melt_ratioV[j][3] = Qatm
    bkg_melt_ratioV[j][4] = loss_bkg[-1]
    bkg_melt_ratioV[j][5] = loss_mcv[-1]
    
    bkg_melt_ratio = loss_bkg[-1]/loss_mcv[-1]
    bkg_melt_ratioV[j][6] = bkg_melt_ratio
    
np.savetxt(mainVfolder+'mass_loss_analysis_'+str(caseStart)+"_"+str(caseEnd)+".dat", bkg_melt_ratioV, fmt='%d %.5f %.5f %.5f %.5f %.5f %.5f')   
print("No data for # cases: ", len(failed_case))
print(failed_case)
    
    