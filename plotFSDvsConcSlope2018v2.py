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
from scipy.integrate import simps
from numpy import trapz

from PIL import Image

#import shapely.geometry

import math
import datetime
import csv
import cv2
import numpy as np
import glob

# objective function for slope Linear
def objective(x, a, b):
    return a * x + b


#Last day coarse
def lastdayfind(timeV, concB):
    last_day = 0
    for i in range(len(timeV)):
        if concB[i] > 0.01:
            last_day += 1
        else:
            break   
    return last_day

#GSD Slope Linear 
def slope_linear(objective, x, y):
    # curve fit for all results
    poptMean, _ = curve_fit(objective, x, y)
    # summarize the parameter values
    aMean = poptMean[0]
    bMean = poptMean[1]
    return aMean

#GSD Slope Linear Removing Flat Parts
def slope_linear(objective, x, y):
    xd = []
    yd = []
    for i in range(len(x)):
        if y[i] > 0:
            xd.append(x[i])
            yd.append(y[i])
            
    # curve fit for all results
    poptMean, _ = curve_fit(objective, xd, yd)
    # summarize the parameter values
    aMean = poptMean[0]
    bMean = poptMean[1]
    return aMean

#GSD Slope Linear Removing Flat Parts
def slope_linear_max(objective, x, y):

    sample_size = 20
    #Clean dataset to use only decaying part and cut flat part (for that get gradient)
    slope_y = np.gradient(y)
    abs_slope_y = np.abs(slope_y)
    max_vals0 = sorted(abs_slope_y, key = lambda x:float(x))
    max_vals = max_vals0[::-1]
    print("GRADIENT SAMPLE")
    print(max_vals)
    max_top = max_vals[0:sample_size] #Choose sample_size top values of slope
    max_slope = sum(max_top)/sample_size
    #max_slope = np.max(np.abs(slope_y))
    return max_slope

#GSD Slope Linear Removing Flat Parts
def slope_linearA(objective, x, y):

    #Clean dataset to use only decaying part and cut flat part (for that get gradient)
    slope_y = np.gradient(y)
    return np.mean(slope_y)

#GSD Slope Linear Removing Flat Parts
def slope_linear_corrected(objective, x, y):

    #Clean dataset to use only decaying part and cut flat part (for that get gradient)
    slope_y = np.gradient(y)
    print("GRADIENT SAMPLE")
    print(slope_y)
    print("Len y:", len(y))
    print("Len dy:", len(slope_y))

    tol = 0.155 #0.105 #Slopes around 0.1 or more seem to be more representative of signficant decline, behavior changes witht his value too low gets odd
    #0.0105 represents coarse better, 0.205 fines 

    #FORM A If three consecutive points keep a very low slope we will cut here as it means ice is no more (probably)
    # real_len = len(y)
    # for i in range(len(y)-2):
    #     if abs(slope_y[i]) < tol and abs(slope_y[i+1]) < tol and abs(slope_y[i+2]) < tol: 
    #         real_len = i+1

    # #Define new arrays to get slope
    # x_f = np.zeros(real_len)
    # y_f = np.zeros(real_len)
    # for i in range(real_len):
    #     x_f[i] = x[i]
    #     y_f[i] = y[i]


    #FORM B (ELIMINATE BELOW TOLERANCE POINTS and reconstruct and find slope of new curve)
    y_f = []
    x_f = []
    ct = 0
    for i in range(len(y)):
        if abs(slope_y[i]) >= tol:
            x_f.append(ct)
            y_f.append(y[i])
            ct = ct + 1



    # curve fit for all results
    poptMean, _ = curve_fit(objective, x_f, y_f)
    # summarize the parameter values
    aMean = poptMean[0]
    bMean = poptMean[1]
    return aMean
    
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

def slope_int(timeCT, Conc):
    slopeVint = 0
    for i in range(len(Conc)-1):
        slopeVint += abs(timeCT[i+1] - timeCT[i]) * 0.5*(Conc[i]+Conc[i+1])  #area = trapz(y, dx=1) area = simps(y, dx=1)        
    #return slopeVint
    return simps(Conc, dx=1)
    
def slope_int_data(timeCT, Conc):
    slopeVint = 0
    for i in range(len(Conc)-1):
        slopeVint += abs(timeCT[i+1] - timeCT[i]) * 0.5*(Conc[i]+Conc[i+1])  #area = trapz(y, dx=1) area = simps(y, dx=1)        
    return slopeVint
    #return simps(Conc, dx=1)
    
def gen_init_conc(caseNo, year, nDays):

    if year == 2018:
        if Prob:
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/"
            #mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
            mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #PROB CHANGE
            #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/"
        else:
            mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
    if year == 2020:
        #mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" 
        mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" 
    testParamsFile             = mainOutputfolder+"testParams0.dat"
    concFile = mainOutputfolder+"normConc0.dat"
    concV  = np.loadtxt(concFile)
    testParamsV      = np.loadtxt(testParamsFile)

    #List of parameters
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
        
    c0C = ConcB[0]
    c0F = ConcF[0]
    c0T = Conc[0]

    return c0C, c0F, c0T 

def gen_values_per_case(caseNo, year, ndays):
    Prob = False
    if year == 2018:
        if Prob:
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/"
            #mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
            mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #PROB CHANGE
            #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/"
        else:
            mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
    if year == 2020:
        #mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" 
        mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" 
    testParamsFile             = mainOutputfolder+"testParams0.dat"
    concFile = mainOutputfolder+"normConc0.dat"
    concV  = np.loadtxt(concFile)
    testParamsV      = np.loadtxt(testParamsFile)

    #List of parameters
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

    #Save Conc. and Normalized Area
    timeCB = []
    timeCF = []
    timeCT = []
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
        timeCB.append(concV[i][0])
        timeCF.append(concV[i][0])
        timeCT.append(concV[i][0])
        nArea.append(concV[i][1])
        Conc.append(concV[i][2]*100)
        ConcB.append(concV[i][3]*100)
        ConcF.append(concV[i][4]*100)
        surfArea.append(concV[i][5])
        aveThick.append(concV[i][6])
        globalOceanTemp.append(concV[i][7])
        aveThickC.append(concV[i][8])
        aveThickF.append(concV[i][9])

    # #Remove zero extreme value
    # timeCB0 = timeC
    # timeCF0 = timeC
    # timeCT0 = timeC
    # ConcB0 = ConcB
    # ConcF0 = ConcF
    # ConcT0 = Conc
    # [timeCB, ConcB] = zeros_filter_slope(timeCB0, ConcB0)
    # [timeCF, ConcF] = zeros_filter_slope(timeCF0, ConcF0)
    # [timeCT, Conc] = zeros_filter_slope(timeCT0, ConcT0)

    #Get Slope
    #slope_coarse = slope_linear_corrected(objective, timeCB, ConcB) 
    #slope_fine = slope_linear_corrected(objective, timeCF, ConcF) 
    #slope_total = slope_linear_corrected(objective, timeCT, Conc) 

    #slope_coarse = slope_linearC(objective, timeCB, ConcB) 
    #slope_fine = slope_linearC(objective, timeCF, ConcF) 
    #slope_total = slope_linearC(objective, timeCT, Conc) 

    #slope_coarse = slope_linearA(objective, timeCB, ConcB) 
    #slope_fine = slope_linearA(objective, timeCF, ConcF) 
    #slope_total = slope_linearA(objective, timeCT, Conc) 
    
    slope_coarse = slope_int(timeCB, ConcB)/ndays
    slope_fine = slope_int(timeCF, ConcF)/ndays
    slope_total = slope_int(timeCT, Conc)/ndays
    
    slopeCoarseVInt = slope_int(timeCT, Conc)/ndays
    last_day = lastdayfind(timeCB, ConcB)

    #slope_coarse = slope_linear_max(objective, timeC, ConcB) 
    #slope_fine = slope_linear_max(objective, timeC, ConcF) 
    #slope_total = slope_linear_max(objective, timeC, Conc) 

    #slope_coarse = slope_linear(objective, timeCB, ConcB) 
    #slope_fine = slope_linear(objective, timeCF, ConcF) 
    #slope_total = slope_linear(objective, timeCT, Conc) 

    return nBreak, qvert, abs(slope_coarse), abs(slope_fine), abs(slope_total), ConcB, surfArea, slopeCoarseVInt, last_day

def get_slope(slope_dataFolder):
    simFile = slope_dataFolder + 'Slope_sim.dat'
    simFileV = np.loadtxt(simFile)
    vecT = []
    slopeV = np.zeros(len(simFileV))
    for i in range(len(simFileV)):
        vecT.append(simFileV[i][0])
        slopeV[i] = max(simFileV[i][1] - 1, 0) #To match FND which is the usual value, we use CFND given gaps in floes due to breakage
    
    alphaV = max(slopeV[0] , 0)
    return alphaV, vecT, slopeV

def get_dataslope(slope_dataFolder):
    dataFile = slope_dataFolder + 'Slope_data.dat'
    dataFileV = np.loadtxt(dataFile)
    vecT = []
    slopeV = np.zeros(len(dataFileV))
    for i in range(len(dataFileV)):
        vecT.append(dataFileV[i][0])
        slopeV[i] = max(dataFileV[i][1] - 1, 0) #To match FND which is the usual value, we use CFND given gaps in floes due to breakage
    
    alphaV = max(slopeV[0] , 0)
    return alphaV, vecT, slopeV
    
def mass_process(caseNo, year):
    if year == 2020:
        mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" #2020
    elif year == 2018:
        #mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #2018
        mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/" #2018
        #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/" #2018
        #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/" #2018

    #Get Mass Loss File from Main_Periodic_Ocean
    massLossFile = mainOutputfolder + "massLoss0.dat"
    massLossV = np.loadtxt(massLossFile)

    #Get number of steps
    numbergFile        =  mainOutputfolder+"numberg0.dat"
    numberg            = np.loadtxt(numbergFile)
    nSteps             = int(len(numberg))
    
    if len(massLossV) >= nSteps-1:
        status=True
    else:
        status=False
        return 0, 0, 0, status

    #Extract value from SIM Output file (Shape = time_plot, tmass, tmasscoarse, tmassfines, loss_mcl, loss_mcv, gain_fines, loss_fines, total_melt_coarse, tmassA, tmasscoarseA )
    time = np.zeros(nSteps) # Time (Days)
    tmass = np.zeros(nSteps) # Total Mass _mass (kg from here on)                     
    tmasscoarse = np.zeros(nSteps) # Total Mass Coarse _mass
    tmassfines = np.zeros(nSteps) # Tota Mass Fines
    loss_mcl = np.zeros(nSteps) # Loss Melt Coarse Lateral
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
        loss_mcl[i]          = massLossV[i][4] 
        loss_mcv[i]          = massLossV[i][5] 
        gain_fines[i]        = massLossV[i][6] 
        loss_fines[i]        = massLossV[i][7] 
        total_melt_coarse[i] = massLossV[i][8] 
        tmassA[i]            = massLossV[i][9]
        tmasscoarseA[i]      = massLossV[i][10]
        fine_net_loss[i]     = massLossV[i][7] - massLossV[i][6] 
        
    pre_mcl = []
    pre_mcv = []
    pre_total = []
    pre_time = []
    tstop = 0 #Count stop time to figure rate. If truncated tstop < nDays, else it will update with nDays
    for i in range(len(total_melt_coarse)):
        if total_melt_coarse[i] > 0:
            pre_time.append(i)
            pre_mcl.append(loss_mcl[i])
            pre_mcv.append(loss_mcv[i])
            pre_total.append(total_melt_coarse[i])
            tstop += 1
        else:
            pre_time.append(i)
            pre_mcl.append(0)
            pre_mcv.append(0)
            pre_total.append(1)
            tstop += 1
        #To avoid flat segments (mostly)
        tol = 0.01*total_melt_coarse[-1]
        if i < len(total_melt_coarse)-2:
            if abs(total_melt_coarse[i] - total_melt_coarse[i+1]) < tol and abs(total_melt_coarse[i] - total_melt_coarse[i+1]) < tol and abs(total_melt_coarse[i] - total_melt_coarse[i+1]) < tol:
                break
    dataV = np.divide(pre_mcl, pre_total)[1:] #Make it pretty removing first point
    dataV2 = np.divide(pre_mcv, pre_total)[1:]
    dataV3 = total_melt_coarse
    pre_time = pre_time[1:]
    pre_time = np.array(pre_time)
    
    if tstop > 0:
        mass_slope =  total_melt_coarse[-1] / tstop #len(pre_time)
    else:
        mass_slope = 0.0
    
    areaBkg = simps(dataV, dx=1)
    areaMelt = simps(dataV2, dx=1)
    tMassC = total_melt_coarse[-1]
    dataV3 = total_melt_coarse/tMassC/(nSteps-1)
    #tMassC = pre_total[-1]
    #dataV3 = pre_total/tMassC/(len(pre_total))
    ratioInt = areaBkg/areaMelt
    ratio = (dataV[-1]/dataV2[-1])
    cum_mass_loss = simps(dataV3, dx = 1)
    
    return ratio, ratioInt, cum_mass_loss, status, mass_slope

#For average alpha of satellite data
def data_flat(dataV0):
    dataV = dataV0
    ave = 0
    count = 0
    for i in range(len(dataV)):
        if dataV[i] > 0.00:  #Remove zero values that bias the results, we just want values before ice disappears, otherwise it doesn't make sense
            ave += dataV[i]
            count += 1
    
    if len(dataV) > 0 and count > 0:
        count = count
    else:
        count = 1
    return ave/count

#For average alpha of simulation data    
def mean_flat(dataV0, i):
    dataV = np.array(dataV0)[i]
    ave = 0
    count = 0
    for i in range(len(dataV)):
        if dataV[i] > 0.00:  #Remove zero values that bias the results, we just want values before ice disappears, otherwise it doesn't make sense
            ave += dataV[i]
            count += 1

    if len(dataV) > 0 and count > 0:
        count = count
    else:
        count = 1
    return ave/count
    
#Find average concentration slope
def conc_slope(dataV):
    tstop = 0
    tol = 0.5 #Define a concentration tolerance we will define as zero (inspect its effect)
    Conc_init = dataV[0]
    Conc_end = dataV[0]
    for i in range(len(dataV)):
        if dataV[i] > tol:
            tstop += 1
            Conc_end = dataV[i]
        else:
            break #Assume non-recoverable sea ice

    if tstop == 0:
        return 0
    
    #return tstop*nDays 
    return abs((Conc_init-Conc_end)/tstop)
    
def conc_slope_data(dataV, timeV):
    tstop = 0
    tol = 0.5 #Define a concentration tolerance we will define as zero (inspect its effect)
    Conc_init = dataV[0]
    Conc_end = dataV[0]
    for i in range(len(dataV)):
        if dataV[i] > tol:
            tstop = timeV[i]
            Conc_end = dataV[i] + 1
        else:
            break #Assume non-recoverable sea ice

    if tstop == 0:
        return 0
    
    #return tstop*nDays 
    return abs((Conc_init-Conc_end)/tstop)

########################################################################################################################################################################################

#alpha = [3.5, 3.1, 3.2, 3.8, 4.1]
cases = [1036, 1035, 315, 1034]
cases = [1035, 315, 1034]
#cases = [315, 210, 280, 245]
#cases = [1036, 1035, 315, 1034, 1033]

#B All cases With extremes
casesB = [ 315,   210,  245,  280,     624, 628,  884, 885, 886, 887, 888, 889, 890,  891,  892 ]
Btag   = [ 1000, 2000, 1700, 1500, 1000000,   1,  200, 250, 400, 500, 600, 700, 800, 1250, 2500 ]
#B All cases No extreme cases
casesB = [ 315,   210,  245,  280, 884, 885, 886, 887, 888, 889, 890,  891,  892 ]
Btag   = [ 1000, 2000, 1700, 1500, 200, 250, 400, 500, 600, 700, 800, 1250, 2500 ]
#B filter cases (qv = 25) (For clearer results) (S = 310)
casesB = [884, 888,  315,  210, 624 ]
Btag   = [200, 600, 1000, 2000, 10000000 ]
invBtag = [1/200, 1/600, 1/1000, 1/2000, 1/10000000]

casesB = [884,  315,  210, 624 ]
Btag   = [200, 1000, 2000, 10000000 ] #5
invBtag = [1/200, 1/1000, 1/2000, 1/10000000] #5

casesB = [315, 210, 280, 245]
casesB = [315, 210, 280, 887, 245]
Btag   = [1000, 2000, 1500,  1700]
Btag   = [1000, 2000, 1500,  500, 1700]
invBtag = [1/1000, 1/2000, 1/1500, 1/1700]
invBtag = [1/1000, 1/2000, 1/1500, 1/500, 1/1700]

# Btag   = [200, 1000 ] #3
# invBtag = [1/200, 1/1000] #3
# casesB = [884, 888,  315,  210, 624]
# Btag   = [200, 600, 1000, 2000, 10000000]
# invBtag = [1/200, 1/600, 1/1000, 1/2000, 1/10000000]



#qv filter cases (B = 1000) (S = 310)
cases_qv = [319, 318, 317, 316, 315, 314, 313,  513, 604]
qv_tag  =  [5,    10,  15,  20,  25,  30,  35,   40,  45]
cases_qv = [318, 316, 315, 314, 313]
qv_tag  =  [10,  20,  25,  30,  35]
cases_qv = [317, 315,  313,  513, 604]
cases_qv = [318,  316,  314,  513 ]
qv_tag  =  [10,    20,   30,   40 ]

# #qv filter cases (B = 2500) (S = 310)
# cases_qv = [671, 722, 773, 827, 892, 957, 1022,  518, 609]
# qv_tag  =  [5,    10,  15,  20,  25,  30,   35,   40,  45]
# cases_qv = [671, 722, 773, 827, 892, 957,  518, 609]
# qv_tag  =  [5,    10,  15,  20,  25,  30,   40,  45]
# cases_qv = [722, 827, 957, 518 ]
# qv_tag  =  [10,  20,  30,   40 ]

#qv filter cases (B = 40000) (S = 310)
cases_qv = [1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355 ]
qv_tag =   [   5,   10,   15,   20,   25,   30,   35,   40,   45,   70,  100,  200 ]
cases_qv = [1344,  1346, 1348, 1350,  1352, 1353, 1354, 1355 ]
qv_tag =   [   5,    15,   25,   35,    45,   70,  100,  200 ]
cases_qv = [1344,   1348,    1354,  1355 ] #5
qv_tag =   [   5,     25,      100,   200 ] #5
# cases_qv = [1344,    1354 ] #3
# qv_tag =   [   5,     100 ] #3


#Comparison for movemebt
cases_no_mov =   [1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387]
cases_mov = [884,  885,  886,  887,  1324, 1336, 1348, 1360]
#High Freq
cases_no_mov =   [1380, 1381, 1382, 1383]
cases_mov = [884,  885,  886,  887]
#Low Freq
#cases_no_mov =   [1384, 1385, 1386, 1387]
#cases_mov = [1324, 1336, 1348, 1360]



movB =         [200,  250,  400,  500, 10000, 20000, 40000,  10000000]

inv_movB =  [1/200,  1/250,  1/400,  1/500, 1/10000, 1/20000, 1/40000,  1/10000000]
inv_movB =  [1/200,  1/250,  1/400,  1/500] #High Freq
#inv_movB =  [1/10000, 1/20000, 1/40000,  1/10000000] #Low Freq


year = 2018

outputFolder = "./Output/"

if year == 2020:
    main_dir = "./V_P_2020/" #2020
    #main_dir = "./V_P_2020old/" #2020
elif year == 2018:
    #main_dir = "./V_P/" #2018 
    main_dir = "./V_P_2018/"
    #main_dir = "./V_P_2018Prob/"  #2018 #PROB. CHANGE
    #main_dir = "./V_P_2018h/" #2018 

Prob = False
if year == 2018:
    ndays = 48
    if Prob:
        #outputFolder2 = "./Output/Analysis/"
        #outputFolder2 = "./Output/Analysis_2018/"
        outputFolder2 = "./Output/AnalysisFSDSlope_2018Prob/" #PROB CHANGE
        #outputFolder2 = "./Output/Analysis_2018h/"
    else:
        outputFolder2 = "./Output/AnalysisFSDSlope_2018_test/"
elif year == 2020:
    ndays = 24
    outputFolder2 = "./Output/AnalysisFSDSlope_2020/"

if os.path.exists(outputFolder2) == False:
	os.makedirs(outputFolder2)


[c0C, c0F, c0T] = gen_init_conc(casesB[0], year, ndays)

alphaV = np.zeros(len(cases))
alphaVB = np.zeros(len(casesB)) #Not really useful, just alpha at t = 0
alphaV_qv = np.zeros(len(cases_qv))
nBreakV = np.zeros(len(cases))
qvertV  = np.zeros(len(cases))
slope_coarseV  = np.zeros(len(cases))
slope_fineV    = np.zeros(len(cases))
slope_totalV   = np.zeros(len(cases))
slopeCoarseVInt = np.zeros(len(cases))
VECT = []
SLOPEV = []
VECTB = []
VECT_qv = []
SLOPEVB = []
SLOPEV_qv = []
CONC_CV = []
RATIOV = []
RATIOVInt = []
RATIOVIntB = []
RATIOVInt_qv = []
cumML = []
mass_slopeV = []
SURF_AREA = []
alphaAveB = []
alphaAve_qv = []
init_surf = np.zeros(len(cases))

alphaInit = []
alphaAve = []
concRateV = []
lastdayV = []

#Mov effect
cum_loss_mov = np.zeros(len(cases_mov))
mass_slope_mov = np.zeros(len(cases_mov))
ConcB_mov = np.zeros(len(cases_mov))
cum_loss_no_mov = np.zeros(len(cases_no_mov))
mass_slope_no_mov = np.zeros(len(cases_no_mov))
ConcB_no_mov = np.zeros(len(cases_no_mov))

#Obtain satellite concentration data
if year == 2018:
    inputFolder = "./Input/SeaIceWGJ2018/"
    dataGSDFile = inputFolder + "GSD_FinalResultsData_2018v2.dat"  #Source: MODIS Data #v2 modified with different days
elif year == 2020:
    #Field Data from Input (2020)
    inputFolder = "./Input/SeaIceWGJ2020/"
    dataGSDFile = inputFolder + "GSD_FinalResults2020.dat" #Source: MODIS Data #v1 modified with different days 2020

dataConcV = np.loadtxt(dataGSDFile)
timeCData = np.zeros(len(dataConcV))
dataConcC  = np.zeros(len(dataConcV))

for i in range(len(timeCData)):
    timeCData[i]    = dataConcV[i][0]
    dataConcC[i]    = dataConcV[i][5]
#data_conc_rate = conc_slope_data(dataConcC, timeCData)
#Data mean concentration
data_conc_rate = slope_int_data(timeCData, dataConcC) /ndays

#Closest RMSE approximation to mass loss
mainVfolder2file_RMSE = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(315) + "/" + "BkMeltratio.dat"
ratioV_RMSE = np.loadtxt(mainVfolder2file_RMSE)
dataRatioInt = ratioV_RMSE[1]

for i in range(len(cases)):
    if i > 0:
        mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases[i]) + "/"
        mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases[i]) + "/" + "BkMeltratio.dat"
    else:
        mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases[i]) + "/"
        mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases[i]) + "/" + "BkMeltratio.dat"
        slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
        [alphaOr, vecOr, slopeOr] = get_dataslope(slope_dataFolder)
        data_init_alpha = slopeOr[0]
        data_ave_alpha = data_flat(slopeOr)
        
    ratioV = np.loadtxt(mainVfolder2file)
    RATIOV.append(ratioV[0])
    RATIOVInt.append(ratioV[1])
    cumML.append(ratioV[2])
    [a, b, aa, d, mass_slopeVi] = mass_process(cases[i], year)
    mass_slopeV.append(mass_slopeVi)
    #Could add mass slope here as  mass_slope.append(ratioV[3])
    slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
    [alphaV[i], vecT, slopeV] = get_slope(slope_dataFolder) #From slope file
    alphaInit.append(slopeV[0])
    alphaAve.append(data_flat(slopeV))
    VECT.append(vecT)
    SLOPEV.append(slopeV)
    [nBreakV[i], qvertV[i], slope_coarseV[i], slope_fineV[i], slope_totalV[i], ConcB, surfArea, slopeCoarseVInt[i], last_day] = gen_values_per_case(cases[i], year, ndays) #From Conc File
    if i == 2:
        ConcB[0] = CONC_CV[0][0]
        ConcB[1] = CONC_CV[0][1]
    CONC_CV.append(ConcB)
    #concRateV.append(conc_slope(ConcB))  slope_coarseV[i]
    concRateV.append(slope_coarseV[i])
    SURF_AREA.append(surfArea)
    lastdayV.append(last_day)
    init_surf[i] = surfArea[0]
    
#Adjust for errors in image segmentation
slopeOr[0] = 2.59 
slopeOr[1] = 2.5


for i in range(len(casesB)):
    if i > 0:
        mainVfolderB = main_dir + "Visualization_" + "Test_qv_nb_alterB_" + str(casesB[i]) + "/"
        mainVfolder2fileB = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(casesB[i]) + "/" + "BkMeltratio.dat"
        slope_dataFolderB = mainVfolderB+'Output_Comparison/GSD_Compare/'
    else:
        mainVfolderB = main_dir + "Visualization_" + "Test_qv_nb_alterB_" + str(casesB[i]) + "/"
        mainVfolder2fileB = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(casesB[i]) + "/" + "BkMeltratio.dat"
        slope_dataFolderB = mainVfolderB+'Output_Comparison/GSD_Compare/'
        [alphaOrB, vecOrB, slopeOrB] = get_dataslope(slope_dataFolderB)
    #print(mainVfolder2fileB)
    ratioVB = np.loadtxt(mainVfolder2fileB)
    RATIOVIntB.append(ratioVB[1])
    [alphaVB[i], vecTB, slopeVB] = get_slope(slope_dataFolderB) #From slope file for breakage impact on alpha
    VECTB.append(vecTB)
    SLOPEVB.append(slopeVB)
    #alphaAveB.append(np.mean(SLOPEVB))
    alphaAveB.append(mean_flat(SLOPEVB, i))

for i in range(len(cases_qv)):
    if i > 0:
        mainVfolder_qv = main_dir + "Visualization_" + "Test_qv_nb_alterqv_" + str(cases_qv[i]) + "/"
        mainVfolder2file_qv = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_qv[i]) + "/" + "BkMeltratio.dat"
        slope_dataFolder_qv = mainVfolder_qv+'Output_Comparison/GSD_Compare/'
    else:
        mainVfolder_qv = main_dir + "Visualization_" + "Test_qv_nb_alterqv_" + str(cases_qv[i]) + "/"
        mainVfolder2file_qv = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_qv[i]) + "/" + "BkMeltratio.dat"
        slope_dataFolder_qv = mainVfolder_qv+'Output_Comparison/GSD_Compare/'
        [alphaOr_qv, vecOr_qv, slopeOr_qv] = get_dataslope(slope_dataFolder_qv)
    print(mainVfolder2file_qv)
    ratioV_qv = np.loadtxt(mainVfolder2file_qv)
    RATIOVInt_qv.append(ratioV_qv[1])
    [alphaV_qv[i], vecT_qv, slopeV_qv] = get_slope(slope_dataFolder_qv) #From slope file for breakage impact on alpha
    VECT_qv.append(vecT_qv)
    SLOPEV_qv.append(slopeV_qv)
    #alphaAve_qv.append(np.mean(SLOPEV_qv))
    alphaAve_qv.append(mean_flat(SLOPEV_qv, i))

#Mov, No Mov
for i in range(len(cases_mov)):
    #if i > 0:
    #   mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases_mov[i]) + "/"
        #mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_mov[i]) + "/" + "BkMeltratio.dat"
    #else:
    #    mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases_mov[i]) + "/"
        #mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_mov[i]) + "/" + "BkMeltratio.dat"
        #slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
        #[alphaOr, vecOr, slopeOr] = get_dataslope(slope_dataFolder)
    #ratioV = np.loadtxt(mainVfolder2file)
    #RATIOV.append(ratioV[0])
    #RATIOVInt.append(ratioV[1])
    #cumML.append(ratioV[2])
    #Could add mass slope here as  mass_slope.append(ratioV[3])
    #slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
    #[alphaV[i], vecT, slopeV] = get_slope(slope_dataFolder) #From slope file
    #VECT.append(vecT)
    #SLOPEV.append(slopeV)
    #print(mainVfolder)
    [a, b, c, d, ConcB_mov[i], e, f, g, h] = gen_values_per_case(cases_mov[i], year, ndays) #From Conc File
    [a, b, cum_loss_mov[i], d, mass_slope_mov[i]] = mass_process(cases_mov[i], year)
    # if i == 2:
    #     ConcB[0] = CONC_CV[0][0]
    #     ConcB[1] = CONC_CV[0][1]
    # CONC_CV.append(ConcB)
    # SURF_AREA.append(surfArea)
    # init_surf[i] = surfArea[0]

for i in range(len(cases_no_mov)):
    #if i > 0:
    #    mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases_no_mov[i]) + "/"
        #mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_no_mov[i]) + "/" + "BkMeltratio.dat"
    #else:
    #    mainVfolder = main_dir + "Visualization_" + "Test_qv_nb_alter_" + str(cases_no_mov[i]) + "/"
        #mainVfolder2file = main_dir + "MassLoss_" + "Test_qv_nb_ML_alter_" + str(cases_no_mov[i]) + "/" + "BkMeltratio.dat"
        #slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
        #[alphaOr, vecOr, slopeOr] = get_dataslope(slope_dataFolder)
    #ratioV = np.loadtxt(mainVfolder2file)
    #RATIOV.append(ratioV[0])
    #RATIOVInt.append(ratioV[1])
    #cumML.append(ratioV[2])
    #Could add mass slope here as  mass_slope.append(ratioV[3])
    #slope_dataFolder = mainVfolder+'Output_Comparison/GSD_Compare/'
    #[alphaV[i], vecT, slopeV] = get_slope(slope_dataFolder) #From slope file
    #VECT.append(vecT)
    #SLOPEV.append(slopeV)
    [a, b, c, d, ConcB_no_mov[i], e, f, g, h] = gen_values_per_case(cases_no_mov[i], year, ndays) #From Conc File
    [a, b, cum_loss_no_mov[i], d, mass_slope_no_mov[i]] = mass_process(cases_no_mov[i], year)
    # if i == 2:
    #     ConcB[0] = CONC_CV[0][0]
    #     ConcB[1] = CONC_CV[0][1]
    # CONC_CV.append(ConcB)
    # SURF_AREA.append(surfArea)
    # init_surf[i] = surfArea[0]


print("Coarse: ", slope_coarseV)
print("Slopes: ", alphaV)

#Plot Slopes (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, slope_coarseV, 'r-*', markersize=15)
plt.axvline(x=alphaOr, color='b', label ='Field Value')
plt.legend()
plt.title('Uniformity Variation vs. Coarse Conc. Rate', size=15)
plt.xlabel('Uniformity Coefficient (FSD Slope)', size=20)
plt.ylabel('Concentration Rate (%/d)', size=20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "FSDSlopevsConcCoarse.png")
plt.close()

#Plot Slopes (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, slopeCoarseVInt, 'r-*', markersize=15)
#plt.legend()
plt.title('FSD Uniformity vs. Cumulative Coarse Concentration', size=20)
plt.xlabel('FSD Slope', size=20)
plt.ylabel('Cumulative Coarse Concentration', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSDSlopevsConcCoarseInt.png")
plt.close()

#Plot Slopes (Coarse, Fine, Total)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, slope_coarseV, 'r*', label='Coarse')
plt.plot(alphaV, slope_fineV, 'g*', label='Fine')
plt.plot(alphaV, slope_totalV, 'k*', label='Total')
plt.axvline(x=alphaOr, color='b', label='Field Value')
plt.legend()
plt.title('FSD Uniformity vs. Concentration Rate', size=15)
plt.xlabel('Uniformity Coefficient (FSD Slope)')
plt.ylabel('Concentration Rate (%/d)')
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "FSDSlopevsConcAll.png")
plt.close()

#Plot Mass Slope (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, mass_slopeV, 'r-*', markersize=15)
plt.title('FSD Uniformity vs. Mass Loss Rate', size=20)
plt.xlabel('FSD Slope', size=20)
plt.ylabel('Mass Loss Rate (kg/day)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSDSlopevsMassSlope.png")
plt.close()

tagInitAlpha = ['~2.0', '~2.5', '~3.0']

#Plot Slope Evolution Graph to Compare
fig = plt.figure(figsize=(10,10))
print(vecOr)
print(slopeOr)
plt.plot(vecOr, slopeOr, 'k--', label = 'Satellite Obs.')
for i in range(len(cases)):
    #plt.plot(VECT[i], SLOPEV[i], '*-', label = 'α: ' + str(round(alphaV[i],2)))
    plt.plot( VECT[i], SLOPEV[i], '*-', label = 'Initial α: ' + tagInitAlpha[i] )
plt.legend(fontsize = 18)
plt.title('FSD Slope (α) over time vs. Initial α', size=20)
plt.ylabel('FSD Slope (α)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSD_Init_SlopeTime.png")
plt.close()


def errorDS(dataX, dataY, simX, simY):
    sum = 0
    length = len(dataX)
    for i in range(len(dataY)):
        for j in range(len(simY)):
            if dataX[i] == simX[j]:
                sum += (simY[i]-dataY[i])**2
    return (sum/length)**0.5


#Plot Slope Evolution Graph to Compare (for BREAKAGE values)
fig = plt.figure(figsize=(10,10))
plt.plot(vecOr, slopeOr, 'ks-', label = 'Satellite Observations')

errorVB = np.zeros(len(casesB))
for i in range(len(casesB)):
    errorVB[i] = errorDS(vecOr, slopeOr, VECTB[i], SLOPEVB[i])

for i in range(len(casesB)):
    #print(i)
    #print(SLOPEVB[i])
    #plt.plot(VECTB[i], SLOPEVB[i], '*-', label = 'Break Frequency: ' + str(round(invBtag[i],9)))
    #plt.plot(VECTB[i], SLOPEVB[i], '*-', label = 'B: ' + str(round(invBtag[i]*24*3600,9)) + ' 1/d')
    #if i == 0 or i == 1:
    plt.plot(VECTB[i], SLOPEVB[i], '*-', label = 'Break Frequency: ' + str(round(invBtag[i]*3600,2)) + ' 1/h RMSE Data: ' + str(round(errorVB[i],3)) )
    print("Bstep:", round(Btag[i]))

plt.legend(fontsize = 18)
plt.title('FSD Uniformity vs. Break Frequency', size=20)
plt.ylabel('Uniformity (α)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSD_SlopeTime_BkgStep.png")
plt.close()

print('B vs. FSD Data')
print('Cases: ', casesB)
for i in range(len(casesB)):
    print(i, casesB[i], invBtag[i])
    print(VECTB[i])
    print(SLOPEVB[i])

#Plot Slope Evolution Graph to Compare (for qv Vertical Melt values)
fig = plt.figure(figsize=(10,10))
plt.plot(vecOr, slopeOr, 'ks-', label = 'Satellite Observations')
for i in range(len(cases_qv)):
    print(i)
    print(SLOPEV_qv[i])
    #plt.plot(VECT_qv[i], SLOPEV_qv[i], '*-', label = 'qv: ' + str(round(qv_tag[i],0)) + ' W/m2/°C')
    #if i == 0 or i == 2:
    plt.plot(VECT_qv[i], SLOPEV_qv[i], '*-', label = 'Melt Rate: ' + str(round(qv_tag[i],0)) + ' W/m2/°C')
plt.legend(fontsize = 18)
plt.title('FSD Uniformity vs. Melt Rate', size=20)
plt.ylabel('Uniformity (α)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSD_SlopeTime_qvMelt.png")
plt.close()


#Plot Slope Evolution Graph to Compare (Melt Break Compare)
#Avoid super high values due to error
for i in range(len(SLOPEVB[2])):
    SLOPEVB[2][i] = min(SLOPEVB[2][i],4.4) 

fig = plt.figure(figsize=(10,10))
plt.plot(vecOr, slopeOr, 'k--', label = 'Satellite Observations')
plt.plot(VECTB[2], SLOPEVB[2], 'k*-', label = 'Obs. Fit Sim.')
plt.plot(VECTB[0], SLOPEVB[0], 'b*-', label = 'High Break Sim.')
plt.plot(VECT_qv[1], SLOPEV_qv[1], 'r*-', label = 'Low Break Sim.')
plt.plot(VECTB[3], SLOPEVB[3], 'g*-', label = 'Higher Break Sim.')
plt.legend(fontsize = 18)
plt.title('FSD Slope vs. Break and Melt Effect', size=20)
plt.ylabel('FSD slope (α)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSD_SlopeTime_BreakvsMelt.png")
plt.close()

print("TIME Ref: ", VECTB[2])
print("Slope Ref: ", SLOPEVB[2])

print('qv vs. FSD Data')
print('Cases: ', cases_qv)
for i in range(len(cases_qv)):
    print(i, cases_qv[i], qv_tag[i])
    print(VECT_qv[i])
    print(SLOPEV_qv[i])

#Plot Surface Area Evolution Graph to Compare
fig = plt.figure(figsize=(10,10))
for i in range(len(cases)):
    plt.plot(VECT[i], SURF_AREA[i], '*-', label = 'α: ' + str(round(alphaV[i],2)))
plt.legend(fontsize = 18)
plt.title('Surface Area versus Uniformity', size=20)
plt.ylabel('Surface Area (km2)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "SurfArea_SlopeTime.png")
plt.close()

#Plot Conc Evolution Graph to Compare
fig = plt.figure(figsize=(10,10))
plt.plot(timeCData, dataConcC, 'ks-', label = 'Observations')
for i in range(len(cases)):
    plt.plot(VECT[i], CONC_CV[i], '*-', label = 'α: ' + str(round(alphaV[i],2)))
plt.legend(fontsize = 18)
plt.title('Coarse Concentration over Time Depending on FSD Slope', size=20)
plt.ylabel('Concentration (%)', size=20)
plt.xlabel('Time (days)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "Conc_SlopeTime.png")
plt.close()

#Plot Initial Surface Area
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, init_surf, 'r-*', markersize=15)
plt.title('Init. Surf versus Slope Coeff.', size=15)
plt.xlabel('Uniformity Coefficient (FSD Slope)')
plt.ylabel('Surface Area (km2)')
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "FSDSlopevsSurfArea.png")
plt.close()

#Plot Ratio Effect (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, RATIOV, 'r*', markersize=15, label = 'Simulation Values')
plt.title('Coarse Floe Breakage Decay vs. Vertical Melt Ratio, FSD Slope Effect', size=15)
plt.xlabel('Uniformity Coefficient (FSD Slope)')
plt.ylabel('Breakage Decay / Vertical Melt Ratio')
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
# plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "FSDSlopevsCoarseRatio.png")
plt.close()

#Plot Ratio Effect Cumulative (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, RATIOVInt, 'rs', markersize=16, label='Simulations')
plt.plot(data_init_alpha, dataRatioInt, 'k*', markersize=20, label='Best fit for Sat. Observations')
plt.axhline(y = 1.0, color = 'r', linestyle = '--', linewidth=4)
plt.title('Initial FSD Uniformity vs. Break / Melt Dominance', size=20)
plt.xlabel('Initial Uniformity (α)', size=20)
plt.ylabel('Breakage / Melt Ratio', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(fontsize = 18)
plt.savefig(outputFolder2 + "FSDSlopevsBMCoarseRatio.png")
plt.close()


print('alphaV: ', alphaV)
print('RATIOVInt: ', RATIOVInt)
print('Sat Data')
print(data_init_alpha)
print(dataRatioInt)


#Plot Ratio Effect Cumulative by different Breakage Frequencies (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(invBtag, RATIOVIntB, 'rs', markersize=20, label = 'Simulation Values')
plt.title('Breakage Step vs. Break/Melt Dominance', size=20)
plt.xlabel('Breakage Step (1/s)', size=20)
plt.ylabel('Breakage/Melt Ratio', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.legend(fontsize = 18)
plt.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "BStepvsCumCoarseRatio.png")
plt.close()

#Ave. Alpha slope by different Breakage Frequencies (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(invBtag, alphaAveB, 'rs', markersize=16, label='Simulations')
plt.axhline(y = data_ave_alpha, color = 'k', linestyle = '--', label='Satellite Observations', linewidth=4)
plt.title('Mean FSD Uniformity vs. Break Frequency', size=20)
plt.xlabel('Break Frequency - B (1/s)', size=20)
plt.ylabel('Mean Uniformity (α)', size=20)
plt.ylim([2.00,4.05])
plt.legend(fontsize = 18)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "BStepvsAveAlpha.png")
plt.close()

#Plot Ratio Effect Cumulative by different qv (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(qv_tag, RATIOVInt_qv, 'r-*', markersize=17)
plt.title('Vertical Melt Rate vs. Break/Melt Dominance', size=20)
plt.xlabel('Vertical Melt Rate', size=20)
plt.ylabel('Breakage/Melt Ratio', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "qvMeltvsCumCoarseRatio.png")
plt.close()

#Ave. Alpha slope by different qv (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(qv_tag, alphaAve_qv, 'rs', markersize=16, label='Simulations')
plt.axhline(y = data_ave_alpha, color = 'k', linestyle = '--', label='Satellite Observations', linewidth=4)
plt.title('Mean FSD Uniformity vs. Melt Rate', size=20)
plt.xlabel('Melt Rate - qv (W/m2/°C)', size=20)
plt.ylabel('Mean Uniformity (α)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.legend(fontsize = 18)
plt.savefig(outputFolder2 + "qvMeltvsAveAlpha.png")
plt.close()


#Ave. Alpha slope by different initial alpha (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaInit, alphaAve, 'rs', markersize=16, label='Simulations')
plt.plot(data_init_alpha, data_ave_alpha, 'k*', markersize=20, label='Satellite Observations')
#plt.axhline(y = data_ave_alpha, color = 'k', linestyle = '--', label='Observations')
plt.title('Mean FSD Uniformity vs. Initial Uniformity', size=20)
plt.xlabel('Initial Uniformity (α)', size=20)
plt.ylabel('Mean Uniformity (α)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(fontsize = 18)
plt.savefig(outputFolder2 + "InitAlphavsAveAlpha.png")
plt.close()

print('AlphaInitCases: ', cases)
print('AlphaInit: ', alphaInit)
print(len(alphaInit))
print('concRateV: ', concRateV)
print('Last day:  ', lastdayV)
print('Sat Data')

#Adjustment
data_init_alpha = alphaInit[1]
print(data_init_alpha, data_ave_alpha)

concRateVN = []
alphaInitN = alphaInit
num = (c0C**2)
for i in range(len(alphaInit)):
    concRateVN.append( num/(2*concRateV[i]*ndays) )
    alphaInitN[i] = alphaInit[i]
    
print('AlphaInitN: ', alphaInitN)
print('concRateVN: ', concRateVN)

#For supplementary material initial alpha
#Conc. Rate versus Initial Alpha
fig = plt.figure(figsize=(10,10))
plt.plot(alphaInitN[:-1], concRateVN[:-1], 'rs', markersize=16, label='Simulations')
plt.plot(alphaInitN[-1], concRateVN[-1], 'bs', markersize=16, label='Simulations: Best Fit')
plt.plot(data_init_alpha, (dataConcC[0]**2)/(2*data_conc_rate*ndays), 'k*', markersize=20, label='Satellite Observations')
#plt.axhline(y = (dataConcC[0]**2)/(2*data_ave_alpha), color = 'k', linestyle = '--', label='Observations')
plt.title(r'Initial FSD $\alpha$ vs. Mean Coarse Concentration Rate', size=20)
plt.xlabel(r'Initial FSD Slope ($\alpha$)', size=20)
plt.ylabel(r'Mean Coarse Concentration Rate (% $d^{-1}$)', size=20)
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(fontsize = 18)
plt.savefig(outputFolder2 + "InitAlphavsConcRateS.png")
plt.close()

#Plot Ratio Effect Cumulative by different qv (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaInit[:-1], RATIOVInt[:-1], 'rs', markersize=16, label='Simulations')
plt.plot(alphaInit[-1], RATIOVInt[-1], 'bs', markersize=16, label='Simulations: Best Fit')
plt.title(r'Initial FSD $\alpha$ vs. Breakage / Melt Ratio', size=20)
plt.xlabel(r'Initial FSD Slope ($\alpha$)', size=20)
plt.ylabel(r'$\mu$', size=20)
plt.axhline(y = 1, color = 'k', linestyle = '--')
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "InitAlphavsCoarseRatioS.png")
plt.close()


#Plot Cumulative Mass Loss Ave (Coarse)
fig = plt.figure(figsize=(10,10))
plt.plot(alphaV, cumML, 'r-*', markersize=17)
plt.title('FSD Uniformity vs. Cumulative Mass Loss', size=20)
plt.xlabel('Uniformity Coefficient (FSD Slope)', size=20)
plt.ylabel('Cumulative Coarse Mass Loss', size=20)
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig(outputFolder2 + "FSDSlopevsCumMassLoss.png")
plt.close()


#Movement and Concentration Area
fig = plt.figure(figsize=(10,10))
plt.plot(inv_movB, ConcB_mov, 'r*', markersize=15, label = 'Movement')
plt.plot(inv_movB, ConcB_no_mov, 'b.', markersize=15, label = 'No Movement')
plt.title('Movement and Concentration Area', size=20)
plt.xlabel('Breakage Freq (1/s)', size=20)
plt.ylabel('Concentration Area', size=20)
plt.legend()
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "Mov_ConcArea_highfreqB.png")
plt.close()

#Movement and Cumulative Mass Loss
fig = plt.figure(figsize=(10,10))
plt.plot(inv_movB, cum_loss_mov, 'r*', markersize=15, label = 'Movement')
plt.plot(inv_movB, cum_loss_no_mov, 'b.', markersize=15, label = 'No Movement')
plt.title('Movement and Cumulative Mass Loss', size=20)
plt.xlabel('Breakage Freq (1/s)', size=20)
plt.ylabel('Cumulative Mass Loss', size=20)
plt.legend()
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "Mov_CumMassLoss_highfreqB.png")
plt.close()

#Movement and Mass Slope
fig = plt.figure(figsize=(10,10))
plt.plot(inv_movB, mass_slope_mov, 'r*', markersize=15, label = 'Movement')
plt.plot(inv_movB, mass_slope_no_mov, 'b.', markersize=15, label = 'No Movement')
plt.title('Movement and Mass Slope', size=20)
plt.xlabel('Breakage Freq (1/s)', size=20)
plt.ylabel('Mass Slope (kg/d)', size=20)
plt.legend()
# plt.xaxis.label.set_size(20)
# plt.yaxis.label.set_size(20)
plt.tick_params(labelsize=15)
plt.savefig(outputFolder2 + "Mov_MassSlope_highfreqB.png")
plt.close()
