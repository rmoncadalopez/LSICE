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
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

from PIL import Image

#import shapely.geometry

import math
import datetime
import csv
import cv2
import numpy as np
import glob

from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

class RemainderFixed(axes_size.Scaled):
    def __init__(self, xsizes, ysizes, divider):
        self.xsizes =xsizes
        self.ysizes =ysizes
        self.div = divider

    def get_size(self, renderer):
        xrel, xabs = axes_size.AddList(self.xsizes).get_size(renderer)
        yrel, yabs = axes_size.AddList(self.ysizes).get_size(renderer)
        bb = Bbox.from_bounds(*self.div.get_position()).transformed(self.div._fig.transFigure)
        w = bb.width/self.div._fig.dpi - xabs
        h = bb.height/self.div._fig.dpi - yabs
        return 0, min([w,h])


def make_square_axes_with_colorbar(ax, size=0.1, pad=0.1):
    """ Make an axes square, add a colorbar axes next to it, 
        Parameters: size: Size of colorbar axes in inches
                    pad : Padding between axes and cbar in inches
        Returns: colorbar axes
    """
    divider = make_axes_locatable(ax)
    margin_size = axes_size.Fixed(size)
    pad_size = axes_size.Fixed(pad)
    xsizes = [pad_size, margin_size]
    yhax = divider.append_axes("right", size=margin_size, pad=pad_size)
    divider.set_horizontal([RemainderFixed(xsizes, [], divider)] + xsizes)
    divider.set_vertical([RemainderFixed(xsizes, [], divider)])
    return yhax


qv_val = True #True use qv, #False use kappa


# objective function for slope Linear
def objective(x, a, b):
    return a * x + b



#GSD Slope Linear Removing Flat Parts
def slope_linear(objective, x, y):
    # curve fit for all results
    poptMean, _ = curve_fit(objective, x, y)
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
    
def slope_int(Conc):
    #slopeVint = 0
    #for i in range(len(Conc)-1):
        #slopeVint += abs(timeCT[i+1] - timeCT[i]) * 0.5*(Conc[i]+Conc[i+1])  #area = trapz(y, dx=1) area = simps(y, dx=1)        
    #return slopeVint
    return simps(Conc, dx=1)


def slope_simple(Conc, nDays):
    tstop = 0
    tol = 1 #Define a concentration tolerance we will define as zero (inspect its effect)
    Conc_init = Conc[0]
    Conc_end = Conc[0]
    for i in range(len(Conc)):
        if Conc[i] > tol:
            tstop += 1
            Conc_end = Conc[i]
        else:
            break #Assume non-recoverable sea ice

    if tstop == 0:
        return 0
    
    #return tstop*nDays 
    return abs((Conc_init-Conc_end)/tstop)*nDays 

Prob = False

def read_ave_velocity(caseNo, year):
    if year == 2018:
        if Prob:
            inputFile = "./V_P_2018Prob/Visualization_Test_qv_nb_" + str(caseNo) + "/velocity_ave.dat" #PROB CHANGE
        else:
            inputFile = "./V_P_2018/Visualization_Test_qv_nb_" + str(caseNo) + "/velocity_ave.dat"
    elif year == 2020:   
        inputFile = "./V_P_2020/Visualization_Test_qv_nb_" + str(caseNo) + "/velocity_ave.dat" 
    dataV = np.loadtxt(inputFile)
    
    ave_vel = 0
    len_med = int(len(dataV)*0.5)+1
    
    for i in range(len_med):
        ave_vel += dataV[i][1]

    if len(dataV) < 1:
        ave_vel = 0
    else:
        ave_vel = ave_vel / len_med
        
    return ave_vel

def read_errors(caseNo, year):
    if year == 2018:
        #inputFile = "./V_P/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat"
        #inputFile = "./V_P_2018/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat"
        if Prob:
            inputFile = "./V_P_2018Prob/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat" #PROB CHANGE
        else:#inputFile = "./V_P_2018h/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat"
            inputFile = "./V_P_2018/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat"
    elif year == 2020:   
        inputFile = "./V_P_2020/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/Param_Errors_" + str(caseNo) + ".dat" 
    dataV = np.loadtxt(inputFile)

    return dataV[4], dataV[5], dataV[6], dataV[7] #C_error, F_error, T_error, Temp_error

def read_params(caseNo, year):
    print(caseNo)	 
    if year == 2018:
        if Prob:
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/"
            #mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
            mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #PROB CHANGE
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + " copy/"   #Just temporal
        else:
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/"
            mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
            #mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #PROB CHANGE
            #mainOutputfolder = "./V_P/Visualization_Test_qv_nb_" + str(caseNo) + "/Output_Comparison/"   #Watch out for missprinting
            #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/"
    elif year == 2020:
        mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/"
        #mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/"
    testParamsFile             = mainOutputfolder+"testParams0.dat"
    testParamsV      = np.loadtxt(testParamsFile)

    if len(testParamsV)<1:
        return 0, 0, 0
    
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

    
    if qv_val:
        return nBreak, qvert, Qatm #For qv
    else:
        return nBreak, Khor, Qatm #For kHor

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
        mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" 
        #mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" 
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

def gen_values_per_case(caseNo, year, nDays):

    if year == 2018:
        if Prob:
            #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/"
            #mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
            mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #PROB CHANGE
            #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/"
        else:
            mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "/"
    if year == 2020:
        mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" 
        #mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" 
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

    # slope_coarse = slope_linear_corrected(objective, timeC, ConcB) 
    # slope_fine = slope_linear_corrected(objective, timeC, ConcF) 
    # slope_total = slope_linear_corrected(objective, timeC, Conc) 

    #Inverse cumulative concentration
    # slope_coarse = 1/slope_int(ConcB) 
    # slope_fine   = 1/slope_int(ConcF) 
    # slope_total  = 1/slope_int(Conc)
    
    #Cumulative Concentration
    slope_coarse = slope_int(ConcB) 
    slope_fine   = slope_int(ConcF) 
    slope_total  = slope_int(Conc)    
    
    #Slope of loss, simply put Conc Lost / days to lose it (units %/d)
    #slope_coarse = slope_simple(ConcB, nDays) 
    #slope_fine   = slope_simple(ConcF, nDays) 
    #slope_total  = slope_simple(Conc, nDays) 

    #Linear fitting techniques
    #slope_coarse = slope_linear_max(objective, timeC, ConcB) 
    #slope_fine = slope_linear_max(objective, timeC, ConcF) 
    #slope_total = slope_linear_max(objective, timeC, Conc) 

    # slope_coarse = slope_linear(objective, timeC, ConcB) 
    # slope_fine = slope_linear(objective, timeC, ConcF) 
    # slope_total = slope_linear(objective, timeC, Conc) 

    if qv_val:
        return nBreak, qvert, slope_coarse, slope_fine, slope_total #qv
    else:
        return nBreak, Khor, slope_coarse, slope_fine, slope_total #kappa

def plot_image(xc, yc, xf, yf, x_label, y_label, title, image_name, year, path):
    fig, ax1 = plt.subplots(figsize=(10,10))
    color = 'tab:red'
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax1.plot(xc, yc, 'r*', label='Coarse')
    ax1.plot(xf, yf, 'g*', label='Fine')
    plt.legend()
    plt.title(title)
    if year == 2018:
        if Prob:
            #plt.savefig("./Output/Analysis/" + image_name +".png", format='png')
            #plt.savefig("./Output/Analysis_2018/" + image_name +".png", format='png')
            plt.savefig(path + image_name +".png", format='png') #PROB CHANGE
            #plt.savefig("./Output/Analysis_2018h/" + image_name +".png", format='png')
        else:
            plt.savefig(path + image_name +".png", format='png')
    elif year == 2020:
        plt.savefig(path + image_name +".png", format='png')    
    plt.close()


def mass_process(caseNo, year):
    if year == 2020:
        #mainOutputfolder = "./Output/SeaIce2020N_" + str(caseNo) + "/" #2020
        mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" #2020
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
        print("BAD mass output")
        return 0, 0, 0, status, 0

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
        mass_slope =  total_melt_coarse[-1] / tstop #len(pre_time) #Total mass lost /  time to lose it
    else:
        mass_slope = 0.0
        
    if mass_slope > 0.0:
        status=True
    else:
        status=False
        print("TOTAL BAD MELT COARSE:", total_melt_coarse)
        return 0, 0, 0, status, 0
    

    #Bkg/Melt Ratio
    #Form A. Using Cumulative
    areaBkg = simps(dataV, dx=1)
    areaMelt = simps(dataV2, dx=1)
    #tMassC = pre_total[-1]
    #dataV3 = pre_total/tMassC/(len(pre_total))
    ratioInt = areaBkg/areaMelt
    ratio = (dataV[-1]/dataV2[-1])
    #Form B. Using instantaneous
    rateBkg = np.gradient(dataV)
    rateMelt = np.gradient(dataV2)
    areaBkg = simps(rateBkg, dx=1)
    areaMelt = simps(rateMelt, dx=1)
    #ratioInt = areaBkg/areaMelt
    ratio = (dataV[-1]/dataV2[-1])
    
    print("MASS DATA OBS")
    print("pre_mcl: ", pre_mcl)
    print("pre_mcv: ", pre_mcv)
    print("dataV: ", dataV)
    print("dataV2: ", dataV2)
    print("rateBkg: ", rateBkg)
    print("rateMelt: ", rateMelt)
    print("total_melt_coarse: ", total_melt_coarse)
    
    #Cumulative Mass Loss
    tMassC = total_melt_coarse[-1]
    dataV3 = total_melt_coarse/tMassC/(nSteps-1)
    print("dataV3: ", dataV3)
    
    #Form A. Integrated Cum Mass Loss again
    cum_mass_loss = simps(dataV3, dx = 1)
    
    #Form B. Derivate Total Loss and then integrate it again
    #mass_loss_rate = np.gradient(dataV3)  #Get mass loss over time
    mass_loss_rate = np.gradient(total_melt_coarse)  #Get mass loss over time
    cum_mass_loss = simps(mass_loss_rate, dx = 1) / tstop  #Integrate to add mass lost and average how fast it took to lose it
    
    
    #print('Mass')
    
    return ratio, ratioInt, cum_mass_loss, status, mass_slope

#MAIN START ##########################################################################
start_case = int(sys.argv[1])
end_case = int(sys.argv[2])
year = int(sys.argv[3])
Satm = int(sys.argv[4])
caseL = (end_case-start_case)+1;
vec_export_text = []

outputFolder = "./Output/"
#year = 2020
#year = 2018

if year == 2018:
    if Prob:
        #outputFolder2 = "./Output/Analysis/"
        #outputFolder2 = "./Output/Analysis_2018/"
        outputFolder2 = "./Output/Analysis_2018Prob/" #PROB CHANGE
        #outputFolder2 = "./Output/Analysis_2018h/"
    else:
        #outputFolder2 = "./Output/Analysis_2018_S_"+str(Satm)+"/"
        if qv_val:
            outputFolder2 = "./Output/Analysis_2018Int_S_"+str(Satm)+"_Test/"
        else:
            outputFolder2 = "./Output/Analysis_2018Int_S_"+str(Satm)+"_kappa/"
elif year == 2020:
    outputFolder2 = "./Output/Analysis_2020/"

if os.path.exists(outputFolder2) == False:
	os.makedirs(outputFolder2)


if year == 2018:
    nDays = 48
    nDaysxy = 1  #For not adimensional days
elif year == 2020:
    nDays = 24 #Confirm!!!
    nDaysxy = 1  #For not adimensional days
	

Qatm_compare = Satm #Originally 290
#nBreakRange = np.linspace(400, 2000, 17)
#nBreakRange = [200, 400, 600, 800, 1000, 1250, 1500, 1700, 2000, 2500]
if Prob:
    nBreakRange = [9200, 9300, 9400, 9500, 9600, 9700]
else:
    nBreakRange = [200, 400, 600, 800, 1000, 1250, 1500, 1700, 2000, 2500]
    nBreakRange = [800,1000, 1250, 1500, 1700, 2000, 2500]
    if qv_val:
        nBreakRange = [800,1000, 1250, 1500, 1700] #FOCUS
        nBreakRange = [1500, 1700, 2000, 2500, 10000, 20000, 40000, 10000000] #Wider for exploring
        nBreakRange = [200, 250, 400, 500, 10000, 20000, 40000, 10000000] #Wider for exploring
        nBreakRange = [200, 250, 400, 500, 800, 1000, 1250, 1500, 2000, 2500, 10000, 20000, 40000, 10000000] #Wider for exploring more
        nBreakRange = [250, 400, 500, 800, 1000, 1250, 1500, 2000, 2500, 10000, 20000, 40000] #Wider for exploring more2
        nBreakRange = [400, 500, 800, 1000, 1250, 1500, 2000, 2500, 10000, 20000, 40000] 
	#nBreakRange = [200, 250, 400, 1000, 1250, 1500, 2000, 2500, 10000, 20000, 40000, 10000000] #Wider for exploring more (FIX)
        #nBreakRange = [10000, 20000, 40000, 10000000] #Wider for exploring #LOW freq
        #nBreakRange = [200, 250, 400, 500] #Wider for exploring #HIGh
    else:
        nBreakRange = [400, 600, 800,1000, 1250, 1500, 1700, 2000] #Kappa try more frequent values
        #nBreakRange = [2000, 6000, 10000, 20000, 40000, 10000000] #Kappa try #New less frequent values  
        #nBreakRange = [20000, 40000] #Kappa try #New less frequent values
    #nBreakRange = [1000, 1250, 1500, 1700] #Just try
    #nBreakRange = [2500, 2000, 1700, 1500, 1250, 1000, 800, 600, 400, 200]
if qv_val:
    #qvertRange = np.linspace(5,35,7)
    #qvertRange = [5, 10, 15, 20, 25, 30, 35, 40, 45]
    qvertRange = [15, 20, 25, 30, 35] #FOCUS
    qvertRange = [5, 10, 15, 20, 25, 30, 35, 40, 45, 70, 100, 200]  #Wider for exploring
    qvertRange = [40, 45, 70, 100, 200]  #Wider for exploring
    qvertRange = [5, 15, 25, 35, 40, 45, 70, 100, 200]  #Wider for exploring more
    qvertRange = [5, 15, 25, 35, 40, 45] #Wider for exploring more2
else:
    qvertRange = [10, 20, 50, 75, 100, 150, 200, 300, 400, 500] #Kappa Try  try more frequent values
    #qvertRange = [10, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000] #Kappa Try higher kappa
    #qvertRange = [10, 50, 100, 200, 500] #Kappa Try higher kappa
#qvertRange = [25, 30, 35] #Just try
#qvertRange = [10,15,20,25,30,35] #earlier more narrow range
#Satm_range = [240,250,260,270,290,310,320]

valid_cases0 = []
bad_cases = []
#Remove non-usable cases for Qatm
for i in range(caseL):
    caseNo = i + start_case
    if (caseNo <= 622 or caseNo>=630):
        if (caseNo < 1037 or caseNo > 1099) and caseNo < 1506:  #Modify for right cases
        	[nBreak, qvert, Qatm] = read_params(caseNo, year)
        	if Qatm == 0:
        	    print("AA")
        	    bad_cases.append(caseNo)
        	elif Qatm == Qatm_compare:
        	    valid_cases0.append(caseNo)

if len(bad_cases) > 0:
    print("Invalid or BAD cases: ", bad_cases)
    #print("Process stopped due to bad cases, results will be distorted!!!")
    #exit(1)

#For MultiD
Vec_MultiD = []

#Plot nB versus Dc/Dt

#Plot qv contours
fig, ax1 = plt.subplots(figsize=(10,10))
ax1.set_xlabel("nBreak")
ax1.set_ylabel("dC/dt")
plt.title("qv Impact")

nBBB = []
qvC = []
qvF = []
qvT= []

#Choose a value of qv first
refT = []
massProcessbad = []
for k in range(len(qvertRange)):
    qvert_compare = qvertRange[k]
    #Loop through all valid cases to get useful cases for this qv
    valid_cases1 = []
    for i in range(len(valid_cases0)):
        [nBreak, qvert, Qatm] = read_params(valid_cases0[i], year);
        if qvert == qvert_compare:
            valid_cases1.append(valid_cases0[i])

    if len(valid_cases1) < 1:
        continue

    #Loop through all useful cases to get nB values
    nB_vec = []
    Cslope_vec = []
    Fslope_vec = []
    Tslope_vec = []
    #Assure order to plot nBreakRange
    for j in range(len(nBreakRange)) :
        nBreak_compare = nBreakRange[j]
        for ii in range(len(valid_cases1)):
            if valid_cases1[ii] == 256:
                [nBreak1, qvert1, Qatm1] = read_params(valid_cases1[ii], year);
                print("TEST: ", nBreak1, qvert1, Qatm1)
            [nBreak, qvert, Qatm] = read_params(valid_cases1[ii], year);
            if nBreak == nBreak_compare:          
                [aa, bb, slope_coarse, slope_fine, slope_total] = gen_values_per_case(valid_cases1[ii], year, nDays)
                nB_vec.append(nBreak)
                Cslope_vec.append(abs(slope_coarse))
                Fslope_vec.append(abs(slope_fine))
                Tslope_vec.append(abs(slope_total))
                print(nBreak)
                print(qvert)
                [C_error, F_error, T_error, Temp_error] = read_errors(valid_cases1[ii], year)
                [ratio, ratioInt, cum_mass_loss, status_good, mass_slope] = mass_process(valid_cases1[ii], year)
                if status_good == False:
                    massProcessbad.append([valid_cases1[ii], qvert_compare, nBreak_compare])
                #multiV = [ qvert_compare, (nBreak_compare), abs(slope_coarse), abs(slope_fine), valid_cases1[ii], C_error, F_error, T_error, Temp_error, abs(slope_total) ];
                multiV = [qvert_compare, (86400 * 1/nBreak_compare), abs(slope_coarse), abs(slope_fine), valid_cases1[ii], C_error, F_error, T_error, Temp_error, abs(slope_total), ratio, ratioInt, cum_mass_loss, mass_slope];
                #multiV = [qvert, (1/nBreak), abs(slope_coarse), abs(slope_fine), valid_cases1[ii], C_error, F_error, T_error, Temp_error];
                refT.append([qvert, nBreak])
                Vec_MultiD.append(multiV)
                break
    print(qvert)
    print(nB_vec)
    print(Cslope_vec)
    plot_image(nB_vec, Cslope_vec, nB_vec, Fslope_vec, 'nBreak', 'dC/dt', 'Concentration slope versus nBreak for qvert = ' + str(qvert), 'C vs nB_qv_' + str(round(qvert,0)), year, outputFolder2)

    nBBB.append(nB_vec)
    qvC.append(Cslope_vec)
    qvF.append(Fslope_vec)
    qvT.append(Tslope_vec)
    ax1.plot(nB_vec, Cslope_vec, '*-', label='Coarse, qv: ' + str(qvert_compare) )
    ax1.plot(nB_vec, Fslope_vec, 'o--', label='Fines, qv: ' + str(qvert_compare) )
    ax1.plot(nB_vec, Tslope_vec, '.', label='Total, qv: ' + str(qvert_compare) )
    
plt.legend()    
plt.savefig(outputFolder2+"qv_impact_coarse_fine.png", format='png')    
plt.close()

print("Need mass again BAD: ")
print(massProcessbad)
#exit(1)

#EXTRA PLOT *************
# #qv Impact by Type
# #C
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("nBreak")
# ax1.set_ylabel("dC/dt")
# plt.title("qv Impact Coarse")

# print(qvertRange)
# print(nBBB)
# print(qvC)

# for k in range(len(qvertRange)):
#      ax1.plot(nBBB[k], qvC[k], '*-', label='Coarse, qv: ' + str(qvertRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"qv_impact_coarse.png", format='png')
# plt.close()

# #F
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("nBreak")
# ax1.set_ylabel("dC/dt")
# plt.title("qv Impact Fines")

# for k in range(len(qvertRange)):
#      ax1.plot(nBBB[k], qvF[k], 'o--', label='Fines, qv: ' + str(qvertRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"qv_impact_fines.png", format='png')
# plt.close()

# #T
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("nBreak")
# ax1.set_ylabel("dC/dt")
# plt.title("qv Impact Total")

# for k in range(len(qvertRange)):
#      ax1.plot(nBBB[k], qvT[k], '.-', label='Total, qv: ' + str(qvertRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"qv_impact_total.png", format='png')
# plt.close()
#EXTRA PLOT *************


#Plot qv versus Dc/Dt

#Plot nB contours
fig, ax1 = plt.subplots(figsize=(10,10))
ax1.set_xlabel("qv")
ax1.set_ylabel("dC/dt")
plt.title("nB Impact")

qVVV = []
nbC = []
nbF = []
nbT = []

print("Second part")
#Choose a value of nB first
for k in range(len(nBreakRange)):
    nBreak_compare = nBreakRange[k]
    #Loop through all valid cases to get useful cases for this nB
    valid_cases1 = []
    for i in range(len(valid_cases0)):
        [nBreak, qvert, Qatm] = read_params(valid_cases0[i], year);
        if nBreak == nBreak_compare:
            valid_cases1.append(valid_cases0[i])

    if len(valid_cases1) < 1:
        continue
 
    #Loop through all useful cases to get qv values
    qv_vec = []
    Cslope_vec = []
    Fslope_vec = []
    Tslope_vec = []
    #Assure order to plot qv
    for j in range(len(qvertRange)): 
        qvert_compare = qvertRange[j]
        for ii in range(len(valid_cases1)):
            [nBreak, qvert, Qatm] = read_params(valid_cases1[ii], year);
            if qvert == qvert_compare:      
                [aa, bb, slope_coarse, slope_fine, slope_total] = gen_values_per_case(valid_cases1[ii], year, nDays)
                qv_vec.append(qvert)
                Cslope_vec.append(abs(slope_coarse))
                Fslope_vec.append(abs(slope_fine))
                Tslope_vec.append(abs(slope_total))
                #print(nBreak)
                #print(qvert)
                break
    
    qVVV.append(qv_vec)
    nbC.append(Cslope_vec)
    nbF.append(Fslope_vec)
    nbT.append(Tslope_vec)
    print(nBreak)
    print(qv_vec)
    print(Cslope_vec)
    plot_image(qv_vec, Cslope_vec, qv_vec, Fslope_vec, 'qvert', 'dC/dt', 'Concentration slope versus qvert for nBreak = ' + str(nBreak), 'C vs qv_nB_' + str(round(nBreak,0)), year, outputFolder2 )

    ax1.plot(qv_vec, Cslope_vec, '*', label='Coarse, nB: ' + str(nBreak_compare) )
    ax1.plot(qv_vec, Fslope_vec, 'o', label='Fines, nB: ' + str(nBreak_compare) )
    ax1.plot(qv_vec, Tslope_vec, '.', label='Total, nB: ' + str(nBreak_compare) )

plt.legend()    
plt.savefig(outputFolder2+"nB_impact_coarse_fine.png", format='png')
plt.close()

#EXTRA PLOT *************
# #nB Impact by Type
# #C
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("qvert")
# ax1.set_ylabel("dC/dt")
# plt.title("nB Impact Coarse")

# for k in range(len(nBreakRange)):
#      ax1.plot(qVVV[k], nbC[k], '*-', label='Coarse, nB: ' + str(nBreakRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"nB_impact_coarse.png", format='png')
# plt.close()

# #F
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("qvert")
# ax1.set_ylabel("dC/dt")
# plt.title("nB Impact Fine")

# for k in range(len(nBreakRange)):
#      ax1.plot(qVVV[k], nbF[k], 'o--', label='Fine, nB: ' + str(nBreakRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"nB_impact_fine.png", format='png')
# plt.close()

# #T
# fig, ax1 = plt.subplots(figsize=(10,10))
# ax1.set_xlabel("qvert")
# ax1.set_ylabel("dC/dt")
# plt.title("nB Impact Total")

# for k in range(len(nBreakRange)):
#      ax1.plot(qVVV[k], nbT[k], '.-', label='Total, nB: ' + str(nBreakRange[k]) )
# plt.legend()    
# plt.savefig(outputFolder2+"nB_impact_total.png", format='png')
# plt.close()
#EXTRA PLOT *************

#MULTID
print('VEC DATA::: ', Vec_MultiD) 
#exit(1)

#TODO
#Contour MultiD Plot for Vec_MultiD

x = np.zeros(( len(qvertRange) * len(nBreakRange) ))
x = np.zeros(( len(Vec_MultiD) ))
y = np.zeros((len(x)))
zSlopeC = np.zeros((len(x)))
zSlopeF = np.zeros((len(x)))
zCases = np.zeros((len(x)))
zErrorC = np.zeros((len(x)))
zErrorF = np.zeros((len(x)))
zSlopeT = np.zeros((len(x)))
zErrorT = np.zeros((len(x)))
zRatioT = np.zeros((len(x)))

ratioV = np.zeros((len(x)))
ratioIntV = np.zeros((len(x)))
cum_mass_lossV = np.zeros((len(x)))
slope_massV = np.zeros((len(x)))

print('VEC DATA 0 ele::: ', Vec_MultiD[0][0]) #CHECK
print("Len A: ", len(qvertRange) * len(nBreakRange)) #CHECK
print("Len B: ", len(Vec_MultiD)) #CHECK
print(refT)


#Convert qB and nB to timescales (units of time)
#For qv timescale
drhoi = 910
dLf = 330000
if year == 2018:
    deltah = 0.5
    deltaT = 1.5
elif year == 2020:
    deltah = 0.3
    deltaT = 2.0
#For nB Timescale
if year == 2018:
    dNf = 493
elif year  == 2020:
    dNf = 300 #Correct!!!!
dDmean = 20
dD0 = 2
Ns = math.log(dDmean**2/dD0**2)/math.log(2)
if year == 2018:
    nDays = 48
    nDaysxy = 1  #For not adimensional days
elif year == 2020:
    nDays = 24 #Confirm!!!
    nDaysxy = 1  #For not adimensional days

#xmin = [((drhoi*dLf*deltah) / (deltaT*25))/86400]
#ymin = [Ns * dNf * (1000)/86400]
xmin = [((drhoi*dLf*deltah) / (deltaT*25))/86400/nDaysxy]
ymin = [Ns * dNf * (1000)/86400/nDaysxy]
xmin = 25 
#Compromise for FSD and concentration
ymin = 86400 * 1/1500
#Best only for concentration
#ymin = 86400 * 1/1000

xvals = []
yvals = []
xvalsr = []
yvalsr = []

qvertRangeF = qvertRange
nBreakRangeF = nBreakRange

#for i in range(len(qvertRange)):  
    #Reverse
    #xvals.append( [((drhoi*dLf*deltah) / (deltaT*qvertRangeF[i]))/86400/nDaysxy] )
    #xvalsr.append( [((drhoi*dLf*deltah) / (deltaT*qvertRangeF[len(qvertRange)-1-i]))/86400/nDaysxy] )
#for j in range(len(nBreakRange)):
    #yvals.append( [Ns * dNf * (nBreakRangeF[i])/86400/nDaysxy] )
    #yvalsr.append( [Ns * dNf * (nBreakRangeF[len(nBreakRange)-1-i])/86400/nDaysxy] )


[c0C, c0F, c0T] = gen_init_conc(valid_cases1[0], year, nDays)
for i in range(len(x)): #CHECK
    xTime = ((drhoi*dLf*deltah) / (deltaT*Vec_MultiD[i][0]))/86400 #For qv in days
    yTime = Ns * dNf * (1/Vec_MultiD[i][1])/86400   #B was inverted here, so we have to return it back
    
    #Using qv and nB
    x[i] = (Vec_MultiD[i][0])
    y[i] = (Vec_MultiD[i][1])
    #Using days
    #x[i] = xTime/nDaysxy
    #y[i] = yTime/nDaysxy
    
    
    zRatioT[i] = x[i] / y[i]
    
    #Proportional to Cbar
    #zSlopeC[i] = 2*(Vec_MultiD[i][2])/(nDays*nDays) #(%/d Units) //Average Approximation
    #zSlopeF[i] = 2*(Vec_MultiD[i][3])/(nDays*nDays) #(%/d Units) //Average Approximation
    #zSlopeT[i] = 2*(Vec_MultiD[i][9])/(nDays*nDays) #(%/d Units) //Average Approximation
    
    cn = 0 #Approximate total loss
    #Inverse Proportional Triangle
    zSlopeC[i] = ( (c0C-cn)**2 ) / ( 2*(Vec_MultiD[i][2]/nDays)*nDays ) #(%/d Units) //Triangle Approximation
    zSlopeF[i] = ( (c0F-cn)**2 ) / ( 2*(Vec_MultiD[i][3]/nDays)*nDays ) #(%/d Units) //Triangle Approximation
    zSlopeT[i] = ( (c0T-cn)**2 ) / ( 2*(Vec_MultiD[i][9]/nDays)*nDays ) #(%/d Units) //Triangle Approximation
    
    zCases[i] = (Vec_MultiD[i][4])
    zErrorC[i] = (Vec_MultiD[i][5])
    zErrorF[i] = (Vec_MultiD[i][6])
    zErrorT[i] = (Vec_MultiD[i][7])
    ratioV[i] = (Vec_MultiD[i][10])
    ratioIntV[i] = (Vec_MultiD[i][11])
    cum_mass_lossV[i] = (Vec_MultiD[i][12]) #Cumulative
    slope_massV[i] = (Vec_MultiD[i][13]) #Mass slope
    

#Dimension of grid for plotting
xgrid = len(qvertRange)
ygrid = len(nBreakRange)
print('xlen: ', xgrid)
print('ylne: ', ygrid)

#Multiple
#Contour Plots
ct = 0
xArr = np.zeros((xgrid,ygrid))
yArr = np.zeros((xgrid,ygrid))
zSlopeCArr = np.zeros((xgrid,ygrid))
zSlopeFArr = np.zeros((xgrid,ygrid))
zErrorCArr = np.zeros((xgrid,ygrid))
zErrorFArr = np.zeros((xgrid,ygrid))
zSlopeTArr = np.zeros((xgrid,ygrid))
zErrorTArr = np.zeros((xgrid,ygrid))
zRatioTArr = np.zeros((xgrid,ygrid))
ratioVArr = np.zeros((xgrid,ygrid))
ratioIntVArr = np.zeros((xgrid,ygrid))
cum_mass_lossVArr = np.zeros((xgrid,ygrid))
slope_massVArr = np.zeros((xgrid,ygrid))

#Save text (results)
dimtx = [len(qvertRange), len(nBreakRange)]

np.savetxt(outputFolder2+'dim_data.dat', dimtx , fmt='%d')
np.savetxt(outputFolder2+'x_data.dat', x , fmt='%.8f')
np.savetxt(outputFolder2+'y_data.dat', y , fmt='%.8f')
np.savetxt(outputFolder2+'zSlopeC_data.dat', zSlopeC , fmt='%.8f')
np.savetxt(outputFolder2+'zSlopeF_data.dat', zSlopeF , fmt='%.8f')
np.savetxt(outputFolder2+'zSlopeT_data.dat', zSlopeT , fmt='%.8f')
np.savetxt(outputFolder2+'zErrorC_data.dat', zErrorC , fmt='%.8f')
np.savetxt(outputFolder2+'zErrorF_data.dat', zErrorF , fmt='%.8f')
np.savetxt(outputFolder2+'zErrorT_data.dat', zErrorT , fmt='%.8f')
np.savetxt(outputFolder2+'ratio_data.dat', ratioV , fmt='%.8f')
np.savetxt(outputFolder2+'ratioInt_data.dat', ratioIntV , fmt='%.8f')
np.savetxt(outputFolder2+'cum_mass_loss_data.dat', cum_mass_lossV , fmt='%.8f')
np.savetxt(outputFolder2+'slope_mass_data.dat', slope_massV , fmt='%.8f')

for i in range(xgrid):
    for j in range(ygrid):
        print(i)
        #Select from MultiV correctly
        xArr[i][j] = x[ct]
        yArr[i][j] = y[ct]
        zSlopeCArr[i][j] = zSlopeC[ct]
        zSlopeFArr[i][j] = zSlopeF[ct]
        zErrorCArr[i][j] = zErrorC[ct]
        zErrorFArr[i][j] = zErrorF[ct]
        zSlopeTArr[i][j] = zSlopeT[ct]
        zErrorTArr[i][j] = zErrorT[ct]
        zRatioTArr[i][j] = zRatioT[ct]
        ratioVArr[i][j] = ratioV[ct]
        ratioIntVArr[i][j] = ratioIntV[ct]
        cum_mass_lossVArr[i][j] = cum_mass_lossV[ct]
        slope_massVArr[i][j] = slope_massV[ct]
        ct += 1

print(xArr)
print(yArr)
print(zSlopeCArr)


#For rounding colorbars

def myfmt(x, pos):
    return '{0:.2f}'.format(x)

def myfmt2(x, pos):
    return '{0:.2f}'.format(x)

#CONTOUR PLOTS
#fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(11,10))

#fig, ax = plt.subplots()
#sc = ax.scatter(x,y, c=y)

#plt.title('Coarse Grain Concentration Inv. Area', size=20)
#plt.title('Concentration Loss Rate: Coarse', size=25)
plt.title('Mean Concentration Rate: Coarse', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, zSlopeCArr, 30)

#CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
#contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
#plt.clabel(contours2, colors = ['black'], inline=True, fontsize=12)

contoursb = plt.contour(xArr, yArr, zSlopeCArr, 20, colors = 'black', linestyles='-')
#contours = plt.contourf(xArr, yArr, zSlopeCArr, 20, cmap='coolwarm')
contours = plt.contourf(xArr, yArr, zSlopeCArr, 20, cmap='coolwarm')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
else: 
    plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt))

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel(r'$\left(\% d^{-1}\right)$', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)

plt.savefig(outputFolder2 + "Results_Plot_AC_ContourR.png")
plt.close()


#Coarse LOG for breakage
fig = plt.figure(figsize=(10,10))
#plt.title('Coarse Grain Concentration Inv. Area', size=20)
#plt.title('Concentration Loss Rate: Coarse LOG', size=25)
plt.title('Mean Concentration: Coarse LOG', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Closest fit to Sat. Observations')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, zSlopeCArr, 30)

#CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
#contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
#plt.clabel(contours2, colors = ['black'], inline=True, fontsize=12)

contoursb = plt.contour(xArr, yArr, zSlopeCArr, 20, colors = 'black', linestyles='-')
contours = plt.contourf(xArr, yArr, zSlopeCArr, 20, cmap='coolwarm_r')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel('Mean Concentration', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=18)
plt.yscale("log")  
if qv_val:
    plt.xlabel('Melt Rate - qv (W/m2/°C)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency - B (1/d)')
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_AC_ContourLOG.png")
plt.close()


#fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(11,10))
#plt.title('Fine Concentration Inv. Area', size=20)
plt.title('Mean Concentration Rate: Fine', size=25)
#plt.title('Concentration Loss Rate: Fine', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
plt.legend(fontsize = 18)

# CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
# plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
# #contours2 = plt.contour(xArr, yArr, zSlopeFArr, 30)
# contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
# #plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
# plt.clabel(contours2, colors = ['k'], inline=True, fontsize=12)

contoursb = plt.contour(xArr, yArr, zSlopeFArr, 19, colors = 'black', linestyles='-')
contours = plt.contourf(xArr, yArr, zSlopeFArr, 19, cmap='coolwarm')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
else: 
    plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt))

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel(r'$\left(\% d^{-1}\right)$', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "Results_Plot_AF_ContourR.png")
plt.close()

#TOTAL
#fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(11,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Mean Concentration Rate: Total', size=25)
#plt.title('Concentration Loss Rate: Total', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
plt.legend(fontsize = 18)

# CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
# plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
# #contours2 = plt.contour(xArr, yArr, zSlopeTArr, 30)
# contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
# #plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
# plt.clabel(contours2, colors = ['k'], inline=True, fontsize=12)

contoursb = plt.contour(xArr, yArr, zSlopeTArr, 20, colors = 'black', linestyles='-')
contours = plt.contourf(xArr, yArr, zSlopeTArr, 20, cmap='coolwarm')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
else: 
    plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt))

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel(r'$\left(\% d^{-1}\right)$', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "Results_Plot_AT_ContourR.png")
plt.close()



# #Analytical expression fit

# # The function to be fit is Z. (DataSet)
# Za = zSlopeTArr

# #define real function to create data from #FUNCTION TO FIT with parameters

# def f(x, y, a, b, c):
#     g = b*(y)*np.exp(a*x)+c
#     return g

# # #Add analytical expression
# # def f(x, y):
# #     alpha = 0.02482;
# #     beta = 0.0064;
# #     gamma = 18.9;
# #     return beta*(y)*np.exp(alpha*x)+gamma

# # Create x and y indices
# # xAN = np.zeros(len(xvals))
# # yAN = np.zeros(len(yvals))
# # for i in range(len(xvals)):
# #     qvertRangeV = qvertRangeF[i]
# #     xAN[i] = [((drhoi*dLf*deltah) / (deltaT*qvertRangeV))/86400/nDaysxy]
# # for i in range(len(yvals)):
# #     nBreakRangeV = nBreakRangeF[i]
# #     yAN[i] = [Ns * dNf * (nBreakRangeV)/86400/nDaysxy]
# #X1, Y1 = np.meshgrid(xvals, yvals)
# #X1r, Y1r = np.meshgrid(xvalsr, yvalsr)

# #Argument to fit
# def fit_f(M, *args):     #ARGUMENT TO FIT with curve fitting package
#     x, y = M
#     arr = np.zeros(x.shape)
#     arr = f(x, y, *args)
#     return arr

# #x_data = np.vstack((X1.ravel(), Y1.ravel()))
# x_data = np.vstack((xArr.ravel(), yArr.ravel())) #coordinates for contours calc.

# flat_data = Za.ravel()

# #initial guess values
# #p0 = [0.01482, 0.01064, 18.9] #For days
# p0 = [1.482, 1.064, 1.9]

# #FIT
# popt, pcov = curve_fit(fit_f, x_data, flat_data, p0)
# fit = f(xArr, yArr, *popt)
# fitp = fit.reshape(len(qvertRange),len(nBreakRange))


# print('Fitted parameters:')
# print(*popt)

# rms = np.sqrt(np.mean((Za - fitp)**2))
# print('RMS residual =', rms)
# print('Data:')
# print(Za)
# print('Fit result:')
# print(fitp)
# # print('DataX:')
# # print(xArr)
# # print('DataY:')
# # print(xArr)
# # print('FitX:')
# # print(X1)
# # print('FitY:')
# # print(Y1)

# ####FIT TRY



# fig = plt.figure(figsize=(10,10))
# #plt.title('Total Concentration Inv. Area', size=20)
# #plt.rcParams['text.usetex'] = True
# plt.title('Total Floes and Approximation: ' +  str(round(popt[1],3)) + '*Tb'+ ' exp(' + str(round(popt[0],3)) + '*Tm) + ' + str(round(popt[2],3)) , size=18)
# #plt.plot(xmin,ymin,'k*',markersize=20)

# # CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
# # plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
# # #contours2 = plt.contour(xArr, yArr, zSlopeTArr, 30)
# # contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
# # #plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
# # plt.clabel(contours2, colors = ['k'], inline=True, fontsize=12)

# contoursb = plt.contour(xArr, yArr, zSlopeTArr, 15, linestyles='-', colors = 'black') #colors = 'black',  cmap='RdGy'
# plt.clabel(contoursb, fmt = '%2.4f', colors = 'k', fontsize=14, inline=True) 
# contours = plt.contourf(xArr, yArr, zSlopeTArr, 15, cmap='coolwarm')  #Reversed coolwarm_r
# #plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
# #plt.clabel(contours, inline=True, fontsize=8, colors='k')
# cbar = plt.colorbar();
# cbar.ax.get_yaxis().labelpad = 15
# #cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
# cbar.ax.set_ylabel('Cumulative Concentration', rotation=270, labelpad=15, size=15)
# cbar.ax.tick_params(labelsize=15)
# #plt.gca().set_aspect('equal', adjustable='box')
# #plt.gca().set_aspect('equal', adjustable='box')
# axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
# axes.tick_params(labelsize=15)
# CAN = plt.contour(xArr, yArr, fitp, 15, linestyles='--', colors='blue');  #colors='blue',   cmap='RdGy'
# plt.clabel(CAN, fmt = '%2.4f', colors = 'b', fontsize=14, inline=True)
# if qv_val:
#     plt.xlabel('qvert (W/m2/Celsius)')    
# else: 
#     plt.xlabel('kHor (m2/s)')
# plt.ylabel('Break Frequency (1/s)')
# # plt.xlabel('Melt Time (days)')
# # plt.ylabel('Break Time (days)')
# # plt.xlabel('Melt Nondimensional Time')
# # plt.ylabel('Break Nondimensional Time')
# #plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
# plt.savefig(outputFolder2 + "Results_Plot_AT_ContourAN.png")
# plt.close()


#Ratio bkg/melt contour
#fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(10,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Breakage vs. Melt Comparison', size=20)
plt.plot(xmin,ymin,'k*',markersize=20, label='Obs. Fit')
contours2 = plt.contour(xArr, yArr, ratioVArr, 30, colors = 'black')
CS = plt.contour(xArr, yArr, ratioVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=14)
contours = plt.contourf(xArr, yArr, ratioVArr, 15, cmap='coolwarm')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel('Melt Rate - qv (W/m2/°C)', fontsize=20)    
else: 
    plt.xlabel('kHor (m2/s)', fontsize=20) 
plt.ylabel('Break Frequency - B (1/d)', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax)

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel('Breakage/Melt Ratio', rotation=270, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "Results_Plot_Ratio_Contour.png")
plt.close()

#Ratio bkg/melt contour and Cum Mass Presence
fig = plt.figure(figsize=(12,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Breakage vs. Melt Ratio', size=25)
#plt.plot(xmin,ymin,'k*',markersize=20)
#contours2 = plt.contour(xArr, yArr, ratioIntVArr, 15, colors = 'black')
CS = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = ['k'], fontsize=12)
CS2 = plt.contour(xArr, yArr, ratioIntVArr, 20, colors=('k',), linestyles=('-',), linewidths=(1.5,))
plt.clabel(CS2, fmt = '%2.2f', colors = ['k'], fontsize=12)

#contours3 = plt.contour(xArr, yArr, zRatioTArr, 30)
#plt.clabel(contours3, colors = 'black', inline=True, fontsize=12)

#contours = plt.contourf(xArr, yArr, ratioIntVArr, 15, cmap='coolwarm') #Area bkg/melt ratio
contours = plt.contourf(xArr, yArr, cum_mass_lossVArr, 15, cmap='coolwarm') #Mass loss
#contour3 = plt.contour(xArr, yArr, cum_mass_lossVArr, 15, colors = 'green') #Mass loss
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 20
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
#cbar.ax.set_ylabel('Breakage/Melt Ratio (Sum Area)', rotation=270, labelpad=20)
cbar.ax.set_ylabel('Cumulative Mass Loss (Normalized)', rotation=270, labelpad=20, size = 18)
cbar.ax.tick_params(labelsize=20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=18)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
#plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_Ratio_SumArea_Contour.png")
plt.close()

#Cumulative Mass Presence contour (or gradient???)
fig = plt.figure(figsize=(10,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Coarse Floe Average Mass Loss', size=24)
plt.plot(xmin,ymin,'k*',markersize=20)

# CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('--',),linewidths=(2,))
# plt.clabel(CS, fmt = '%2.1d', colors = 'black', fontsize=15)
# contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black', linestyles='--')
# #contours2 = plt.contour(xArr, yArr, cum_mass_lossVArr, 30)
# #plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
# plt.clabel(contours2, colors = ['k'], inline=True, fontsize=12)

contoursb = plt.contour(xArr, yArr, cum_mass_lossVArr, 25, colors = 'black', linestyles='-')
contours = plt.contourf(xArr, yArr, cum_mass_lossVArr, 25, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel('Nondimensional Mass Loss', rotation=270, labelpad=15, size=15)
cbar.ax.tick_params(labelsize=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_CumMassLoss_Contour.png")
plt.close()


#Ratio bkg/melt contour and Mass Slope
fig = plt.figure(figsize=(12,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Breakage vs. Melt Ratio and Mass Slope', size=25)
#plt.plot(xmin,ymin,'k*',markersize=20)
#contours2 = plt.contour(xArr, yArr, ratioIntVArr, 15, colors = 'black')
CS = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = ['k'], fontsize=12)
CS2 = plt.contour(xArr, yArr, ratioIntVArr, 20, colors=('k',), linestyles=('-',), linewidths=(1.5,))
plt.clabel(CS2, fmt = '%2.2f', colors = ['k'], fontsize=12)

#contours3 = plt.contour(xArr, yArr, zRatioTArr, 30)
#plt.clabel(contours3, colors = 'black', inline=True, fontsize=12)

#contours = plt.contourf(xArr, yArr, ratioIntVArr, 15, cmap='coolwarm') #Area bkg/melt ratio
contours = plt.contourf(xArr, yArr, slope_massVArr, 15, cmap='coolwarm') #Mass loss
#contour3 = plt.contour(xArr, yArr, cum_mass_lossVArr, 15, colors = 'green') #Mass loss
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 20
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
#cbar.ax.set_ylabel('Breakage/Melt Ratio (Sum Area)', rotation=270, labelpad=20)
cbar.ax.set_ylabel('Slope Mass Loss (kg/day)', rotation=270, labelpad=20, size = 18)
cbar.ax.tick_params(labelsize=20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=18)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
#plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_Ratio_SlopeMass_Contour.png")
plt.close()


#Mass Slope only
fig = plt.figure(figsize=(11,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Mass Loss Rate: Coarse', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Closest fit to Sat. Observations')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, ratioIntVArr, 15, colors = 'black')
#CS = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = ['k'], fontsize=12)
CS2 = plt.contour(xArr, yArr, slope_massVArr, 20, colors=('k',), linestyles=('-',), linewidths=(1.5,))
#plt.clabel(CS2, fmt = '%2.2f', colors = ['k'], fontsize=12)

#contours3 = plt.contour(xArr, yArr, zRatioTArr, 30)
#plt.clabel(contours3, colors = 'black', inline=True, fontsize=12)

#contours = plt.contourf(xArr, yArr, ratioIntVArr, 15, cmap='coolwarm') #Area bkg/melt ratio
contours = plt.contourf(xArr, yArr, slope_massVArr, 20, cmap='coolwarm') #Mass loss
#contour3 = plt.contour(xArr, yArr, cum_mass_lossVArr, 15, colors = 'green') #Mass loss
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 20
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
#cbar.ax.set_ylabel('Breakage/Melt Ratio (Sum Area)', rotation=270, labelpad=20)
cbar.ax.set_ylabel('Mass Loss Rate (kg/day)', rotation=270, labelpad=20, size = 18)
cbar.ax.tick_params(labelsize=20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=18)
if qv_val:
    plt.xlabel('qv (W/m2/°C)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('B (1/d)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
#plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_Ratio_SlopeMass_ONLY.png")
plt.close()

#Mass Slope only LOG
fig = plt.figure(figsize=(11,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Mass Loss Rate LOG', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Closest fit to Sat. Observations')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, ratioIntVArr, 15, colors = 'black')
#CS = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = ['k'], fontsize=12)
CS2 = plt.contour(xArr, yArr, slope_massVArr, 20, colors=('k',), linestyles=('-',), linewidths=(1.5,))
#plt.clabel(CS2, fmt = '%2.2f', colors = ['k'], fontsize=12)

#contours3 = plt.contour(xArr, yArr, zRatioTArr, 30)
#plt.clabel(contours3, colors = 'black', inline=True, fontsize=12)

#contours = plt.contourf(xArr, yArr, ratioIntVArr, 15, cmap='coolwarm') #Area bkg/melt ratio
contours = plt.contourf(xArr, yArr, slope_massVArr, 20, cmap='coolwarm') #Mass loss
#contour3 = plt.contour(xArr, yArr, cum_mass_lossVArr, 15, colors = 'green') #Mass loss
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 20
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
#cbar.ax.set_ylabel('Breakage/Melt Ratio (Sum Area)', rotation=270, labelpad=20)
cbar.ax.set_ylabel('Mass Loss Rate (kg/day)', rotation=270, labelpad=20, size = 18)
cbar.ax.tick_params(labelsize=20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=18)
plt.yscale("log") 
if qv_val:
    plt.xlabel('qv (W/m2/°C)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('B (1/d)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
#plt.savefig(outputFolder2 + "Results_Plot_InvAT_Contour.png")
plt.savefig(outputFolder2 + "Results_Plot_Ratio_SlopeMass_ONLY_LOG.png")
plt.close()

#Ratio bkg/melt contour
#fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(11,10))
#plt.title('Total Concentration Inv. Area', size=20)
plt.title('Breakage vs. Melt Dominance: Coarse', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, ratioIntVArr, 15, colors = 'black')
CS = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(3.5,))
plt.clabel(CS, fmt = '%2.2f', colors = ['k'], fontsize=16)
CS2 = plt.contour(xArr, yArr, ratioIntVArr, 30, colors=('k',), linestyles=('-',), linewidths=(1.5,))
plt.clabel(CS2, fmt = '%2.2f', colors = ['k'], fontsize=14)
#CS3 = plt.contour(xArr, yArr, ratioIntVArr,levels = [1.1200000,1.12000001], colors=('k',),linestyles=('--',),linewidths=(2.0,))
#plt.clabel(CS3, fmt = '%2.2f', colors = ['k'], fontsize=14)

#contours3 = plt.contour(xArr, yArr, zRatioTArr, 30)
#plt.clabel(contours3, colors = 'black', inline=True, fontsize=12)

#contours = plt.contourf(xArr, yArr, ratioIntVArr, 15, cmap='coolwarm') #Area bkg/melt ratio
contours = plt.contourf(xArr, yArr, ratioIntVArr, 30, cmap='coolwarm_r') #Mass loss
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
else: 
    plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt2))

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel(r'$\mu$', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)
plt.savefig(outputFolder2 + "Results_Plot_RatioBM_ONLY.png")
plt.close()



#Point Plots

#Coarse
fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = zSlopeC, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Coarse Grain Slope', size=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(zSlopeC):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_Plot_SlopeC_Points.png")
plt.close()

#Fine
fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = zSlopeF, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Fine Slope', size=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(zSlopeF):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_Plot_SlopeF_Points.png")
plt.close()

#Cases
fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = zCases, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Case Code', size=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(zCases):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_Cases_Points.png")
plt.close()


#####
# #Same as above but for Error Contours
fig, ax = plt.subplots(figsize=(11,10))
#plt.title('Coarse Grain Concentration Inv. Area', size=20)
#plt.title('Concentration Loss Rate: Coarse', size=25)
plt.title('Coarse Floe Concentration RMSE', size=25)
plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
plt.legend(fontsize = 18)
#contours2 = plt.contour(xArr, yArr, zSlopeCArr, 30)

#CS = plt.contour(xArr, yArr, zRatioTArr,levels = [1.000000,1.0000001], colors=('k',),linestyles=('-',),linewidths=(2,))
#plt.clabel(CS, fmt = '%2.1d', colors = 'k', fontsize=15)
#contours2 = plt.contour(xArr, yArr, zRatioTArr, 30, colors = 'black')
#plt.clabel(contours2, colors = ['black'], inline=True, fontsize=12)

contours2 = plt.contour(xArr, yArr, zErrorCArr, 35, colors=('k',), linestyles=('-',), linewidths=(1.0,))
contours = plt.contourf(xArr, yArr, zErrorCArr, 35, cmap='coolwarm')
#plt.clabel(contours2, contours2.levels, inline=True, fontsize=15)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
#axes = plt.gca()
# axes.xaxis.label.set_size(20)
# axes.yaxis.label.set_size(20)
#axes.tick_params(labelsize=18)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
if qv_val:
    plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
else: 
    plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# plt.xlabel('Melt Time (days)')
# plt.ylabel('Break Time (days)')
# plt.xlabel('Melt Nondimensional Time')
# plt.ylabel('Break Nondimensional Time')
#plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
#cbar = plt.colorbar();
cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt))

#cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
cbar.ax.set_ylabel('RMSE', rotation=90, labelpad=18, size=18)
#cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
cbar.ax.tick_params(labelsize=18)

plt.savefig(outputFolder2 + "Error_C_Contour.png")
plt.close()

# fig = plt.figure(figsize=(12,10))
# plt.title('Coarse Floe Concentration RMSE', size=25)
# plt.plot(xmin,ymin,'k*',markersize=18,label='Obs. Fit')
# contours2 = plt.contour(xArr, yArr, zErrorCArr, 35, colors=('k',), linestyles=('-',), linewidths=(1.5,))
# contours = plt.contourf(xArr, yArr, zErrorCArr, 35, cmap='coolwarm')
# #plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
# #plt.clabel(contours, inline=True, fontsize=8, colors='k')
# ax.xaxis.set_tick_params(labelsize=18)
# ax.yaxis.set_tick_params(labelsize=18)
# if qv_val:
#     plt.xlabel(r'Melt Rate - $q_{v} \left(W m^{-2} °C^{-1}\right)$', fontsize=20)    
# else: 
#     plt.xlabel(r'kHor ($\frac{m^{2}}{s}$)', fontsize=20) 
# plt.ylabel(r'Break Frequency - $B \left(d^{-1}\right)$', fontsize=20) 
# # plt.xlabel('Melt Time (days)')
# # plt.ylabel('Break Time (days)')
# # plt.xlabel('Melt Nondimensional Time')
# # plt.ylabel('Break Nondimensional Time')
# #plt.savefig(outputFolder2 + "Results_Plot_InvAC_Contour.png")
# #cbar = plt.colorbar();
# cax = make_square_axes_with_colorbar(ax, size=0.23, pad=0.15)
# cbar = fig.colorbar(contours, cax=cax, format=matplotlib.ticker.FuncFormatter(myfmt2))

# #cbar.ax.get_yaxis().labelpad = 15
# #cbar.ax.set_ylabel('Concentration Inv. Area (1/d)', rotation=270, labelpad=15)
# cbar.ax.set_ylabel('RMSE', rotation=90, labelpad=18, size=18)
# #cbar.ax.set_ylabel('Concentration Loss Rate (%/d)', rotation=270, labelpad=18, size=18)
# cbar.ax.tick_params(labelsize=18)
# plt.savefig(outputFolder2 + "Error_C_Contour.png")
# plt.close()

fig = plt.figure(figsize=(10,10))
plt.title('Fine Error Concentration', size=20)
contours = plt.contourf(xArr, yArr, zErrorFArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
#cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Error', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
plt.savefig(outputFolder2 + "Error_F_Contour.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.title('Total Error Concentration', size=20)
contours = plt.contourf(xArr, yArr, zErrorTArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
#cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Error', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
plt.savefig(outputFolder2 + "Error_T_Contour.png")
plt.close()

#Point Plots

#Coarse
fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = zErrorC, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Coarse Grain Error', size=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(zErrorC):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2+ "Error_C_Points.png")
plt.close()

#Fine
fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = zErrorF, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Fine Error', size=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(zErrorF):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Error_F_Points.png")
plt.close()


fig11 = plt.figure(figsize=(10,10))
plt.title('Coarse Grain Concentration Slope', size=20)
contours = plt.contourf(xArr, yArr, zSlopeCArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
#cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Concentration Slope (1/d)', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')


fig12 = plt.figure(figsize=(10,10))
plt.title('Fine Concentration Slope', size=20)
contours = plt.contourf(xArr, yArr, zSlopeFArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Area (1/d)', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
if qv_val:
    plt.xlabel('qvert (W/m2/Celsius)')    
else: 
    plt.xlabel('kHor (m2/s)')
plt.ylabel('Break Frequency (1/s)')

#show_image_list(list_images=[fig11, fig12], num_cols=2, figsize=(20, 10), grid=False, title_fontsize=20)


plt.savefig(outputFolder2 + "Results_Plot_SlopeCF_Contour_Join.png")
plt.close()

print("Initial concentrations C/F/T: ", c0C, " ", c0F, " ", c0T)

##Get velocity data
vel_data = []

for k in range(caseL):
    ave_vel = read_ave_velocity(start_case+k, year)
    [nBreak, qvert, Qatm] = read_params(start_case+k, year)
    vel_data.append([start_case+k, ave_vel, nBreak])  #Output (Case Number, Velocity, Breakage Step)

np.savetxt(outputFolder2 + 'velocity_ave_Bkg.dat', vel_data, fmt='%d %.8f %.8f')
