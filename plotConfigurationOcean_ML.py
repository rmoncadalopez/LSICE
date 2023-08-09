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

#Basic Plot Functions
def plot_simple(y_label, time, dataV, data_label, plot_title, File_Name, mainVfolder):
    fig = plt.figure(figsize=(10,10))
    plt.xlabel('Time (days)', fontsize = 20)
    plt.ylabel(y_label, fontsize = 20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
    plt.legend(prop={'size': 20})
    plt.title(plot_title, fontsize = 25)
    plt.savefig(mainVfolder + File_Name, format='png')
    plt.close()
    
#Basic Plot Functions (x2)
def plot_double(y_label, time, dataV, data_label, dataV2, data_label2, plot_title, File_Name, mainVfolder):
    fig = plt.figure(figsize=(10,10))
    plt.xlabel('Time (days)', fontsize = 20)
    plt.ylabel(y_label, fontsize = 20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
    plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
    plt.legend(prop={'size': 20})
    plt.title(plot_title, fontsize = 25)
    plt.savefig(mainVfolder + File_Name, format='png')
    plt.close()    

#(x3)
def plot_double(y_label, time, dataV, data_label, dataV2, data_label2, dataV3, data_label3, plot_title, File_Name, mainVfolder):
    fig = plt.figure(figsize=(10,10))
    plt.xlabel('Time (days)', fontsize = 20)
    plt.ylabel(y_label, fontsize = 20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
    plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
    plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
    plt.legend(prop={'size': 20})
    plt.title(plot_title, fontsize = 25)
    plt.savefig(mainVfolder + File_Name, format='png')
    plt.close()  


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

if create_folder:
    caseNo = sys.argv[1]
    folder_title = sys.argv[2]
    mainVfolder = main_dir + "MassLoss_" + str(folder_title) + "/"
    #mainVfolder = "./Visualizaton_" + str(folder_title) + "/"
    if year == 2020:
        mainOutputfolder = "./Output/SeaIce2020_" + str(caseNo) + "/" #2020
    elif year == 2018:
        #mainOutputfolder = "./Output/SeaIce2018Prob_" + str(caseNo) + "/" #2018
        mainOutputfolder = "./Output/SeaIce2018_" + str(caseNo) + "_22/" #2018
        #mainOutputfolder = "./Output/SeaIce_" + str(caseNo) + "/" #2018
        #mainOutputfolder = "./Output/SeaIce2018h_" + str(caseNo) + "/" #2018
    
    #Check if folder already exists
    if os.path.exists(mainVfolder) == False:
        #Create folders
        os.mkdir(mainVfolder)
else:
    mainVfolder = main_dir    
    #mainVfolder = "./Visualization/"
    mainOutputfolder = "./Output/SeaIce/"

#Get Mass Loss File from Main_Periodic_Ocean
massLossFile = mainOutputfolder + "massLoss0.dat"
massLossV = np.loadtxt(massLossFile)

#Get number of steps
numbergFile        =  mainOutputfolder+"numberg0.dat"
numberg            = np.loadtxt(numbergFile)
nSteps             = int(len(numberg))
print(nSteps)


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
loss_mcv_solar = np.zeros(nSteps) # Loss Melt Coarse Vertical (Thickness) due to sun
loss_mcv_ocean = np.zeros(nSteps) # Loss Melt Coarse Vertical (Thickness) due to ocean

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
    loss_mcv_solar[i]          = massLossV[i][11]
    loss_mcv_ocean[i]          = massLossV[i][12]

loss_Total = np.zeros(nSteps)
loss_Tfines = np.zeros(nSteps)
loss_Tcoarse = np.zeros(nSteps)
loss_TotalA = np.zeros(nSteps)
loss_TcoarseA = np.zeros(nSteps)
for i in range(nSteps):
    if i == 0:
        loss_Total[i] = 0
        loss_Tfines[i] = 0
        loss_Tcoarse[i] = 0
    else:
        loss_Total[i] = - tmass[i] + tmass[i-1]  
        loss_Tcoarse[i] = - tmasscoarse[i] + tmasscoarse[i-1] 
        loss_Tfines[i] = - tmassfines[i] + tmassfines[i-1]  
        loss_TotalA[i] = - tmassA[i] + tmassA[i-1]  
        loss_TcoarseA[i] = - tmasscoarseA[i] + tmasscoarseA[i-1] 
        
#Plot total_melt_coarse
plot_simple('Ice Mass Loss (kg)', time, total_melt_coarse, 'Cumulative Total_melt_coarse', 'Cumulative Total_melt_coarse',  'Deriv_Melt_Coarse.png', mainVfolder)
#Plot Melt Coarse Lateral
plot_simple('Ice Mass Loss (kg)', time, loss_mcl, 'Cumulative Coarse Breakage Loss', 'Cumulative Coarse Breakage Loss',  'Coarse_Bkg_Loss.png', mainVfolder)
#Plot Melt Coarse Vertical
plot_simple('Ice Mass Loss (kg)', time, loss_mcv, 'Cumulative Coarse Vertical Melt', 'Cumulative Coarse Vertical Melt',  'Coarse_Vertical_Melt.png', mainVfolder)

#Compare all Melt Coarse
y_label = 'Ice Mass Loss (kg)'
data_label = 'Coarse Breakage Loss'            
data_label2 = 'Coarse Vertical Melt' 
data_label3 = 'Total_melt_coarse'
plot_title = 'Evolution of Cumulative Coarse Mass Loss (by Type)'
File_Name = 'MassLossCoarse_AllTypes.png'
dataV = loss_mcl
dataV2 = loss_mcv
dataV3 = total_melt_coarse
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  

#Compare all Melt Coarse (for PPT)
y_label = 'Fraction of Loss'
data_label = 'Coarse Breakage Decay Fraction'            
data_label2 = 'Coarse Vertical Melt Fraction' 
#data_label3 = 'Total_melt_coarse'
plot_title = 'Comparison of Coarse Mass Loss Mechanisms'
File_Name = 'MassLossCoarse_Compare.png'
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

if tstop > 0:
    mass_slope =  total_melt_coarse[-1] / tstop #len(pre_time)
else:
    mass_slope = 0.0
    
dataV = np.divide(pre_mcl, pre_total)[1:] #Make it pretty removing first point
dataV2 = np.divide(pre_mcv, pre_total)[1:]
dataV3 = total_melt_coarse
pre_time = pre_time[1:]
pre_time = np.array(pre_time)
fig = plt.figure(figsize=(10,10))
#plt.xlabel('Nondimensional Time', fontsize = 20)
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(pre_time/(nSteps-1), dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(pre_time/(nSteps-1), dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
#plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  

areaBkg = simps(dataV, dx=1)
areaMelt = simps(dataV2, dx=1)
# tMassC = total_melt_coarse[-1]
# dataV3 = total_melt_coarse/tMassC/(nSteps-1)
tMassC = pre_total[-1]
dataV3 = pre_total/tMassC/(len(pre_total))
ratioInt = areaBkg/areaMelt
areaMassLoss = simps(dataV3, dx=1)

text_exp = []
text_exp.append(dataV[-1]/dataV2[-1])
text_exp.append(ratioInt)
text_exp.append(areaMassLoss)
text_exp.append(mass_slope)
np.savetxt(mainVfolder+"BkMeltratio.dat", text_exp, fmt='%.8f')

#Plot Bkg/Melt Ratio as Bar Plot
fig = plt.figure(figsize=(10,10))
#plt.xlabel('Nondimensional Time', fontsize = 20)
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('Melt or Break Fraction', fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.bar(range(len(dataV2)), dataV2, label='Melt Loss Fraction', alpha=0.99, color='r', width = 1.0)
plt.bar(range(len(dataV)), dataV, bottom=dataV2, label='Breakage Loss Fraction', alpha=0.99, color='b', width = 1.0)
plt.axhline(0.5, linestyle='--', color='k')
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + 'MassLossCoarse_CompareBAR.png', format='png')
plt.close()


#Compare all Melt Coarse Rate
y_label = 'Ice Mass Loss (kg/day)'
data_label = 'Coarse Breakage Loss'            
data_label2 = 'Coarse Vertical Melt' 
data_label3 = 'Total_melt_coarse'
plot_title = 'Evolution of Rate of Coarse Mass Loss (by Type)'
File_Name = 'MassLossCoarseRate_AllTypes.png'
dataV = np.gradient(loss_mcl)
dataV2 = np.gradient(loss_mcv)
dataV3 = np.gradient(total_melt_coarse)
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  

#Melt Ocean vs. Solar
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('Cumulative Mass Loss (kg)', fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plot_title = 'Evolution of Coarse Mass Loss (by Melt Type)'
plt.title(plot_title, fontsize = 25)
File_Name = 'Sep_MassLossCoarseRate_MeltAllTypes.png'
data_label = 'Total Coarse Melt'            
data_label2 = 'Solar Coarse Melt' 
data_label3 = 'Ocean Coarse Melt'
dataV = loss_mcv
dataV2 = loss_mcv_solar
dataV3 = loss_mcv_ocean
plt.plot(time, dataV, 'g', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'b', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  


#Plot Fine Gains
plot_simple('Ice Mass Loss (kg)', time, gain_fines, 'Cumulative Fine Gains', 'Cumulative Fine Gains',  'Fine_Gains.png', mainVfolder)
#Plot Fine Loss
plot_simple('Ice Mass Loss (kg)', time, loss_fines, 'Cumulative Fine Loss', 'Cumulative Fine Loss',  'Fine_Loss.png', mainVfolder)    
#Plot Fine Net Loss
plot_simple('Ice Mass Loss (kg)', time, fine_net_loss, 'Cumulative Fine Net Loss', 'Cumulative Fine Net Loss',  'Fine_Net_Loss.png', mainVfolder)    

#Compare all Melt Fine
y_label = 'Ice Mass Loss (kg)'
data_label = 'Fine Gains'            
data_label2 = 'Fine Loss' 
data_label3 = 'Fine Net Loss'
plot_title = 'Evolution of Cumulative Fine Mass Loss (by Type)'
File_Name = 'MassLossFine_AllTypes.png'
dataV = gain_fines
dataV2 = loss_fines
dataV3 = fine_net_loss
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  

#Compare all Melt Fine
y_label = 'Ice Mass Loss (kg/day)'
data_label = 'Fine Gains'            
data_label2 = 'Fine Loss' 
data_label3 = 'Fine Net Loss'
plot_title = 'Evolution of Rate of Fine Mass Loss (by Type)'
File_Name = 'MassLossFineRate_AllTypes.png'
dataV = np.gradient(gain_fines)
dataV2 = np.gradient(loss_fines)
dataV3 = np.gradient(fine_net_loss)
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close() 

# #Plot total mass versus time _mass
# plot_simple('Total Ice Mass (kg)', time, tmass, 'Total Ice Mass', 'Evolution of Ice Mass',  'TotIceMass_Evo.png', mainVfolder)

#Plot total mass versus time points
plot_simple('Total Ice Mass (kg)', time, tmassA, 'Total Ice Mass', 'Evolution of Ice Mass',  'TotIceMassP_Evo.png', mainVfolder)

# #Compare Total Mass Methods
# y_label = 'Total Ice Mass (kg)'
# data_label = 'Total Ice Mass (_mass)'              
# data_label2 = 'Total Ice Mass (points)'                     
# plot_title = 'Evolution of Ice Mass (Comparison)'
# File_Name = 'TotIceMassAvsP_Evo.png'
# dataV = tmass
# dataV2 = tmassA
# fig = plt.figure(figsize=(10,10))
# plt.xlabel('Time (days)', fontsize = 20)
# plt.ylabel(y_label, fontsize = 20)
# plt.tick_params(axis='x', labelsize=15)
# plt.tick_params(axis='y', labelsize=15)
# plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
# plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
# plt.legend(prop={'size': 20})
# plt.title(plot_title, fontsize = 25)
# plt.savefig(mainVfolder + File_Name, format='png')
# plt.close()   

# #Compare Coarse Mass Methods
# y_label = 'Coarse Ice Mass (kg)'
# data_label = 'Coarse Ice Mass (_mass)'              
# data_label2 = 'Coarse Ice Mass (points)'                    
# plot_title = 'Evolution of Coarse Ice Mass (Comparison)'
# File_Name = 'TotIceMassAvsP_Coarse.png'
# dataV = tmasscoarse
# dataV2 = tmasscoarseA
# fig = plt.figure(figsize=(10,10))
# plt.xlabel('Time (days)', fontsize = 20)
# plt.ylabel(y_label, fontsize = 20)
# plt.tick_params(axis='x', labelsize=15)
# plt.tick_params(axis='y', labelsize=15)
# plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
# plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
# plt.legend(prop={'size': 20})
# plt.title(plot_title, fontsize = 25)
# plt.savefig(mainVfolder + File_Name, format='png')
# plt.close()  

#Plot total mass versus time points
plot_simple('Total Fine Ice Mass (kg)', time, tmassfines, 'Total Fine Ice Mass', 'Evolution of Fine Ice Mass',  'TotIceMassFine_Evo.png', mainVfolder)

#Plot total mass versus time points
plot_simple('Total Coarse Ice Mass (kg)', time, tmasscoarseA, 'Total Ice Mass Coarse Points', 'Evolution Coarse of Ice Mass Points',  'TotIceMassPCoarse_Evo.png', mainVfolder)


# # #Plot total mass versus time (total, coarse, fine) _mass
# # plot_triple('Ice Mass (kg)', time, tmass, 'Total Ice Mass', tmasscoarse, 'Coarse Ice Mass', tmassfines, 'Fine Ice Mass', 'Evolution of Ice Mass (by Type)',  'TotIceMass_Types.png', mainVfolder)
# y_label = 'Ice Mass (kg)'
# data_label = 'Total Ice Mass'            
# data_label2 = 'Coarse Ice Mass' 
# data_label3 = 'Fine Ice Mass'
# plot_title = 'Evolution of Ice Mass (by Type)'
# File_Name = 'TotIceMass_Types.png'
# dataV = tmass
# dataV2 = tmasscoarse
# dataV3 = tmassfines
# fig = plt.figure(figsize=(10,10))
# plt.xlabel('Time (days)', fontsize = 20)
# plt.ylabel(y_label, fontsize = 20)
# plt.tick_params(axis='x', labelsize=15)
# plt.tick_params(axis='y', labelsize=15)
# plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
# plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
# plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
# plt.legend(prop={'size': 20})
# plt.title(plot_title, fontsize = 25)
# plt.savefig(mainVfolder + File_Name, format='png')
# plt.close()  

# #Plot total mass versus time (total, coarse, fine) Points
# plot_triple('Ice Mass (kg)', time, tmassA, 'Total Ice Mass', tmasscoarseA, 'Coarse Ice Mass', tmassfines, 'Fine Ice Mass', 'Evolution of Ice Mass (by Type) using Points',  'TotIceMass_Types.png', mainVfolder)
y_label = 'Ice Mass (kg)'
data_label = 'Total Ice Mass'            
data_label2 = 'Coarse Ice Mass' 
data_label3 = 'Fine Ice Mass'
plot_title = 'Evolution of Ice Mass (by Type) using Points'
File_Name = 'TotIceMass_TypesPts.png'
dataV = tmassA
dataV2 = tmasscoarseA
dataV3 = tmassfines
fig = plt.figure(figsize=(10,10))
plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel(y_label, fontsize = 20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.plot(time, dataV, 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
plt.plot(time, dataV2, 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
plt.plot(time, dataV3, 'g', marker='s', fillstyle='none', label=data_label3) #Same Size Data
plt.legend(prop={'size': 20})
plt.title(plot_title, fontsize = 25)
plt.savefig(mainVfolder + File_Name, format='png')
plt.close()  

#Plot Fine Loss
plot_simple('Ice Mass Lost (kg/day)', time[1:], loss_Tfines[1:], 'Fine Mass Lost', 'Evolution of Fine Ice Mass Loss',  'TotIceMassLoss_Fines.png', mainVfolder)

# #Plot Coarse Loss Mass
# plot_simple('Ice Mass Lost (kg)', time[1:], loss_Tcoarse[1:], 'Coarse Mass Lost', 'Evolution of Coarse Ice Mass Loss',  'TotIceMassLoss_Coarse.png', mainVfolder)

#Plot Coarse Loss Points
plot_simple('Ice Mass Lost (kg/day)', time[1:], loss_TcoarseA[1:], 'Coarse Mass Lost', 'Evolution of Coarse Ice Mass Loss Pts.',  'TotIceMassLoss_CoarseP.png', mainVfolder)

# #Plot Total Loss Mass
# plot_simple('Ice Mass Lost (kg)', time[1:], loss_Total[1:], 'Total Mass Lost', 'Evolution of Total Ice Mass Loss',  'TotIceMassLoss_Evo.png', mainVfolder)

#Plot Total Loss Points
plot_simple('Ice Mass Lost (kg/day)', time[1:], loss_TotalA[1:], 'Total Mass Lost', 'Evolution of Total Ice Mass Loss Pts.',  'TotIceMassLoss_EvoP.png', mainVfolder)

# #Compare Coarse Loss
# y_label = 'Ice Mass Lost (kg)'
# data_label = 'Coarse Mass Lost (_mass)'            
# data_label2 =  'Coarse Mass Lost (points)'                    
# plot_title = 'Evolution of Coarse Ice Mass Loss'
# File_Name = 'TotIceMassLoss_Coarse.png'
# dataV = loss_Tcoarse
# dataV2 = loss_TcoarseA
# fig = plt.figure(figsize=(10,10))
# plt.xlabel('Time (days)', fontsize = 20)
# plt.ylabel(y_label, fontsize = 20)
# plt.tick_params(axis='x', labelsize=15)
# plt.tick_params(axis='y', labelsize=15)
# plt.plot(time[1:], dataV[1:], 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
# plt.plot(time[1:], dataV2[1:], 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
# plt.legend(prop={'size': 20})
# plt.title(plot_title, fontsize = 25)
# plt.savefig(mainVfolder + File_Name, format='png')
# plt.close()  

# # #Compare Total Loss
# # plot_double('Ice Mass Lost (kg)', time, loss_Total, 'Total Mass Lost (_mass)', loss_TotalA, 'Total Mass Lost (points)', 'Evolution of Total Ice Mass Loss',  'TotIceMassLoss_Evo.png', mainVfolder)
# y_label = 'Ice Mass Lost (kg)'
# data_label = 'Total Mass Lost (_mass)'           
# data_label2 =  'Total Mass Lost (points)'                   
# plot_title = 'Evolution of Total Ice Mass Loss'
# File_Name = 'TotIceMassLoss_Evo.png'
# dataV = loss_Total
# dataV2 = loss_TotalA
# fig = plt.figure(figsize=(10,10))
# plt.xlabel('Time (days)', fontsize = 20)
# plt.ylabel(y_label, fontsize = 20)
# plt.tick_params(axis='x', labelsize=15)
# plt.tick_params(axis='y', labelsize=15)
# plt.plot(time[1:], dataV[1:], 'b', marker='s', fillstyle='none', label=data_label) #Same Size Data
# plt.plot(time[1:], dataV2[1:], 'r', marker='s', fillstyle='none', label=data_label2) #Same Size Data
# plt.legend(prop={'size': 20})
# plt.title(plot_title, fontsize = 25)
# plt.savefig(mainVfolder + File_Name, format='png')
# plt.close()  

# #Get coarse mass lost plus decayed coarse mass lost
# rho = 910e9
# decayMass = np.zeros((nSteps))
# for i in range(nSteps):
#     tsdFileLoc = mainOutputfolder + 'Vertical_Thickness/' + 'TSD_Fine_step_' + str(i) + '.dat'
#     tsdV = np.loadtxt(tsdFileLoc)
#     if len(tsdV) <= 100:
#         decayMass[i] = 0
#     else:
#         decaySum = tsdV[100:]
#         for j in range(len(tsdV) - 100):
#             decayMass[i] += decaySum[j][0] * decaySum[j][0] * 0.001 * rho
#     # if i > 0:
#     #     decayMass[i] = decayMass[i-1] + decayMass[i] #For accumulated mass loss

# lostDecayMass = np.zeros((nSteps))
# for i in range(nSteps):
#     if i > 0:
#         lostDecayMass[i] = decayMass[i] - decayMass[i-1]

# #Sum of massLoss_mcv and decayed (do not used mcl since it just goes to decayed)
# alter_mass_coarse_loss_cum =  lostDecayMass + loss_mcv

# #Plot Melt Alter
# plot_simple('Ice Mass Loss (kg)', time, alter_mass_coarse_loss_cum, 'Cumulative Coarse ONLY Vertical Melt + Decay', 'Cumulative Coarse ONLY Vertical Melt + Decay',  'Coarse_Vertical_Melt_and_Decay.png', mainVfolder)
