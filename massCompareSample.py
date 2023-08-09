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



nBreakRange = [1/200, 1/250, 1/400, 1/1000, 1/1250, 1/1500, 1/1700, 1/2000]
qvertRange = [5, 10, 15, 20, 25, 30, 35, 40, 45]

#Example for S = 310 for 2018
inputFile = "./V_P_2018/MassLossAnalysis_mass_coarse_loss/results_310.dat"
inputFileV = np.loadtxt(inputFile)

#outputFolder2 = "./Output/Analysis/"
outputFolder2 = "./Output/Analysis_2018_mass_loss/"
#outputFolder2 = "./Output/Analysis_2018Prob/" #PROB CHANGE
#outputFolder2 = "./Output/Analysis_2018h/"
if os.path.exists(outputFolder2) == False:
	os.makedirs(outputFolder2)

#For MultiD
Vec_MultiD = []

ct = 0
for k in range(len(nBreakRange)):
    for j in range(len(qvertRange)) :
        multiV = [inputFileV[ct][2], 1/inputFileV[ct][1], inputFileV[ct][4]/1e12, inputFileV[ct][5]/1e12, inputFileV[ct][6], inputFileV[ct][0]];  #qv, nb, break_loss, melt_loss, ratio_bk_melt
        Vec_MultiD.append(multiV)
        ct = ct + 1

#MULTID
print(Vec_MultiD) 
#exit(1)

#TODO
#Contour MultiD Plot for Vec_MultiD

x = np.zeros(( len(qvertRange) * len(nBreakRange) ))
x = np.zeros(( len(Vec_MultiD) ))
y = np.zeros((len(x)))
bkloss = np.zeros((len(x)))
meltloss = np.zeros((len(x)))
ratio = np.zeros((len(x)))
case_grid  = np.zeros((len(x)))

print(Vec_MultiD[0][0]) #CHECK
print("Len A: ", len(qvertRange) * len(nBreakRange)) #CHECK
print("Len B: ", len(Vec_MultiD)) #CHECK

for i in range(len(x)): #CHECK
    x[i] = (Vec_MultiD[i][0])
    y[i] = (Vec_MultiD[i][1])
    bkloss[i] = (Vec_MultiD[i][2])
    meltloss[i] = (Vec_MultiD[i][3])
    ratio[i] = (Vec_MultiD[i][4])
    case_grid[i] = (Vec_MultiD[i][5])


#Dimension of grid for plotting
xgrid = len(qvertRange)
ygrid = len(nBreakRange)

#Multiple
#Contour Plots
ct = 0
xArr = np.zeros((xgrid,ygrid))
yArr = np.zeros((xgrid,ygrid))
bklossArr = np.zeros((xgrid,ygrid))
meltlossArr = np.zeros((xgrid,ygrid))
ratioArr = np.zeros((xgrid,ygrid))
caseArr = np.zeros((xgrid,ygrid))

for i in range(xgrid):
    for j in range(ygrid):

        #Select from MultiV correctly
        xArr[i][j] = x[ct]
        yArr[i][j] = y[ct]
        bklossArr[i][j] = bkloss[ct]
        meltlossArr[i][j] = meltloss[ct]
        ratioArr[i][j] = ratio[ct]
        caseArr[i][j] = case_grid[ct]
        ct += 1

print(xArr)
print(yArr)
print(ratioArr)

fig = plt.figure(figsize=(10,10))
plt.title('Coarse Bkg/Melt Ratio', size=20)
bx = [25]
by = [1/1000]
plt.plot(bx,by,'k*', markersize=20) #Optimal case for field
contours = plt.contourf(xArr, yArr, ratioArr, 15, cmap='coolwarm')
cbar = plt.colorbar();
#cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Bkg/Melt Mass Loss Ratio', rotation=270, labelpad=15)
CS = plt.contour(xArr, yArr, ratioArr)
plt.clabel(CS, fontsize=10)
CS2 = plt.contour(xArr, yArr, ratioArr, levels = [1.00000001, 1.00001], colors=('k',),linestyles=('-',),linewidths=(2,))
plt.clabel(CS2, fmt = '%2.1d', colors = 'k', fontsize=14)
#plt.imshow(zSlopeCArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
plt.savefig(outputFolder2 + "Results_massloss_bkgmeltratio_Contour.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.title('Breakage Loss', size=20)
contours = plt.contourf(xArr, yArr, bklossArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Breakage Loss (1e12 kg)', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
plt.savefig(outputFolder2 + "Results_bkgmassloss_Contour.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.title('Vertical melt Loss', size=20)
contours = plt.contourf(xArr, yArr, meltlossArr, 15, cmap='coolwarm')
#plt.imshow(zSlopeFArr, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdGy', alpha=0.5)
#plt.clabel(contours, inline=True, fontsize=8, colors='k')
cbar = plt.colorbar();
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Melt Loss (1e12 kg)', rotation=270, labelpad=15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.gca().set_aspect('equal', adjustable='box')
axes = plt.gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
plt.savefig(outputFolder2 + "Results_meltmassloss_Contour.png")
plt.close()


#Point Plots

fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = ratio, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Coarse Bkg/Melt Ratio', size=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(ratio):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_massloss_bkgmeltratio_Points.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = bkloss, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Breakage Loss', size=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(bkloss):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_bkgmassloss_Points.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = meltloss, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Vertical Melt Loss', size=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(meltloss):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_meltmassloss_Points.png")
plt.close()

fig = plt.figure(figsize=(10,10))
plt.scatter(x,y, c = case_grid, s = 50)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.title('Vertical Melt Loss', size=15)
plt.xlabel('qvert (W/m2/Celsius)')
plt.ylabel('Break Frequency (1/s)')
axes = plt. gca()
axes.xaxis.label.set_size(20)
axes.yaxis.label.set_size(20)
axes.tick_params(labelsize=15)
for i, txt in enumerate(case_grid):
    axes.annotate(round(txt,3), (x[i], y[i]))
plt.savefig(outputFolder2 + "Results_meltmassloss_Points_casev.png")
plt.close()
