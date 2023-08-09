import numpy as np
import h5py

import time

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import LineCollection
from sys import argv

from scipy.optimize import curve_fit
from scipy import fftpack
import os.path
from alive_progress import alive_bar

# objective function for linear fitting
def objective(x, a, b):
    return a * x + b

#Interpret all different seed mass results
def process_mass(massV):

    #Prepare data
    case = np.zeros(len(massV))
    mass_loss = np.zeros(len(massV))
    for i in range(len(massV)):
        case[i] = massV[i][0]
        mass_loss[i] = massV[i][1]

    #Define if jammed or not jammed

    jam_state = 5 #Declare jamming if last 5 results are jammed
    tol = 0.0000000000001 #Declare jammed if flow is below this value

    gradV = np.gradient(mass_loss) #Derivative of mass flow
    tail = gradV[-jam_state:] #Last values to evaluate
    jammed = 1

    #Check tail derivative and if not all are jammed declare unjammed
    for i in range(len(tail)):
        if tail[i]>tol:
            jammed = 0
            break

    #Simply average all derivatives to get an average mass loss rate        
    meanD = np.mean(gradV)       

    #Get derivative using linear fit of line as mx + b = y
    poptMean, _ = curve_fit(objective, case, mass_loss)
    aMean = poptMean[0] #Slope
    bMean = poptMean[1]
    
    meanF = aMean

    return jammed, meanF, meanD
    
#Interpret all different seed KE results
def process_KE(keV):

    #Prepare data
    case = np.zeros(len(keV))
    KE = np.zeros(len(keV))
    for i in range(len(keV)):
        case[i] = keV[i][0]
        KE[i] = keV[i][1]

    #Define if jammed or not jammed

    jam_state = 20 #Declare jamming if last N results are jammed
    tol = 5 #Declare jammed if KE is below this value

    gradV = np.gradient(KE) #Derivative of mass flow
    tail = gradV[-jam_state:] #Last values to evaluate
    jammed = 1

    #Check tail derivative and if not all are jammed declare unjammed
    for i in range(len(tail)):
        if tail[i]>tol:
            jammed = 0
            break

    #Simply average all derivatives to get an average KE 
    meanD = np.mean(gradV)       

    #Get derivative using linear fit of line as mx + b = y
    #poptMean, _ = curve_fit(objective, case, mass_loss)
    #aMean = poptMean[0] #Slope
    #bMean = poptMean[1]
    
    #meanF = aMean

    return jammed, meanD



#Main folder of results
main_folder = str(argv[1])+'/'        #e.g. #mass_loss_Hratio_""2.5""
sub_loc = 'mass_loss_'+str(argv[2])   #e.g. 'mass_loss_r-6e-3_s-'
output_folder = str(argv[3])+'/'      #e.g. mass_loss_KG_0.5
H_ratio = str(argv[4])
angle = str(argv[5])

#Result list
seed_list = np.linspace(1, 100, 100) #Remember to conver to integer

#Output Results for Average and Plotting of the N seeds for a specific case
cases = np.zeros(len(seed_list))
jam = np.zeros(len(seed_list))
jamKE = np.zeros(len(seed_list))
ave_rateF = np.zeros(len(seed_list))
ave_rateD = np.zeros(len(seed_list))
ave_rateDKE = np.zeros(len(seed_list))

#Process each seed output mass file
miss_countA = 0
missVA = []
miss_countB = 0
missVB = []
with alive_bar(len(seed_list)) as bar: ##PROGRESS BAR
    for i in range(len(seed_list)):
        cases[i] = round(int(seed_list[i]),0)
        
        #file_loc = main_folder+sub_loc+str(int(cases[i]))+'.dat'  #15
        #file_loc = main_folder+sub_loc+str(int(cases[i]))+'_a.dat'   #30
        file_loc = main_folder+sub_loc+str(int(cases[i]))+'_a-'+angle+'.dat'   #45
        if os.path.exists(file_loc) == False or len(np.loadtxt(file_loc, delimiter = " ")) < 1 :
            miss_countA += 1
            missVA.append(cases[i])
            file_loc = main_folder+sub_loc+str(1)+'_a-'+angle+'.dat'
            
        file_loc2 = 'examples_output/hopper/nofrac__r-6e-3_h-'+H_ratio+'_s-'+str(int(cases[i]))+'_a-'+angle+'/ke.dat'   #45
        if os.path.exists(file_loc2) == False or  len(np.loadtxt(file_loc2, delimiter = " ")) < 1 :
            miss_countB += 1
            missVB.append(cases[i])
            file_loc2 = 'examples_output/hopper/nofrac__r-6e-3_h-'+H_ratio+'_s-'+str(1)+'_a-'+angle+'/ke.dat'   #45
    
        massV = np.loadtxt(file_loc, delimiter = " ")
        keV = np.loadtxt(file_loc2, delimiter = " ")
        [jam[i], ave_rateF[i], ave_rateD[i]] = process_mass(massV)
        [jamKE[i], ave_rateDKE[i]] = process_KE(keV)
        bar() ##PROGRESS BAR

#For bar plots
#Average all seeds for output
mJam = np.mean(jam)
mJamKE = np.mean(jamKE)
mRateF = np.mean(ave_rateF)
mRateD = np.mean(ave_rateD)

#Standard deviation
sJam = np.std(jam)
sRateF = np.std(ave_rateF)
sRateD = np.std(ave_rateD)

#Min Value
minRateF = np.min(ave_rateF)
minRateD = np.min(ave_rateD)

#Max Value
maxRateF = np.max(ave_rateF)
maxRateD = np.max(ave_rateD)

#Plot all cases for comparison

#Jam or Not
out_png = output_folder+'Seed_Jam_Hratio_'+str(argv[4])+'.png'
plt.plot(cases, jam, 'b*')
plt.xlabel('Case')
plt.ylabel('Jammed')
plt.title('Jam for KG: ' + str(argv[4]))
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

#Average Loss Rate (Linear Fit)
out_png = output_folder+'Seed_AveRF_Hratio_'+str(argv[4])+'.png'
plt.plot(cases, ave_rateF, 'b*')
plt.xlabel('Case')
plt.ylabel('Linear Fit Av. Rate')
plt.title('Average Rate (Linear Fit) for Hratio: ' + str(argv[4]))
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

#Average Loss Rate (Derivative Average)
out_png = output_folder+'Seed_AveRD_Hratio_'+str(argv[4])+'.png'
plt.plot(cases, ave_rateD, 'b*')
plt.xlabel('Case')
plt.ylabel('Derivative Av. Rate')
plt.title('Average Rate (Derivative) for Hratio: ' + str(argv[4]))
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

#Output results for additional processing
dataV = []
#for ii in range(len(tt)):
dataV.append([float(argv[4]), mJam, mRateF, mRateD, sJam, sRateF, sRateD, minRateF, minRateD, maxRateF, maxRateD, mJamKE])
np.savetxt(output_folder+'Results_Hratio_'+str(argv[4])+'.dat', dataV, fmt='%8.8f')

#KE
#Jam or Not using KE
out_png = output_folder+'Seed_KE_Jam_Hratio_'+str(argv[4])+'.png'
plt.plot(cases, jamKE, 'b*')
plt.xlabel('Case')
plt.ylabel('Jammed')
plt.title('Jam using KE for Hratio: ' + str(argv[4]))
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

#Average Loss Rate (Derivative Average)
out_png = output_folder+'Seed_KE_AveRD_Hratio_'+str(argv[4])+'.png'
plt.plot(cases, ave_rateDKE, 'b*')
plt.xlabel('Case')
plt.ylabel('Derivative Av. Rate')
plt.title('Average Rate (Derivative) KE for Hratio: ' + str(argv[4]))
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

if miss_countA > 0:
    print('mass loss WARNING MISSING DATA:', missVA)
    print('# MISSING DATA:', miss_countA)
if miss_countB > 0:
    print('KE WARNING MISSING DATA:', missVB)
    print('# MISSING DATA:', miss_countB)
