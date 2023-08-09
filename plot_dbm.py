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
import scipy
from scipy.optimize import curve_fit
import sys
import shutil

import matplotlib.ticker as mticker
import statistics

from PIL import Image

#import shapely.geometry

import math
import datetime
import csv
import cv2
import numpy as np
import glob

dataV = np.loadtxt("dbm_data.dat")

lenB = 7
lenq = 6

Bvec = np.zeros(lenB)
qvec = np.zeros(lenq)
dbmvec = np.zeros(len(dataV))
dbmvec2 = np.zeros(len(dataV))

for i in range(lenB):
    Bvec[i] = dataV[i][1]
for i in range(lenq):
    qvec[i] = dataV[i*(lenB)][2]
ct = 0
for i in range(lenB):
    for j in range(lenq):
        dbmvec[ct] = dataV[ct][3]
        dbmvec2[ct] = dataV[ct][4]
        ct += 1

bf_ind = 17
hb_ind = 18
lb_ind = 21

#DBM Plot LOG
plt.figure(figsize=(10,10)) #18,18 ok
plt.xlabel(r'Breakage rate ($d^{-1}$)', size=20)
plt.ylabel(r'$d_{bm}$ (km)', size=20)
plt.title(r'Transition Diameter vs. B and $q_v$', size=20)
#plt.plot(gsdD2r, NLMeltr, 'rs-', label='Number Loss Melt')
#plt.plot(gsdD2r, NLBreakr, 'bs-', label='Number Loss Break')

for i in range(len(qvec)-1,-1,-1):
    if qvec[i]>1:
        plt.plot(Bvec, dbmvec[i*lenB:(i+1)*lenB], 's-', label=r'$q_{v} =$ ' + str(qvec[i]) )
    else:
        plt.plot(Bvec, dbmvec[i*lenB:(i+1)*lenB], 'k--', label=r'$q_{v} =$ ' + str(qvec[i]) + " (No Solar Term)" )

plt.plot(Bvec[3], dbmvec[bf_ind], 'k*', label='Best Fit', markersize=20)
plt.plot(Bvec[4], dbmvec[hb_ind], 'b*', label='High Break', markersize=20)
plt.plot(Bvec[0], dbmvec[lb_ind], 'r*', label='Low Break', markersize=20)
#plt.axvline(x=dbm_1, color='k', label=r'$d_{bm}$ = ' + str(round(dbm_1,1)) + " km", ls='--')
#plt.xlim( math.log(1), math.log(30) )
#plt.xlim(0,30)
#plt.ylim(-200,200)
#plt.ylim( math.log(1), 200 )
# plt.ylim(tolV,200)
# # plt.ylim(-200,200)
#plt.ylim(-2,2)
#plt.ylim(0.00007,100)
plt.yscale('log')
# plt.xscale('log')
#plt.ylim(0,math.log(200))
#plt.legend(fontsize = 18, framealpha=1, loc='upper right')
plt.legend(fontsize = 18, loc='upper right')
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax=plt.gca()
#ax.set_xticks([5,10,20,30,40,50])
ax.set_yticks([5,10,20,30,40,50])
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
# ax.xaxis.get_major_formatter().set_scientific(False)
# ax.xaxis.get_major_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
# ax.yaxis.get_major_formatter().set_scientific(False)
# ax.yaxis.get_major_formatter().set_useOffset(False)
#ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

plt.savefig("dbm_B_qvLOG.png", format='png')
plt.close()

#DBM Plot
plt.figure(figsize=(10,10)) #18,18 ok
plt.xlabel(r'Breakage rate ($d^{-1}$)', size=20)
plt.ylabel(r'$d_{bm}$ (km)', size=20)
plt.title(r'Transition Diameter vs. B and $q_v$', size=20)
#plt.plot(gsdD2r, NLMeltr, 'rs-', label='Number Loss Melt')
#plt.plot(gsdD2r, NLBreakr, 'bs-', label='Number Loss Break')

for i in range(len(qvec)):
    if qvec[i]>1:
        plt.plot(Bvec, dbmvec[i*lenB:(i+1)*lenB], '*-', label=r'$q_{v} =$ ' + str(qvec[i]) )
    else:
        plt.plot(Bvec, dbmvec[i*lenB:(i+1)*lenB], '*-', label=r'$q_{v} =$ ' + str(qvec[i]) + " (No Solar Term)" )

plt.plot(Bvec[3], dbmvec[bf_ind], 'k*', label='Best Fit', markersize=20)
plt.plot(Bvec[4], dbmvec[hb_ind], 'b*', label='High Break', markersize=20)
plt.plot(Bvec[0], dbmvec[lb_ind], 'r*', label='Low Break', markersize=20)
#plt.axvline(x=dbm_1, color='k', label=r'$d_{bm}$ = ' + str(round(dbm_1,1)) + " km", ls='--')
#plt.xlim( math.log(1), math.log(30) )
#plt.xlim(0,30)
#plt.ylim(-200,200)
#plt.ylim( math.log(1), 200 )
# plt.ylim(tolV,200)
# # plt.ylim(-200,200)
#plt.ylim(-2,2)
#plt.ylim(0.00007,100)
#plt.yscale('log')
# plt.xscale('log')
#plt.ylim(0,math.log(200))
#plt.legend(fontsize = 18, framealpha=1, loc='upper right')
plt.legend(fontsize = 18, loc='upper right')
plt.tick_params(labelsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax=plt.gca()
#ax.set_xticks([5,10,20,30,40,50])
ax.set_yticks([5,10,20,30,40,50])
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
# ax.xaxis.get_major_formatter().set_scientific(False)
# ax.xaxis.get_major_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
# ax.yaxis.get_major_formatter().set_scientific(False)
# ax.yaxis.get_major_formatter().set_useOffset(False)
#ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

plt.savefig("dbm_B_qv.png", format='png')
plt.close()


