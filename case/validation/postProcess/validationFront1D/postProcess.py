# %%
################################################################################
'''
Purpose: Data processing & Plotting 1D line figures
Requirement: Anaconda
Author: Fengchang Yang
Date: 09/04/2018
'''
################################################################################
# Import useful modules
from __future__ import division
import numpy as np
import numba as nb
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from scipy.integrate import simps
from scipy.integrate import trapz

################################################################################
# %%
# import data
legendList = ['simulation', 'analytical']
fileList = ['validationFront2D']
solidMass = []

for idx1, x in enumerate(fileList):
    df = pd.read_csv('/home/fy0/OpenFOAM/fy0-v1912/run/'+x+'/outputResults.csv', skiprows=0)
    solidMass.append(np.transpose(np.array(df.loc[:,['Time[s]',' solidMass[kg]']])))
solidMass = solidMass[0]

# %%
# check energy
tAna = np.linspace(0,10000,101)
D = 1.0
k = 0.001
L0 = 1.0
S0 = 1.0
C0 = 16
Ceq = 2
rhom = 16
solidMassAna = (D/k+L0)-np.sqrt((D/k+L0-S0)**2-2*D*(C0-Ceq)*tAna/rhom)
solidMass[1] = solidMass[1]/rhom/4

df = pd.DataFrame({'Time[s]':tAna, 'solidMass[kg]':solidMassAna})
df.to_csv('solidMassValidationAnalytical.csv')

# %%
# Custom greek font in label
# plt.rc('font', family='serif')
font = {'family': 'Times New Roman',
        'weight': 'normal'}
plt.rc('font', **font)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'

# RGB color schemes from matlab
colors = [(0, 0.447, 0.7410), (0.8500, 0.3250, 0.0980),
          (0.9290, 0.6940, 0.1250)]

# Figure parameters
fsize = 20
hsize = 16
fwidth = 6.5
fwidthlong = 6.5 * 2
fheight = 5
fwidthlong = 6.5 * 2
markStep = [200, 20]
markerList = ['o', '^', 's', 'v', '<', '>', '*', 'p']
axisList = [700, 600, 200, 250]

# %%
# plot figures
# figure 1
fig = plt.figure(figsize=(fwidth, fheight))
plt.plot(tAna, solidMassAna, 'k-', fillstyle='none', label=legendList[1])
plt.plot(solidMass[0][::2], solidMass[1][::2], markerList[0], fillstyle='none', label=legendList[0])
plt.xlabel('Time (s)', fontname='Times New Roman', fontsize=fsize)
plt.ylabel(r'Front position, S',
            fontname='Times New Roman', fontsize=fsize)
H = plt.legend(loc='lower right', prop={'size': 14}, numpoints=1,
               frameon=False)
plt.axis([0, 10000, 1, 10])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
# save figure
plt.tight_layout()
plt.savefig('mpFoamValidation1_2D' + '.png', dpi=300)
plt.show()
