#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 18:06:31 2019

@author: macbookpro
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

fig=plt.figure(figsize=(20,20))
fig.subplots_adjust(top=0.92, bottom=0.01, hspace=0.4, wspace=0.25)

def PintaFluence(FeatureSize, Fluence, Sample, NumberPlot, labelsTag):
    
    ax=fig.add_subplot(3,3,NumberPlot) #number plot 1-2-3
    ax.set_title(Sample +'\n', fontsize=15)
    ax.set_ylabel(r'Fluence [ph/$\mu$m$^{2}$]', fontsize=12)  
    ax.set_xlabel(r'Feature size [nm]', fontsize=12) 
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim([10, 100])
    ax.set_ylim([10**11, 10**18])
    plt.plot(FeatureSize*10**(7), Fluence, linestyle ='-', label = labelsTag) #color = colors[i]
    plt.grid(True)
    plt.legend()
    
    return fig

def PintaDose(FeatureSize, Dose, Sample, NumberPlot, labelsTag):
    
    ax=fig.add_subplot(3,3,NumberPlot) #number plot 4-5-6
    ax.set_ylabel(r'Dose [Gy]', fontsize=12)  
    ax.set_xlabel(r'Feature size [nm]', fontsize=12) 
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim([10, 100])
    plt.plot(FeatureSize*10**(7), Dose, linestyle='-', label = labelsTag) #color = colors[i]
    plt.grid(True)
    plt.legend()
    
    return fig
  

Samples = ['Biomolecule', 'DNA', 'PMMA']  

labelsTag = ['Ideal det. (5um cell)', 'Optimistic det. (5um cell)', 'Conservative det. (5um cell)']
labelsTag_complete_bkgs = ['Ideal det. (5um cell + medium)', 'Optimistic det. (5um cell + medium)', 'Conservative (5um cell + medium)']

######################### With water equivalent for background ###########################################

file_path_name_biomolecule = ['./OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Ideal.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Optimistic.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Realistic.txt']

file_path_name_DNA = ['./OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Ideal.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Optimistic.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Realistic.txt']

file_path_name_PMMA = ['./OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Ideal.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Optimistic.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Realistic.txt']

######################### With water equivalent and air for background --> Complete situation ###############

file_path_name_biomolecule_complete_bkgs = ['./OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Ideal_bkg_complete.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Optimistic_bkg_complete.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_Realistic_bkg_complete.txt']

file_path_name_DNA_complete_bkgs = ['./OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Ideal_bkg_complete.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Optimistic_bkg_complete.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_Realistic_bkg_complete.txt']

file_path_name_PMMA_complete_bkgs = ['./OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Ideal_bkg_complete.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Optimistic_bkg_complete.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_Realistic_bkg_complete.txt']

'''
########################## BIOMOLECULE Air bkg ####################################

array_txt = np.loadtxt(file_path_name_biomolecule[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[0])

array_txt = np.loadtxt(file_path_name_biomolecule[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[1])

array_txt = np.loadtxt(file_path_name_biomolecule[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[2])


############################ DNA Air bkg #######################################

array_txt = np.loadtxt(file_path_name_DNA[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[0])

array_txt = np.loadtxt(file_path_name_DNA[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[1])

array_txt = np.loadtxt(file_path_name_DNA[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[2])


########################### PMMA Air bkg ######################################

array_txt = np.loadtxt(file_path_name_PMMA[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[0])

array_txt = np.loadtxt(file_path_name_PMMA[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[1])

array_txt = np.loadtxt(file_path_name_PMMA[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[2])
'''
##############################################################################################
##############################################################################################

########################## BIOMOLECULE Complete situation ####################################
'''
array_txt = np.loadtxt(file_path_name_biomolecule_complete_bkgs[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[0])

array_txt = np.loadtxt(file_path_name_biomolecule_complete_bkgs[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[1])

array_txt = np.loadtxt(file_path_name_biomolecule_complete_bkgs[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[0], 1, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 4, labelsTag[2])


############################ DNA Complete situation #######################################

array_txt = np.loadtxt(file_path_name_DNA_complete_bkgs[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[0])

array_txt = np.loadtxt(file_path_name_DNA_complete_bkgs[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[1])

array_txt = np.loadtxt(file_path_name_DNA_complete_bkgs[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[1], 2, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 5, labelsTag[2])


########################### PMMA Complete situation ######################################

array_txt = np.loadtxt(file_path_name_PMMA_complete_bkgs[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[0])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[0])

array_txt = np.loadtxt(file_path_name_PMMA_complete_bkgs[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[1])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[1])

array_txt = np.loadtxt(file_path_name_PMMA_complete_bkgs[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Fluence = array_txt[:,1] #ph/mum2
Dose = array_txt[:,2] #Gy
figura = PintaFluence(FeatureSize, Fluence, Samples[2], 3, labelsTag[2])
figura = PintaDose(FeatureSize, Dose, Samples[0], 6, labelsTag[2])
'''

######################## Result independent plot with Villanueva data and ZOOM #################################

array_txt = np.loadtxt('./Data/Villanueva_values.txt', skiprows=1) #Tendria que volver a escanear los datos de Pablo para corregir el inicio
FeatureSize = array_txt[:,0] #Size (cm)
Dose = array_txt[:,1] #Gy

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

A,B = curve_fit(f, FeatureSize, Dose)[0] # your data x, y to fit

Array_x = np.linspace(1, 1000, 1000, endpoint = True)
Array_f = f(Array_x, A, B)

fig, ax1 = plt.subplots()
ax1.plot(Array_x, Array_f, 'k', label = 'Maximun tolerable dose (Villanueva et al.)')

####

array_txt_1 = np.loadtxt(file_path_name_biomolecule[0], skiprows=1)
FeatureSize_1 = array_txt_1[:,0] #Size (cm)
Dose_1 = array_txt_1[:,2] #Gy
ax1.plot(FeatureSize_1*10**(7), Dose_1, linestyle='-', label = labelsTag[0])

array_txt_2 = np.loadtxt(file_path_name_biomolecule[1], skiprows=1)
FeatureSize_2 = array_txt_2[:,0] #Size (cm)
Dose_2 = array_txt_2[:,2] #Gy
ax1.plot(FeatureSize_2*10**(7), Dose_2, linestyle='-', label = labelsTag[1])

array_txt_3 = np.loadtxt(file_path_name_biomolecule[2], skiprows=1)
FeatureSize_3 = array_txt_3[:,0] #Size (cm)
Dose_3 = array_txt_3[:,2] #Gy
ax1.plot(FeatureSize_3*10**(7), Dose_3, linestyle='-', label = labelsTag[2])

####

array_txt_1_bkg = np.loadtxt(file_path_name_biomolecule_complete_bkgs[0], skiprows=1)
FeatureSize_1_bkg = array_txt_1_bkg[:,0] #Size (cm)
Dose_1_bkg = array_txt_1_bkg[:,2] #Gy
ax1.plot(FeatureSize_1_bkg*10**(7), Dose_1_bkg, linestyle='--', label = labelsTag_complete_bkgs[0])

array_txt_2_bkg = np.loadtxt(file_path_name_biomolecule_complete_bkgs[1], skiprows=1)
FeatureSize_2_bkg = array_txt_2_bkg[:,0] #Size (cm)
Dose_2_bkg = array_txt_2_bkg[:,2] #Gy
ax1.plot(FeatureSize_2_bkg*10**(7), Dose_2_bkg, linestyle='--', label = labelsTag_complete_bkgs[1])

array_txt_3_bkg = np.loadtxt(file_path_name_biomolecule_complete_bkgs[2], skiprows=1)
FeatureSize_3_bkg = array_txt_3_bkg[:,0] #Size (cm)
Dose_3_bkg = array_txt_3_bkg[:,2] #Gy
ax1.plot(FeatureSize_3_bkg*10**(7), Dose_3_bkg, linestyle='--', label = labelsTag_complete_bkgs[2])

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel(r'Sample: biomolecule - Dose [Gy]', fontsize=12)  
ax1.set_xlabel(r'Feature size [nm]', fontsize=12) 
ax1.grid(True)
plt.gca().set_xticks(np.arange(20.0, 80.0, 10.0), minor=True)
ax1.xaxis.grid(True, which='minor') 
ax1.legend(loc=2, fontsize = 12)
ax1.set_xlim([10.0, 100])
ax1.set_ylim([10**8, 10**13])

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.55,0.62,0.5,0.5])
ax2.set_axes_locator(ip)
# Mark the region corresponding to the inset axes on ax1 and draw lines
# in grey linking the two axes.
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

# The data into the ax2
ax2.plot(Array_x, Array_f, 'k', label = 'Maximun tolerable dose (Villanueva et al.)')
ax2.plot(FeatureSize_1*10**(7), Dose_1, linestyle='-', label = labelsTag[0])
ax2.plot(FeatureSize_2*10**(7), Dose_2, linestyle='-', label = labelsTag[1])
ax2.plot(FeatureSize_3*10**(7), Dose_3, linestyle='-', label = labelsTag[2])
ax2.plot(FeatureSize_1_bkg*10**(7), Dose_1_bkg, marker='o', label = labelsTag[0])
ax2.plot(FeatureSize_2_bkg*10**(7), Dose_2_bkg, marker='o', label = labelsTag[1])
ax2.plot(FeatureSize_3_bkg*10**(7), Dose_3_bkg, marker='o', label = labelsTag[2])

#ax2.set_yscale('log')
#ax2.set_xscale('log')
ax2.set_xlim([30.0, 50.0])
ax2.set_ylim([2*10**9, 6*10**9])
ax2.set_xticks([30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0])
ax2.grid(True, which='both')

plt.show()
