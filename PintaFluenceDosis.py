#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 18:06:31 2019

@author: macbookpro
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

'''
fig=plt.figure(figsize=(20,20))
fig.subplots_adjust(top=0.92, bottom=0.01, hspace=0.4, wspace=0.25)

def PintaFluence(FeatureSize, Fluence, Sample, NumberPlot, labelsTag):
    
    ax=fig.add_subplot(3,3,NumberPlot) #number plot 1-2-3
    ax.set_title(Sample +'\n', fontsize=15)
    ax.set_ylabel(r'Fluence [ph/$\mu$m$^{2}$]', fontsize=12)  
    ax.set_xlabel(r'Feature size [nm]', fontsize=12) 
    ax.set_yscale('log')
    ax.set_xlim([10, 100])
    #ax.set_xscale('log')
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
'''    

Samples = ['Biomolecule', 'DNA', 'PMMA']  
  
file_path_name_biomolecule = ['./OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_IdealDetector.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_100-45-5.txt', \
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_50-25-5.txt'\
                              './OutputFiles/Biomolecule_Fluence_Dose_vs_FeatureSize_fake.txt']


file_path_name_DNA = ['./OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_IdealDetector.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_100-45-5.txt', \
                              './OutputFiles/DNA_Fluence_Dose_vs_FeatureSize_50-25-5.txt']

file_path_name_PMMA = ['./OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_IdealDetector.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_100-45-5.txt', \
                              './OutputFiles/PMMA_Fluence_Dose_vs_FeatureSize_50-25-5.txt']

labelsTag = ['Ideal detector', 'Realistic detector', 'Conservative detector']

'''
########################## BIOMOLECULE ####################################

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


############################ DNA #######################################

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


########################### PMMA ######################################

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
######################## Result independent plot with Villanueva data #################################

array_txt = np.loadtxt('Villanueva_values.txt', skiprows=1) #Tendria que volver a escanear los datos de Pablo para corregir el inicio
FeatureSize = array_txt[:,0] #Size (cm)
Dose = array_txt[:,1] #Gy

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

A,B = curve_fit(f, FeatureSize, Dose)[0] # your data x, y to fit

Array_x = np.linspace(1, 1000, 1000, endpoint = True)
Array_f = f(Array_x, A, B)


figureFinal=plt.figure(figsize=(10,5))
figureFinal.subplots_adjust(top=0.90, bottom=0.1, hspace=0.4, wspace=0.25)
plt.title('Biomolecule', fontsize = 15)
axF=figureFinal.add_subplot(111)
axF.set_ylabel(r'Dose [Gy]', fontsize=12)  
axF.set_xlabel(r'Feature size [nm]', fontsize=12) 
plt.plot(Array_x, Array_f, 'k', label = 'Maximun tolerable dose (Villanueva et al.)')
axF.set_yscale('log')
axF.set_xscale('log')
axF.set_xlim([10.0, 100])
axF.set_ylim([10**8, 10**13])

array_txt = np.loadtxt(file_path_name_biomolecule[0], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Dose = array_txt[:,2] #Gy
plt.plot(FeatureSize*10**(7), Dose, linestyle='-', label = labelsTag[0])

array_txt = np.loadtxt(file_path_name_biomolecule[1], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Dose = array_txt[:,2] #Gy
plt.plot(FeatureSize*10**(7), Dose, linestyle='-', label = labelsTag[1])

array_txt = np.loadtxt(file_path_name_biomolecule[2], skiprows=1)
FeatureSize = array_txt[:,0] #Size (cm)
Dose = array_txt[:,2] #Gy
plt.plot(FeatureSize*10**(7), Dose, linestyle='-', label = labelsTag[2])


plt.grid(True, which='both')
plt.legend()



 