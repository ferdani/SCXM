#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:44:36 2018

@author: macbookpro

GENERATOR OF RANDOM NUMBERS ACORDING TO DISTRIBUTION FUNCTION
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rd
import sys
import Functions as Fun

'''    
#---------------------------- Generating random numbers----------------------------------------------------------------
'''
def RandomNumberDistribution(C, n, InitEnergy, thetaRad, Z, SinterpolantSmall, SinterpolantBig, Finterpolant):
    #define the interval (theta_min < theta_max)
    theta_min = np.min(thetaRad)
    theta_max = np.max(thetaRad)

    #find the maximun and the minimun of the distribution function
    ymin = 0.0
    ymax = 0.0

    numSteps = 100000
    XR=np.array([])
    YR=np.array([])
    ymin = 0.0 #start value
    ymax = 0.0 #start value

    for i in range(numSteps):
        x = theta_min + (theta_max - theta_min) * float(i) /numSteps
        try:
            if C==1: y = Fun.KleinNishinaCrossSection(InitEnergy, x)
            if C==2: y = Fun.IncoherentCrossSection(InitEnergy, x, SinterpolantSmall, SinterpolantBig)
            if C==3: y = Fun.ThompsonCrossSection(x)
            if C==4: y = Fun.CoherentCrossSection(InitEnergy, x, Z, Finterpolant)
            if C==5: y = Fun.KNpolarization(InitEnergy, x, 0.0) #plane phi = 0.0 degrees
            if C==6: y = Fun.THpolarization(x, 0.0) #plane phi = 0.0 degrees
        except:
            continue  
        if y < ymin: ymin = y
        if y > ymax: ymax = y 

    ymax = ymax
    ymin = ymin 

    #using an uniform random numbers distribution
    XR=np.array([])
    YR=np.array([])

    print('\n') 
    print('Constructing the random photon events...')

    for i in range(n):
        while True:
            xr = rd.uniform(0,1)
            yr = rd.uniform(0,1) 
            x = theta_min + (theta_max - theta_min) * xr     
            y = ymin + (ymax - ymin) * yr 
            try:
                if C==1: Fun.KleinNishinaCrossSection(InitEnergy, x)
                if C==2: Fun.IncoherentCrossSection(InitEnergy, x, SinterpolantSmall, SinterpolantBig)
                if C==3: Fun.ThompsonCrossSection(x)
                if C==4: Fun.CoherentCrossSection(InitEnergy, x, Z, Finterpolant)
                if C==5: Fun.KNpolarization(InitEnergy, x, 0.0) #plane phi = 0.0 degrees
                if C==6: Fun.THpolarization(x, 0.0) #plane phi = 0.0 degrees

            except:
                continue
        
            if C==1: 
                if y <= Fun.KleinNishinaCrossSection(InitEnergy, x):
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break 
            if C==2:        
                if y <= Fun.IncoherentCrossSection(InitEnergy, x, SinterpolantSmall, SinterpolantBig):
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break 
            if C==3:        
                if y <= Fun.ThompsonCrossSection(x):
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break 
            if C==4:        
                if y <= Fun.CoherentCrossSection(InitEnergy, x, Z, Finterpolant):
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break
            if C==5:
                if y <= Fun.KNpolarization(InitEnergy, x, 0.0): #plane phi = 0.0 degrees
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break       
            if C==6:
                if y <= Fun.THpolarization(x, 0.0): #plane phi = 0.0 degrees
                    if y >= 0:
                        XR=np.append(XR,x)
                        YR=np.append(YR,y*10**24) #units in barns
                        break        
                
            sys.stdout.write("\r{0}%".format(round(((float(i)/n)*100),1)))
            sys.stdout.flush()

    return [XR, YR]

'''
---------------------------------- Build the real distribution functions to plot then -------------------------------
'''
def RealDistributionFunctionNumbers(C, InitEnergy, thetaRad, Z, SinterpolantSmall, SinterpolantBig, Finterpolant):
    #Build data with the real functions to plot it
    SigmaArray = np.array([])

    for i in range(0, len(thetaRad)):
        if C == 1:
            value = Fun.KleinNishinaCrossSection(InitEnergy, thetaRad[i])
            SigmaArray = np.append(SigmaArray, value*10**24)
        if C == 2:
            value = Fun.IncoherentCrossSection(InitEnergy, thetaRad[i], SinterpolantSmall, SinterpolantBig)
            SigmaArray = np.append(SigmaArray, value*10**24)
        if C == 3:
            value = Fun.ThompsonCrossSection(thetaRad[i])
            SigmaArray = np.append(SigmaArray, value*10**24)
        if C == 4:
            value = Fun.CoherentCrossSection(InitEnergy, thetaRad[i], Z, Finterpolant)
            SigmaArray = np.append(SigmaArray, value*10**24)
        if C == 5:
            value = Fun.KNpolarization(InitEnergy, thetaRad[i], 0.0) #plane phi = 0.0 degrees
            SigmaArray = np.append(SigmaArray, value*10**24)
        if C == 6:
            value = Fun.THpolarization(thetaRad[i], 0.0) #plane phi = 0.0 degrees
            SigmaArray = np.append(SigmaArray, value*10**24)
            
    return SigmaArray   
     
'''
#---------------------------------- Plotter ----------------------------------------
'''
def PlotterRandomNumbersDistributions(C, n, InitEnergy, thetaRad, SigmaArray, XR, YR):        
    fig=plt.figure(figsize=(20,10))
    fig.subplots_adjust(top=0.90, bottom=0.10, hspace=0.2, wspace=0.2)
    ax=fig.add_subplot(2,2,(1,2))
    if (C==1): 
        plt.title(r'%i X-ray photons with Klein-Nishina distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{KN}$ [barn]', fontsize=15)
    if (C==2): 
        plt.title(r'%i X-ray photons with Incoherent distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{KN}$ [barn]', fontsize=15)    
    if (C==3): 
        plt.title(r'%i X-ray photons with Thomson distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{Th}$ [barn]', fontsize=15)
    if (C==4): 
        plt.title(r'%i X-ray photons with Coherent distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{Th}$ [barn]', fontsize=15)   
    if (C==5): 
        plt.title(r'%i X-ray photons with K-N polarized distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{Th}$ [barn]', fontsize=15)   
    if (C==6): 
        plt.title(r'%i X-ray photons with Thomson polarized distribution with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
        ax.set_ylabel(r'$\sigma_{Th}$ [barn]', fontsize=15)  
    
    ax.set_xlabel(r'$\theta$ [radians]', fontsize=15)
    ax.scatter(XR,YR, s=4, color='b')   
    if C==1: tag = r'Klein-Nishina $\sigma$'
    if C==2: tag = r'Incoherent $\sigma_{inco}$'
    if C==3: tag = r'Thomson $\sigma$'
    if C==4: tag = r'Coherent $\sigma_{coh}$'
    if C==5: tag = r'Klein-Nishina polarized $\sigma_{KN}^{pol}$'
    if C==6: tag = r'Thomson polarized $\sigma_{TH}^{pol}$'   
    plt.plot(thetaRad, SigmaArray, color='r', linestyle='-', label= str(tag))
    plt.legend()
    plt.grid(True)

    ax1=fig.add_subplot(2,2,3)
    plt.hist(XR, bins=20, alpha=1, rwidth=0.8, linewidth=1, color='b', label=r'E = %f MeV' %InitEnergy)
    ax1.set_xlabel(r'$\theta$ [radians]', fontsize=15)
    ax1.set_ylabel('distribution x-axis', fontsize=15)
    plt.legend()
    plt.grid(True)

    ax2=fig.add_subplot(2,2,4)
    plt.hist(YR, bins=20, alpha=1, rwidth=0.8, linewidth=1, color='b', label=r'E = %f MeV' %InitEnergy)
    ax2.set_xlabel(r'$\sigma$ [barn]', fontsize=15)
    ax2.set_ylabel('distribution y-axis', fontsize=15)
    plt.legend()
    plt.grid(True)

    fig2=plt.figure(figsize=(20,10))
    fig2.subplots_adjust(top=0.90, bottom=0.10, hspace=0.2, wspace=0.2)

    ax3=fig2.add_subplot(111, polar=True)
    if (C==1): 
        plt.title(r'%i X-ray photons Klein-Nishina polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
    if (C==2): 
        plt.title(r'%i X-ray photons Incoherent polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
    if (C==3):
        plt.title(r'%i X-ray photons Thomson polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
    if (C==4): 
        plt.title(r'%i X-ray photons Coherent polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
    if (C==5): 
        plt.title(r'%i X-ray photons KN-polarized polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)
    if (C==6): 
        plt.title(r'%i X-ray photons TH-polarized polar distribution in [barn] with initial energy E = %f MeV' %(n, InitEnergy), fontsize=15)

    plt.plot(thetaRad, SigmaArray, linestyle='-', color = 'r', linewidth=2, label='E = %f MeV' %InitEnergy) #los angulos se meten en radianes siempre
    ax3.scatter(XR, YR, s=4, color='b')
    ax3.set_xticks(np.pi/180. * np.linspace(0,  360, 12, endpoint=False))
    ax3.set_rmin(0.0)
    plt.legend(bbox_to_anchor=(1.04, 0.65), loc=2, borderaxespad=0. , prop={'size': 11})
    plt.show()
    
    return fig


#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    print('\n')
    print('Ejecutando como programa principal \n')