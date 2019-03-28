#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:10:17 2019

@author: macbookpro

FLUENCE AND DOSIS FUNCTIONS
"""

import numpy as np
from scipy import interpolate
from scipy.integrate import quad
import Functions as Fun

'''
--------------------------- Parameters (cm units) ------------------------------------------------
'''

#Beam parameters:
#L = 0.1 + 0.0005 + 2.0*10**(-5) + 0.0005 + 0.1
L = 5.0*10**(-4)
d_prime = 34.0*10**(-7)

#Energy of the photons:
InitEnergy = 0.064 #MeV

#feauter size:
d = d_prime

#densities:
rho_matter = 1.70 #DNA
#rho_matter = 1.35 #Biomolecule
#rho_matter = 1.18 #PMMA
rho_bkg = 0.997 #water H2O

#effective mass coeficient:
M_molec_matter = 739.6116 * 1/(6.022*10**(23)) #[g/mol] * Navogadro = [g/molecula] for DNA
#M_molec_matter = 728.8353 * 1/(6.022*10**(23)) #for Biomolecule
#M_molec_matter = 100.1155 * 1/(6.022*10**(23)) #for PMMA
M_molec_bkg =  18.01528 * 1/(6.022*10**(23)) #[g/mol] * Navogadro = [g/molecula]

#Efficiencies and aceptances:
'''
Ageometric_matter = 0.886247 #DNA
Ageometric_bkg = 0.874459 #water H2O
Effi_matter = 0.734107  #DNA 
Effi_bkg = 0.750495 #water H2O
'''

Ageometric_matter = 1.0 #DNA
Ageometric_bkg =1.0 #water H2O
Effi_matter = 1.0  #DNA 
Effi_bkg = 1.0 #water H2O

#cross sections:
Sigma_matter = 214.3950246*10**(-24) #total cross section for DNA [cm2/atom]
#Sigma_matter = 217.867780712*10**(-24) #total cross section for biomolecule [cm2/atom]
#Sigma_matter = 30.073974187*10**(-24) #total cross section for PMMA [cm2/atom]
Sigma_bkg = 5.6225453216509171*10**(-24) #total cross section for water [cm2/atom]

#cross sections per unit mass functions
def Sigma_Mass_Omega(Sigma, Ageometric, Effi, M_molec):
    sigma = Sigma * Ageometric * Effi / M_molec
    return sigma

Sigma_Mass_Omega_matter = Sigma_Mass_Omega(Sigma_matter, Ageometric_matter, Effi_matter, M_molec_matter)
Sigma_Mass_Omega_bkg = Sigma_Mass_Omega(Sigma_bkg, Ageometric_bkg, Effi_bkg, M_molec_bkg)


'''
CALCULO DE LA SIGMA DE PHOTOABSORCION
'''

def Photoabsorption_sigma_interpolant(Element, energy):
    
    if (Element == 'Hydrogen') or (Element == 'hydrogen'):
        #the path to the data
        file_path_name = './Data/H_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Carbon') or (Element == 'carbon'):
        #the path to the data
        file_path_name = './Data/C_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Nitrogen') or (Element == 'nitrogen'):
        #the path to the data
        file_path_name = './Data/N_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Oxygen') or (Element == 'oxygen'):
        #the path to the data
        file_path_name = './Data/O_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Fluorine') or (Element == 'fluorine'):
        #the path to the data
        file_path_name = './Data/F_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Sulfur') or (Element == 'sulfur'):
        #the path to the data
        file_path_name = './Data/S_NIST_PhotelAbsorb.txt'
        
    if (Element == 'Gold') or (Element == 'gold'):
        #the path to the data
        file_path_name = './Data/Au_NIST_PhotelAbsorb.txt'
    
    #extract the tow rows
    array_txt = np.loadtxt(file_path_name, usecols=(0,1), skiprows=3)
    E = array_txt[:,0] #Energy
    sigma = array_txt[:,1] #PhotoElectric cross section [b/atom]
    sigma_interpolant = interpolate.interp1d(E, sigma, kind='linear')

    Photoabsorption_sigma = sigma_interpolant(energy)
    
    return Photoabsorption_sigma #[b/atom]


MoleculeElements = ['Hydrogen', 'Carbon', 'Nitrogen', 'Oxygen'] #DNA
MoleculeAtoms = np.array([35.5, 30.8, 1.7, 18.9]) #DNA

#MoleculeElements = ['Hydrogen', 'Carbon', 'Nitrogen', 'Oxygen', 'Sulfur'] #Biomolecule
#MoleculeAtoms = np.array([50.0, 30.0, 9.0, 10.0, 1.0]) #Biomolecule

#MoleculeElements = ['Carbon', 'Oxygen', 'Hydrogen'] #PMMA
#MoleculeAtoms = np.array([5.0, 2.0, 8.0]) #PMMA

Sigma_photo_elements = np.array([])
for i in range(0, len(MoleculeElements)):
    print('Calculating for: ' + str(MoleculeElements[i]))
    s = Photoabsorption_sigma_interpolant(MoleculeElements[i], InitEnergy)  #for DNA
    SdotAtom = s*MoleculeAtoms[i]
    Sigma_photo_elements = np.append(Sigma_photo_elements, SdotAtom)

Sigma_photo = sum(Sigma_photo_elements)*10**(-24) #[cm2/atom]

Sigma_photo_matter = Sigma_Mass_Omega(Sigma_photo, 1.0, 1.0, M_molec_matter)
print(Sigma_photo_matter)
'''
--------------------------- Fluence calcultaion ------------------------------------------------
'''

def Fluence(L, d_prime, d, Sigma_Mass_Omega_matter, Sigma_Mass_Omega_bkg):
    #that is the minimun fluence
    #fluence = 25.0 * (2.0 * L * Sigma_Mass_Omega_bkg * rho_bkg) / (d**2 * d_prime**2 * (Sigma_Mass_Omega_matter * rho_matter - Sigma_Mass_Omega_bkg * rho_bkg)**2.0)
    # this is the complete equation without approximations
    fluence = 25.0 * (2.0 * Sigma_Mass_Omega_bkg * rho_bkg * (L - d) + Sigma_Mass_Omega_matter * rho_matter * d)/(d**2 * d_prime**2 * (Sigma_Mass_Omega_matter * rho_matter - Sigma_Mass_Omega_bkg * rho_bkg)**2)
    return fluence
    

F = Fluence(L, d_prime, d, Sigma_Mass_Omega_matter, Sigma_Mass_Omega_bkg)
print('minimun fluence: ', format(F*10**(-8), ".5E"), 'particles/mum2$')

'''
-------------------------- Dose calculation ----------------------------------------------------
'''

def Incoherent_Integrated_EnergyCorrection(InitEnergy, dSigmaArray, theta):
    Sigma_incoh = np.array([])
    Sigma_incoh_Error = np.array([])
    #integramos en cos(theta)
    Xp = np.cos(theta)
    Energy_difference = (1 - 1/(1 + (InitEnergy/0.510998928)*(1 - np.cos(theta))))
    Yaxis = dSigmaArray*Energy_difference
    Interp = interpolate.interp1d(Xp, Yaxis, kind='linear')
    a, err = quad(Interp, -0.9999, 1.0)
    Sigma_incoh = np.append(Sigma_incoh, a*2*np.pi)
    Sigma_incoh_Error = np.append(Sigma_incoh, err)
    return [Sigma_incoh, Sigma_incoh_Error]

#def Incoherent_integral_for_Dose
def Dose(fluence, InitEnergy, Sigma_photo, Sigma_incoh):
    D = fluence * InitEnergy * (Sigma_photo + Sigma_incoh)
    return D

'''
thetaRad=np.linspace(0.0, 2*np.pi, num=1000, endpoint=True) #0 degrees to 360 degrees

dSigmaArray = Fun.KleinNishinaCrossSection(InitEnergy, thetaRad)

[Sigma_incoh, Sigma_incoh_Error] = Incoherent_Integrated_EnergyCorrection(InitEnergy, dSigmaArray, thetaRad)
'''
#las calculé con el main cambiando la funcion IntegratedCrossSection en functions.py
sigmaincoh = 20.320737033*10**(-24)  #DNA
#sigmaincoh = 20.6080376852*10**(-24) #Biomolecule
#sigmaincoh = 2.86335061489*10**(-24) #PMMA

Sigma_incoh_matter = Sigma_Mass_Omega(sigmaincoh, 1.0, 1.0, M_molec_matter)
print(Sigma_incoh_matter)
D = Dose(F, InitEnergy, Sigma_photo_matter, Sigma_incoh_matter)
print(D)
DoseGy = D*1.60213*10**(-13)*1000
print('Dose in Gy: ', format(DoseGy, ".5E"), 'Gy')


'''
CODIGO ABSURDO QUE SAQUÉ DEL MAIN
############################## DEPENDENCE WITH THE Initial Energy OF THE PHOTONS ################################

InitEnergyArray = np.linspace(0.001, 0.1, num=100, endpoint=True) #MeV
#InitEnergyArray = [0.064]
Fluence_E_array = np.array([])
DoseGy_E_array = np.array([])

for i in range(0, len(InitEnergyArray)):
    #Here the energy changes so we need to recalculate the cross section
    SigmaTotalMatterE = Fun.SigmaTotalMatter(InitEnergyArray[i], thetaRad, ArrayZ, ArraySinterpolantSmall, \
                                      ArraySinterpolantBig, ArrayFinterpolant, CellPhotonsArray, MoleculeAtoms)
    
    #cell matter parameters using the sigma evaluated previusly
    Sigma_matter = SigmaTotalMatterE*10**(-24) #total cross section for cell [cm2/atom]
    F = Fun.Fluence(L, d_prime, d, rho_matter, rho_bkg, BkgElements, BkgAtoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter, Effi_matter, Sigma_bkg, Ageometric_bkg, Effi_bkg)
    
    Sigma_photo = Fun.PhotoabsorptionSigma(MoleculeElements, MoleculeAtoms, InitEnergyArray[i])
    
    [SigmaTotalIncoherentMoleculeCorrectedEnergy, SigmaTotalCoherentMoleculeCorrectedEnergy] = Fun.PhotonScatteringMoleculeCrossSectionsCorrected(InitEnergyArray[i], thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeElements, MoleculeAtoms)
    DoseGy = Fun.Dose(F, InitEnergyArray[i], Sigma_photo, SigmaTotalIncoherentMoleculeCorrectedEnergy*10**(-24))

    Fluence_E_array = np.append(Fluence_E_array, F*10**(-8)) #ph/mu m^{2}
    DoseGy_E_array = np.append(DoseGy_E_array, DoseGy)
    
    sys.stdout.write("\r{0}%".format(round(((float(i)/len(InitEnergyArray))*100),1)))
    sys.stdout.flush()


data1 = np.array([InitEnergyArray, Fluence_E_array, DoseGy_E_array])  
data1 = data1.T #here you transpose your data, so to have it in two columns 
np.savetxt(tag + '_Fluence_Dose_vs_Energy.txt', data1, fmt=['%1.4f','%.5E','%.5E'], header='     Energy(MeV)                  Fluence(mum2)              Dose(Gy)')
 
fig2 = Fun.FluenceDosePlotterEnergy(Fluence_E_array, DoseGy_E_array, InitEnergyArray, tag)
'''

'''
CODIGO ABSURDO QUE SAQUÉ DEL FUNCTIONS

def FluenceDosePlotterEnergy(Fluence, DoseGy, EnergyArray, Sample):
    fig=plt.figure(figsize=(20,10))
    fig.subplots_adjust(top=0.90, bottom=0.10, hspace=0.2, wspace=0.2)
    
    ax=fig.add_subplot(2,2,1)
    plt.title(str(Sample) + r' sample', fontsize=15)
    ax.set_ylabel(r'Fluence [ph/$\mu$m$^{2}$]', fontsize=15)  
    ax.set_xlabel(r'Energy [MeV]', fontsize=15)
    plt.plot(EnergyArray, Fluence, color='b', linestyle='-', label='Ideal detector')
    ax.set_yscale('log')
    #ax.set_xscale('log')
    plt.grid(True)
    plt.legend()
    
    ax1=fig.add_subplot(2,2,3)
    plt.title(str(Sample) + r' sample', fontsize=15)
    ax1.set_ylabel(r'Dose [Gy]', fontsize=15)  
    ax1.set_xlabel(r'Energy [MeV]', fontsize=15)
    plt.plot(EnergyArray, DoseGy, color='b', marker='o', label='Ideal detector')
    ax1.set_yscale('log')
    #ax1.set_xscale('log')
    plt.legend()
    plt.grid(True)
    
    plt.show()
    
    return fig
'''

'''
############################## DEPENDENCE WITH THE d SIZE INCOHERENT ideal detector ################################
darray = np.linspace(1.0*10**(-7), 0.0001, num=100, endpoint=True) #1 nm to 1000 nm in cm case

Fluence_d_array = np.array([])
DoseGy_d_array = np.array([])

#The Init Energy is conStant so only is neccessary to evaluate one time:
SigmaTotalMatter = Fun.SigmaTotalMatter(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, \
                                        ArraySinterpolantBig, ArrayFinterpolant, CellPhotonsArray, MoleculeAtoms)
#cell matter parameters using the sigma evaluated previusly
Sigma_matter = SigmaTotalMatter*10**(-24) #total cross section for cell [cm2/atom]
    
for i in range(0, len(darray)):
    
    F = Fun.Fluence(L, darray[i], darray[i], rho_matter, rho_bkg, BkgElements, BkgAtoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter, Effi_matter, Sigma_bkg, Ageometric_bkg, Effi_bkg)
    Sigma_photo = Fun.PhotoabsorptionSigma(MoleculeElements, MoleculeAtoms, InitEnergy)
    [SigmaTotalIncoherentMoleculeCorrectedEnergy, SigmaTotalCoherentMoleculeCorrectedEnergy] = Fun.PhotonScatteringMoleculeCrossSectionsCorrected(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeElements, MoleculeAtoms)
    DoseGy = Fun.Dose(F, InitEnergy, Sigma_photo, SigmaTotalIncoherentMoleculeCorrectedEnergy*10**(-24))

    Fluence_d_array = np.append(Fluence_d_array, F*10**(-8)) #ph/mu m^{2}
    DoseGy_d_array = np.append(DoseGy_d_array, DoseGy)
    
    sys.stdout.write("\r{0}%".format(round(((float(i)/len(darray))*100),1)))
    sys.stdout.flush()
    
data = np.array([darray, Fluence_d_array, DoseGy_d_array])  
data = data.T #here you transpose your data, so to have it in two columns 
np.savetxt(tag + '_Fluence_Dose_vs_FeatureSize_IdealDetector.txt', data, fmt=['%.5E','%.5E','%.5E'], header=' FeatureSize(cm)    Fluence(mum2)   Dose(Gy)')

fig = Fun.FluenceDosePlotterFeatureSize(InitEnergy, Fluence_d_array, DoseGy_d_array, darray*10**(7), tag)
'''

'''
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    print('\n')
    print('Ejecutando como programa principal \n')
''' 