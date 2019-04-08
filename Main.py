#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:08:52 2018

@author: Daniel Fernandez Fernandez
daniel.fernandez.fernandez.94@gmail.com

THE MAIN, READ THE DATA AND RUN THE SIMULATION OF THE DETECTOR
"""
import Functions as Fun
import numpy as np
import RandomGenerator as RG
import DetectorGeometry as DG
import matplotlib.pyplot as plt
import random as rd
import sys

'''
#---------------------------- Input variables--------------------------------------------------------------------
'''

thetaRad=np.linspace(0.0, 2*np.pi, num=1000, endpoint=True) #0 degrees to 360 degrees
phiRad=np.linspace(0.0, np.pi, num=1000, endpoint=True) #0 degrees to 180 degrees

#initial energy?
InitEnergy = 0.064 #MeV 

#Distribution function: 1 for Klein-Nishina, 
#                       2 for Incoherent; Klein-Nishina with S function,
#                       3 for Thomson,
#                       4 for Coherent; Thomson with F form factor,
#                       5 for Klein-Nishina with polarization,
#                       6 for Thomson with polarization.
ArrayC = [2, 4] #un array incoherente y coherente, pero se puede correr individualemente poniendo un número dentro

#how many photons are in the beam pipe?
n = 10**4
CellPhotonsArray = Fun.BeamPhotons(n)

TableName = './Data/TablasHubbelFuncionesForma.xlsx'


ArrayZ = np.array([])
ArraySinterpolantBig = np.array([])
ArraySinterpolantSmall = np.array([])
ArrayFinterpolant = np.array([])

ResultsEfficiencies_inco = np.zeros(5)
ResultsEfficiencies_coh = np.zeros(5)

'''
-------------- choose atom or molecule ------------
'''
'''
#for example H2O - water
MoleculeElements = ['Hydrogen', 'Oxygen']
MoleculeSymbols = ['H', 'O']
MoleculeAtoms = np.array([2, 1])
print('Calculating for the stechiometric formula: ')
formula = map(''.join, zip(MoleculeSymbols, map(str, MoleculeAtoms)))
print(''.join(formula))
'''

#for example biomolecule-protein
MoleculeElements = ['Hydrogen', 'Carbon', 'Nitrogen', 'Oxygen', 'Sulfur']
MoleculeSymbols = ['H', 'C', 'N', 'O', 'S']
MoleculeAtoms = np.array([50.0, 30.0, 9.0, 10.0, 1.0])
print('Calculating for the stechiometric formula: ')
formula = map(''.join, zip(MoleculeSymbols, map(str, MoleculeAtoms)))
print(''.join(formula))
rho_matter = 1.35 #Biomolecule
tag = 'Biomolecule'

'''
#for example DNA
MoleculeElements = ['Hydrogen', 'Carbon', 'Nitrogen', 'Oxygen']
MoleculeSymbols = ['H', 'C', 'N', 'O']
MoleculeAtoms = np.array([35.5, 30.8, 1.7, 18.9])
print('Calculating for the stechiometric formula: ')
formula = map(''.join, zip(MoleculeSymbols, map(str, MoleculeAtoms)))
print(''.join(formula))
rho_matter = 1.70 #DNA
tag = 'DNA'
'''
'''
#for example PMMA
MoleculeElements = ['Carbon', 'Oxygen', 'Hydrogen']
MoleculeSymbols = ['C', 'O', 'H']
MoleculeAtoms = np.array([5.0, 2.0, 8.0])
print('Calculating for the stechiometric formula: ')
formula = map(''.join, zip(MoleculeSymbols, map(str, MoleculeAtoms)))
print(''.join(formula))
rho_matter = 1.18 #PMMA
tag = 'PMMA'
'''
'''
#for example Air
MoleculeElements = ['Nitrogen', 'Oxygen']
MoleculeSymbols = ['N', 'O']
MoleculeAtoms = np.array([1.6, 0.4])
print('Calculating for the stechiometric formula: ')
formula = map(''.join, zip(MoleculeSymbols, map(str, MoleculeAtoms)))
print(''.join(formula))
rho_matter = 0.0011838999999999999 #Air at 25 degrees temperature
tag = 'Air'
'''
'''
#for example Carbon
MoleculeElements = ['Carbon']
MoleculeSymbols = ['C']
MoleculeAtoms = np.array([1.0])
MoleculeAtoms = np.array([1])
'''

'''
#----------------------- Estimate the backgrounds -------------------------------
'''

#------------------------------change the lambda of the detector here:
#Lambda_det = 0.04285 #the mean free path of the detector material Bromuro Lantano (cm)
Lambda_det = 25.0 #cm for Xenon
#with the energy of the photons we need to correct the mean free path lambda for xenon (esto no está corregido para los fondos)
Xe_MFP_interpolant = Fun.Xe_MeanFreePath_interpolant()

#change here the parameters of the backgorudns line block (cm)
l1_air = 0.1
l1_water = 0.0005
l_cell = 2.0*10**(-5)
l2_water = l1_water
l2_air = l1_air

#change here the parameters of the detector's size
a_det = 45.0
b_det = 100.0
c_det = 5.0

#evaluate the beam photons with the line-block's interaction
#[AirLine_1_Bkg, WaterLine_1_Bkg, CellPhotonsArray, WaterLine_2_Bkg, AirLine_2_Bkg] = DG.GenerateBackgrounds(l1_air, l1_water, l_cell, l2_water, l2_air, n)
#[MeasuredPhoton_l1_air, MeasuredPhoton_l1_water, MeasuredPhoton_l2_water, MeasuredPhoton_l2_air] = DG.InteractDetectorBackgrounds(a_det, b_det, c_det, l1_air, l2_air, AirLine_1_Bkg, WaterLine_1_Bkg, WaterLine_2_Bkg, AirLine_2_Bkg, Lambda_det)

'''
CellPhotonsArray = Fun.BeamPhotons(n)
'''

'''
#----------------------------- Probabilities on cell/molecule for coherent and incoherent -------------------------
#------------------------------------ with the CellPhotonsArray of the before step --------------------------------
'''

for ind, Element in enumerate(MoleculeElements):
    [Z, SinterpolantBig, SinterpolantSmall, Finterpolant] = Fun.ReadGetData(Element, TableName)
    ArrayZ = np.append(ArrayZ, Z)
    ArraySinterpolantBig = np.append(ArraySinterpolantBig, SinterpolantBig)
    ArraySinterpolantSmall = np.append(ArraySinterpolantSmall, SinterpolantSmall)
    ArrayFinterpolant = np.append(ArrayFinterpolant, Finterpolant)


[SigmaTotalMolecule, SigmaTotalIncoherentMolecule, SigmaTotalCoherentMolecule, \
            ProbabilityMoleculeElement, ProbabilityIncoherentElement, ProbabilityCoherentElement, \
            ProbabilityIncoherentMolecule, ProbabilityCoherentMolecule] = \
Fun.PhotonScatteringMoleculeProbabilities(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, \
                                      ArraySinterpolantBig, ArrayFinterpolant, CellPhotonsArray, MoleculeAtoms)

'''
[nMoleculeIncoherent, nMoleculeCoherent] = Fun.PhotonScatteringMoleculeNumbers(CellPhotonsArray, ProbabilityIncoherentMolecule, ProbabilityCoherentMolecule)
[nMoleculeElementIncoherent, nMoleculeElementCoherent] = Fun.PhotonScatteringElementNumbers(nMoleculeIncoherent, nMoleculeCoherent, ProbabilityIncoherentElement, ProbabilityCoherentElement)
'''

'''
#----------------------------- generate the random events for incoherent and coherent---------------------------------------------------------
'''
'''
ThetaAnglesElementIncoherent = np.array([])
ThetaAnglesElementCoherent = np.array([])
DiffCrossSectionsElementIncoherent = np.array([])
DiffCrossSectionsElementCoherent = np.array([])

print('\n')
print('-----------------random numbers for incoherent scattering--------------')
for ind, Element in enumerate(MoleculeElements):
    #incoherent scattering
    [XR_ElementIncoherent, YR_ElementIncoherent] = RG.RandomNumberDistribution(ArrayC[0], nMoleculeElementIncoherent[ind], InitEnergy, thetaRad, ArrayZ[ind], ArraySinterpolantSmall[ind], ArraySinterpolantBig[ind], ArrayFinterpolant[ind])
    ThetaAnglesElementIncoherent = np.concatenate([ThetaAnglesElementIncoherent, XR_ElementIncoherent])
    DiffCrossSectionsElementIncoherent = np.concatenate([DiffCrossSectionsElementIncoherent, YR_ElementIncoherent])


print('\n')
print('----------------random numbers for coherent scattering-----------------')     
for ind, Element in enumerate(MoleculeElements):
    #coherent scattering
    [XR_ElementCoherent, YR_ElementCoherent] = RG.RandomNumberDistribution(ArrayC[1], nMoleculeElementCoherent[ind], InitEnergy, thetaRad, ArrayZ[ind], ArraySinterpolantSmall[ind], ArraySinterpolantBig[ind], ArrayFinterpolant[ind])
    ThetaAnglesElementCoherent = np.concatenate([ThetaAnglesElementCoherent, XR_ElementCoherent])
    DiffCrossSectionsElementCoherent = np.concatenate([DiffCrossSectionsElementCoherent, YR_ElementCoherent])
'''

'''
#------------------------ plot the cross sections for incoherent and coherent and with random numbers ------------------------------------
'''

#fig = RG.PlotterRandomNumbersDistributions(C, n, InitEnergy, thetaRad, SigmaArray, XR, YR)
 

'''
#------------------------ Interaction in the detector and its geometry ----------------------------------------
'''
'''
#water block dimensions (cm units):
a_water = 0.001
b_water = 0.001

#air block dimensions (cm units):
a_air = c_det 
b_air = b_det

#--------------------------------------------- incoherente --------------------------------------------
UpperPhotons_XR_inco = DG.UpperPhotons(ThetaAnglesElementIncoherent) #angles distribution
UpperPhotons_FinalEnergy_inco = Fun.FinalEnergy(InitEnergy, UpperPhotons_XR_inco) #Final energy of the photons, the photons against the cell have got the initial energy 0.064
Lambda_inco = Fun.Xe_MeanFreePath(UpperPhotons_FinalEnergy_inco, Xe_MFP_interpolant)

Upper_NoWaterInteract_XR_inco = DG.WaterBlock(a_water, b_water, UpperPhotons_XR_inco)
Upper_NoAirInteract_XR_inco = DG.AirBlock(a_air, b_air, Upper_NoWaterInteract_XR_inco)

[permitted_XR_inco, MeasuredPhoton_XR_inco] = DG.DetectorBox(a_det, b_det, c_det, Upper_NoAirInteract_XR_inco, Lambda_inco)
MeasuredPhoton_FinalEnergy_inco = Fun.FinalEnergy(InitEnergy, MeasuredPhoton_XR_inco)

[Eff_water_inco, Eff_air_inco, A_geo_inco, Eff_det_inco, Eff_tot_inco] = DG.Efficiencies(UpperPhotons_XR_inco, Upper_NoWaterInteract_XR_inco, Upper_NoAirInteract_XR_inco, permitted_XR_inco, MeasuredPhoton_XR_inco)
ResultsEfficiencies_inco = np.vstack([ResultsEfficiencies_inco, [Eff_water_inco, Eff_air_inco, A_geo_inco, Eff_det_inco, Eff_tot_inco]])

fig0 = DG.PlotterEfficiencyDistributions(CellPhotonsArray, InitEnergy, UpperPhotons_XR_inco, Upper_NoAirInteract_XR_inco, permitted_XR_inco, MeasuredPhoton_XR_inco)
    
#delete the 0 rows due to
ResultsEfficiencies_inco = ResultsEfficiencies_inco[~np.all(ResultsEfficiencies_inco == 0, axis=1)]


#--------------------------------------------- coherente --------------------------------------------
UpperPhotons_XR_coh = DG.UpperPhotons(ThetaAnglesElementCoherent)
UpperPhotons_FinalEnergy_coh = Fun.FinalEnergy(InitEnergy, UpperPhotons_XR_coh)
Lambda_coh = Fun.Xe_MeanFreePath(UpperPhotons_FinalEnergy_coh, Xe_MFP_interpolant)

Upper_NoWaterInteract_XR_coh = DG.WaterBlock(a_water, b_water, UpperPhotons_XR_coh)
Upper_NoAirInteract_XR_coh = DG.AirBlock(a_air, b_air, Upper_NoWaterInteract_XR_coh)

[permitted_XR_coh, MeasuredPhoton_XR_coh] = DG.DetectorBox(a_det, b_det, c_det, Upper_NoAirInteract_XR_coh, Lambda_coh)
MeasuredPhoton_FinalEnergy_coh = Fun.FinalEnergy(InitEnergy, MeasuredPhoton_XR_coh)

[Eff_water_coh, Eff_air_coh, A_geo_coh, Eff_det_coh, Eff_tot_coh] = DG.Efficiencies(UpperPhotons_XR_coh, Upper_NoWaterInteract_XR_coh, Upper_NoAirInteract_XR_coh, permitted_XR_coh, MeasuredPhoton_XR_coh)
ResultsEfficiencies_coh = np.vstack([ResultsEfficiencies_coh, [Eff_water_coh, Eff_air_coh, A_geo_coh, Eff_det_coh, Eff_tot_coh]])

#fig1 = DG.PlotterEfficiencyDistributions(CellPhotonsArray, InitEnergy, UpperPhotons_XR_coh, Upper_NoAirInteract_XR_coh, permitted_XR_coh, MeasuredPhoton_XR_coh)
    
#delete the 0 rows
ResultsEfficiencies_coh = ResultsEfficiencies_coh[~np.all(ResultsEfficiencies_coh == 0, axis=1)]
'''

'''
#----------------------------- Fluence and Dose -----------------------------------------------------------------------
'''

#voxel parameters: ¡¡¡Son distintos de los parametros de la celula y agua para evaluar las eficiencias, realmente da igual !!!
#L = 5.0*10**(-4)
L = 10.0*10**(-4)
d_prime = 34.0*10**(-7)

#feauture size:
d = d_prime

#air size:
a = 0.1

#bkg water parameters:
rho_bkg_water = 0.997 #water H2O
Bkg_water_Elements = ['Hydrogen', 'Oxygen']
Bkg_water_Atoms = [2.0, 1.0]
Sigma_water_bkg = 5.6225453216509171*10**(-24) #total cross section for water [cm2/atom] at 0.064 MeV initial energy

#bkg air parameters
rho_bkg_air = 0.0011838999999999999 #air NO
Bkg_air_Elements = ['Nitrogen', 'Oxygen']
Bkg_air_Atoms = [1.6, 0.4]
Sigma_air_bkg = 8.10279297127*10**(-24) #total cross section for air [cm2/atom] at 0.064 MeV initial energy

############################## DEPENDENCE WITH THE d SIZE INCOHERENT FOR DIFFERENT DETECTORS ################################
#Ideal Detector // Detector 1: 100*45*5 // Detector 2: 100*45*1 // Detector 3: 50*25*5

Ageometric_matter_array = np.array([1.0, 0.95, 0.88])
Effi_matter_array       = np.array([1.0, 0.91, 0.74])
Ageometric_bkg_array    = np.array([1.0, 0.95, 0.88]) #The water and air have got the same efficiency and aceptance
Effi_bkg_array          = np.array([1.0, 0.91, 0.74])

# 'IdealDetector', '100-45-5', '50-25-5'
#DetectorTag = np.array(['Ideal', 'Optimistic', 'Realistic'])
DetectorTag = np.array(['Ideal_bkg_complete', 'Optimistic_bkg_complete', 'Realistic_bkg_complete'])

darray = np.linspace(1.0*10**(-7), 0.0001, num=1000, endpoint=True) #1 nm to 1000 nm in cm case
        
#The Init Energy is conStant so only is neccessary to evaluate one time:
SigmaTotalMatter = Fun.SigmaTotalMatter(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, \
                                        ArraySinterpolantBig, ArrayFinterpolant, CellPhotonsArray, MoleculeAtoms)
#cell matter parameters using the sigma evaluated previusly
Sigma_matter = SigmaTotalMatter*10**(-24) #total cross section for cell [cm2/atom]

for j in range(0, len(Ageometric_matter_array)):   
    if (j==0):    
        print('\n')
        print('Doing the loop of Fluence and Dose ...')
        print('\n')    
   
    data = []
    Fluence_d_array = np.array([])
    DoseGy_d_array = np.array([])
    
    for i in range(0, len(darray)):    
        #F = Fun.Fluence_bkg_water(L, darray[i], darray[i], rho_matter, rho_bkg_water, Bkg_water_Elements, Bkg_water_Atoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter_array[j], Effi_matter_array[j], Sigma_water_bkg, Ageometric_bkg_array[j], Effi_bkg_array[j])
        F = Fun.Fluence_bkg_air_water(L, darray[i], darray[i], a, rho_matter, rho_bkg_air, rho_bkg_water, Bkg_air_Elements, Bkg_water_Elements, \
                          Bkg_air_Atoms, Bkg_water_Atoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter_array[j], Effi_matter_array[j], \
                          Sigma_air_bkg, Sigma_water_bkg, Ageometric_bkg_array[j], Ageometric_bkg_array[j], Effi_bkg_array[j], Effi_bkg_array[j])

        Sigma_photo = Fun.PhotoabsorptionSigma(MoleculeElements, MoleculeAtoms, InitEnergy)
        [SigmaTotalIncoherentMoleculeCorrectedEnergy, SigmaTotalCoherentMoleculeCorrectedEnergy] = Fun.PhotonScatteringMoleculeCrossSectionsCorrected(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeElements, MoleculeAtoms)
        DoseGy = Fun.Dose(F, InitEnergy, Sigma_photo, SigmaTotalIncoherentMoleculeCorrectedEnergy*10**(-24))
    
        Fluence_d_array = np.append(Fluence_d_array, F*10**(-8)) #ph/mu m^{2}
        DoseGy_d_array = np.append(DoseGy_d_array, DoseGy)
        
        sys.stdout.write("\r{0}%".format(round(((float(i)/len(darray))*100),1)))
        sys.stdout.flush()
        
    print('\n' + 'index of sample: ', j)
    
    data = np.array([darray, Fluence_d_array, DoseGy_d_array])  
    data = data.T #here you transpose your data, so to have it in two columns 
    np.savetxt('./OutputFiles/' + tag + '_Fluence_Dose_vs_FeatureSize_' + DetectorTag[j] + '.txt', data, fmt= ['%.5E','%.5E','%.5E'] , header=' FeatureSize(cm)    Fluence(mum2)   Dose(Gy)')
   
    
#fig = Fun.FluenceDosePlotterFeatureSize(InitEnergy, Fluence_d_array, DoseGy_d_array, darray*10**(7), tag)


'''
#------------------------------- run the code independently for example to calculate random numbers -------------------------------------------------
'''
'''
thetaRandomArray = np.array([])
phiRandomArray = np.array([])
THETA_RANDOM_ARRAY_DNA = np.array([])
PHI_RANDOM_ARRAY_DNA = np.array([])
SIGMA_RANDOM_ARRAY_DNA = np.array([])


for k in range(0, 10): #haciendo esto añadimos ordenes de magnitud a los eventos que buscamos
    print('index: ', k)
    for j in range(0, len(nMoleculeElementIncoherent)):
        for i in range(0, nMoleculeElementIncoherent[j]):
            thetaRandom = rd.uniform(0,1)*np.pi*2.0
            thetaRandomArray = np.append(thetaRandomArray, thetaRandom)
            phiRandom = rd.uniform(0,1)*np.pi
            phiRandomArray = np.append(phiRandomArray, phiRandom)
        [THETA, PHI, RS, XS, YS, ZS] = Fun.IncoherentPolarization3D(InitEnergy, thetaRandomArray, phiRandomArray, ArraySinterpolantSmall[j], ArraySinterpolantBig[j])
        THETA_RANDOM_ARRAY_DNA = np.append(THETA_RANDOM_ARRAY_DNA, THETA[0])
        p = [i[0] for i in PHI] #extrayendo los valores distintos, esto se debe al meshgrid de dentro de la funcion Incoherent3D
        PHI_RANDOM_ARRAY_DNA = np.append(PHI_RANDOM_ARRAY_DNA, p)
        SIGMA_RANDOM_ARRAY_DNA = np.append(SIGMA_RANDOM_ARRAY_DNA, RS)
        
        THETA = []
        PHI = []
        p = []
        RS = []
        XS = []
        YS = []
        ZS = []
        thetaRandomArray = []
        phiRandomArray = []

data = np.array([THETA_RANDOM_ARRAY_DNA, PHI_RANDOM_ARRAY_DNA])  
data = data.T #here you transpose your data, so to have it in two columns 
np.savetxt('Random_3D_polarized_incoherent_distribution.txt', data, fmt=['%1.18f','%1.18f'], header='   theta                  phi')
'''

'''
SigmaInco = np.array([])
    
for i in range(0, len(thetaRad)):
    SigmaInco = np.append(SigmaInco, Fun.IncoherentCrossSection(InitEnergy, thetaRad[i]))
'''