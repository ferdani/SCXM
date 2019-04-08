#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:29:56 2018

@author: macbookpro

DEFINE THE READ-GET DATA FUNCTION, DEFINE THE DISTRIBUTION FUNCTIONS 
AND DEFINE THE DOSE, FLUENCE AND AUXILIAR FUNCTIONS. 
"""
import numpy as np
from scipy import interpolate
from scipy.integrate import quad
import openpyxl
import random as rd
import matplotlib.pyplot as plt

#definitions of constants:
hbar = 6.582119514*10**(-22) #Planck Constant (MeV*s)
re = 2.8179402894*10**(-13) #Electron radio (cm)
m = 0.510998928 #Electron mass (MeV/c^2)
c = 1.0 #Fundamental light constant (cm/s)

'''
------------------------------ Extract numbers for choosen Element --------------------------------------------------------
'''

def ReadGetData(Element, TableName):
    """
    Read the data into the Excel and load the values for S and F functions per each element (the choosen Element as variable)
    """
    #----------------------------- Read the data of excel of the elements------------------------------------------------
    # Should be necessary to save previously the Excel in csv and then save the csv as Excel. Only for security.        
    doc = openpyxl.load_workbook(TableName) #Load the Excel
    
    Hoja = doc.get_sheet_by_name('Hoja1') #Choose the sheet one and called Hoja

    print('\n' + 'Calculating element by element: ' + Element) 
    
    if (Element == 'Hydrogen') or (Element == 'hydrogen'):
        F = Hoja['B2':'B46']
        S = Hoja['C2':'C46']
        Z = 1.0
    if (Element == 'Carbon') or (Element == 'carbon'):
        F = Hoja['D2':'D46']
        S = Hoja['E2':'E46']
        Z = 6.0
    if (Element == 'Nitrogen') or (Element == 'nitrogen'):
        F = Hoja['F2':'F46']
        S = Hoja['G2':'G46']
        Z = 7.0
    if (Element == 'Oxygen') or (Element == 'oxygen'):
        F = Hoja['H2':'H46']
        S = Hoja['I2':'I46']
        Z = 8.0
    if (Element == 'Fluorine') or (Element == 'fluorine'):
        F = Hoja['J2':'J46']
        S = Hoja['K2':'K46']
        Z = 9.0
    if (Element == 'Sulfur') or (Element == 'sulfur'):
        F = Hoja['L2':'L46']
        S = Hoja['M2':'M46']
        Z = 16.0    
    if (Element == 'Gold') or (Element == 'gold'):
        F = Hoja['N2':'N46']
        S = Hoja['O2':'O46']
        Z = 79.0    
        
    X = Hoja['A2':'A46']
    
    Xarray = np.array([])
    Farray = np.array([])
    Sarray = np.array([])
    
    for row in X:
        for cell in row:
            Xarray = np.append(Xarray, cell.value)
    print()
    for row in F:
        for cell in row:
            Farray = np.append(Farray, cell.value)
    
    for row in S:
        for cell in row:
            Sarray = np.append(Sarray, cell.value)        
    
    '''
    ------------------------ Extract the data ---------------------------------------------------------------------------
    '''   
    #works in log space 
    XarrayLogScale = np.array([])
    SarrayLogScale = np.array([])
    FarrayLogScale = np.array([])
    
    for i in range(0,len(Xarray)):
        if ((Xarray[i] > 0.0) & (Farray[i] > 0.0) & (Sarray[i] > 0.0)):
            XarrayLogScale = np.append(XarrayLogScale, np.log(Xarray[i]))   
            FarrayLogScale = np.append(FarrayLogScale, np.log(Farray[i]))
            SarrayLogScale = np.append(SarrayLogScale, np.log(Sarray[i]))
    
    #working in small X with other interpolant cubic
    XarraySmall = np.array([])
    SarraySmall = np.array([])
    for i in range(0, 10):
        XarraySmall = np.append(XarraySmall, Xarray[i])
        SarraySmall = np.append(SarraySmall, Sarray[i])
            
    '''
    ---------------------- Build the interpolants with this data (X, S, F) -----------------------------------------------
    '''
    # S scattering function:
    #For big X one linear interpolant working in log space:
    SinterpolantBig = interpolate.interp1d(XarrayLogScale, SarrayLogScale, kind='linear')
    #For small X one cubic interpolant working in linear space:
    SinterpolantSmall = interpolate.interp1d(XarraySmall, SarraySmall, kind='cubic')        
    # F form factor
    Finterpolant = interpolate.interp1d(XarrayLogScale, FarrayLogScale, kind='linear')

    return [Z, SinterpolantBig, SinterpolantSmall, Finterpolant]


'''
------------------------------------ Interpolants for S and F functions -----------------------------
'''

def SinterpolantFunction(X, SinterpolantSmall, SinterpolantBig):
    """
    Using the data from S values per each element an interpolant is build, with this interpolant it's posible
    to extract the correction value
    """
    #using the interpolant of S function extract the correction value
    try:
        ScateringFunctionValueS = SinterpolantSmall(X)
    except:
        ScateringFunctionValueS = np.exp(SinterpolantBig(np.log(X)))

    return ScateringFunctionValueS


def FinterpolantFunction(X, Z, Finterpolant):
    """
    Using the data from F values per each element an interpolant is build, with this interpolant it's posible
    to extract the correction value
    """
    #using the interpolant of F function extract the correction value
    try:
        ScateringFunctionValueF = np.exp(Finterpolant(np.log(X)))
    except:
        ScateringFunctionValueF = Z

    return ScateringFunctionValueF

'''
------------------------------------- Beam photons--------------------------------------------------
'''

def BeamPhotons(n):
    """
    Generate the empty array of the beam pipe line with angle 0 (the angle of the beam line)
    """
    BeamPhotonsArray = [0]*n
    print('\n')
    print('Initializing the beam with ' + str(n) + ' X-rays ...')
    return BeamPhotonsArray

'''
-------------------------------------- Compton Energy ---------------------------------------------
'''

def FinalEnergy(E, theta):
    """
    Calculate the final energy of the photon with the Compton equation
    """
    FinalEnergyArray = np.array([])
    for i in range(0, len(theta)):
        Efinal = E/(1 + (E/m)*(1 - np.cos(theta[i])))
        FinalEnergyArray = np.append(FinalEnergyArray, Efinal)
    return FinalEnergyArray


'''
----------------------------------- Mean free path dependence with energy ------------------------
'''

def Xe_MeanFreePath_interpolant():
    """
    Using the data of photoelectric absorption from the NIST an interpolant is created to build the 
    mean free path for Xenon
    """
    #the path to the data
    file_path_name = './Data/Xe_NIST_PhotelAbsorb.txt'
    #extract the tow rows
    array_txt = np.loadtxt(file_path_name, usecols=(0,1), skiprows=3)
    E = array_txt[:,0] #Energy
    sigma = array_txt[:,1] #PhotoElectric mu coefficient
    #construct the interpolant lambda meen free path
    Xe_rho = 0.005851 #g/cm3
    Lamb = 1/(Xe_rho * sigma)
    Xe_MFP_interpolant = interpolate.interp1d(E, Lamb, kind='linear')

    return Xe_MFP_interpolant

def Xe_MeanFreePath(energy, Xe_MFP_interpolant):
    """
    Use the interpolant for Xenon with the mean free path and evaluate it for an input-variable energy
    """
    #with the previus interpolant make the lambdas for each energy
    Lambda = np.array([])
    for i in range(0, len(energy)):
        L = Xe_MFP_interpolant(energy[i]) #the energy in MeV
        Lambda = np.append(Lambda, L)
    return Lambda
    
'''
-------------------------------------- Ditribution Functions --------------------------------------
'''

def KleinNishinaCrossSection(InitEnergy, theta):
    """
    The differential Compton cross section called Klein-Nishina
    """
    #definition of Klein Nishina cross section. InitEnergy in MeV. Theta in rads
    P = 1.0/(1 + (InitEnergy/(m*c**2))*(1-np.cos(theta)))
    return ((re**2/2)*P**2*(P+P**-1-np.sin(theta)**2))
 
    
def ThompsonCrossSection(theta):
    """
    The differential Thomson cross section
    """
    #definition of Thompson cross section.
    return (re**2/2)*(1+np.cos(theta)**2)


def IncoherentCrossSection(InitEnergy, theta, SinterpolantSmall, SinterpolantBig):
    """
    The incoherent cross section, by definition: is the Klein-Nishina cross section corrected with S
    scattering function
    """
    #definition of the distribution Klein-Nishina corrected with S
    KN=KleinNishinaCrossSection(InitEnergy, theta)
    X=np.sin(theta/2.0)*(InitEnergy*10**6)/12398.520 #la energía en eV
    ScateringFunctionValueS = SinterpolantFunction(X, SinterpolantSmall, SinterpolantBig)
    IncoherentValue=KN*ScateringFunctionValueS
    return IncoherentValue

       
def CoherentCrossSection(InitEnergy, theta, Z, Finterpolant):
    """
    The coherent cross section, by definition: is the Klein-Nishina cross section corrected with F
    atomic form factor
    """
    #definition of the distribution Thomson corrected with F**2
    TH=ThompsonCrossSection(theta)
    X=np.sin(theta/2.0)*(InitEnergy*10**6)/12398.520 #la energía en eV
    ScateringFunctionValueF = FinterpolantFunction(X, Z, Finterpolant)
    CoherentValue=TH*ScateringFunctionValueF**2
    return CoherentValue


def KNpolarization(InitEnergy, theta, phi):
    """
    The Compton cross section with polarization, it means without phi average in Klein-Nishina
    """
    #definition of the differential cross section KleinNishina with polarization.
    FinalEnergy = InitEnergy/(1 + (InitEnergy/(m*c**2))*(1-np.cos(theta)))
    ratio = FinalEnergy/InitEnergy
    return ((re**2/2)*(ratio)**2*(ratio + ratio**-1 - 2*(np.sin(theta)**2)*(np.cos(phi)**2)))


def THpolarization(theta, phi):
    """
    The Thomson cross section with polarization
    """
    #definition of Thompson cross section with polarization.
    return (re**2)*(1-(np.sin(theta)**2)*(np.cos(phi)**2))


def IncoherentPolarization3D(InitEnergy, theta, phi, SinterpolantSmall, SinterpolantBig):
    """
    The incoherent cross section in three dimensions with the polarization
    """
    #definition of incoherent cross section with polarization in three dimensions
    THETA, PHI = np.meshgrid(theta, phi)
    sfactor = np.array([])
    for i in range(len(theta)):
        X=np.sin(theta[i]/2.0)*(InitEnergy*10**6)/12398.520 #la energía en eV
        sfactor = np.append(sfactor, SinterpolantFunction(X, SinterpolantSmall, SinterpolantBig))

    RS = KNpolarization(0.064, THETA, PHI) * sfactor
    #cartesian coordinates
    XS = RS * np.sin(THETA) * np.cos(PHI)
    YS = RS * np.sin(THETA) * np.sin(PHI)
    ZS = RS * np.cos(THETA)
    return [THETA, PHI, RS, XS, YS, ZS]


def CoherentPolarization3D(InitEnergy, theta, phi, Z, Finterpolant):
    """
    The coherent cross section in three dimensions with the polarization
    """
    #definition of incoherent cross section with polarization in three dimensions
    THETA, PHI = np.meshgrid(theta, phi)
    ffactor = np.array([])
    for i in range(len(theta)):
        X=np.sin(theta[i]/2.0)*(InitEnergy*10**6)/12398.520 #la energía en eV
        ffactor = np.append(ffactor, FinterpolantFunction(X, Z, Finterpolant)**2)

    RF = THpolarization(THETA, PHI) * ffactor
    #cartesian coordinates
    XF = RF * np.sin(THETA) * np.cos(PHI)
    YF = RF * np.sin(THETA) * np.sin(PHI)
    ZF = RF * np.cos(THETA)
    return [THETA, PHI, RF, XF, YF, ZF]


def IntegratedCrossSection(theta, dSigmaArray):
    """
    Integrated any differential cross section array without phi dependence, returns Sigma total and 
    the error in the integration
    """
    Sigma = []
    SigmaError = []
    #integramos en cos(theta)
    Xp = np.cos(theta)
    Interp = interpolate.interp1d(Xp, dSigmaArray, kind='linear')
    a, err = quad(Interp, -0.9999, 1.0)
    Sigma = np.append(Sigma, a*2*np.pi)
    SigmaError = np.append(SigmaError, err)
    return [Sigma, SigmaError]
    

def PhotonScatteringMoleculeProbabilities(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeAtoms):
    """
    Estimate the cross section per element in a molecule and the total cross-section, both for the incoherent process 
    and for the coherent process. The probabilities associated with scatterings on each element of the molecule 
    are also estimated.
    """
    SigmaIncoherentMoleculeElement = np.array([])
    SigmaCoherentMoleculeElement = np.array([])
    for ind, Z in enumerate(ArrayZ):
        #Scattering incoherente
        DiffSigmaArrayIncoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = IncoherentCrossSection(InitEnergy, thetaRad[i], ArraySinterpolantSmall[ind], ArraySinterpolantBig[ind])
            DiffSigmaArrayIncoherentElement = np.append(DiffSigmaArrayIncoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        [SigmaIncoherentElement, SigmaIncoherentErrorElement] = IntegratedCrossSection(thetaRad, DiffSigmaArrayIncoherentElement)
        SigmaIncoherentMoleculeElement = np.append(SigmaIncoherentMoleculeElement, SigmaIncoherentElement)
        #SigmaIncoherentMoleculeElement es el array con la sección eficaz incoherente integrada para cada elemento
    
        #Scattering coherente
        DiffSigmaArrayCoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = CoherentCrossSection(InitEnergy, thetaRad[i], ArrayZ[ind], ArrayFinterpolant[ind])
            DiffSigmaArrayCoherentElement = np.append(DiffSigmaArrayCoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        [SigmaCoherentElement, SigmaCoherentError] = IntegratedCrossSection(thetaRad, DiffSigmaArrayCoherentElement)
        SigmaCoherentMoleculeElement = np.append(SigmaCoherentMoleculeElement, SigmaCoherentElement)
        #SigmaCoherentMoleculeElement es el array con la sección eficaz coherente integrada para cada elemento


    #secciones eficaces totales de la molecula incoherente y coherente
    SigmaTotalIncoherentMolecule = (np.sum(SigmaIncoherentMoleculeElement * MoleculeAtoms))
    SigmaTotalCoherentMolecule = (np.sum(SigmaCoherentMoleculeElement * MoleculeAtoms))
    SigmaTotalMolecule = SigmaTotalIncoherentMolecule + SigmaTotalCoherentMolecule
    
    #Probabilities:
    ProbabilityMoleculeElement = np.array([])
    for i in range(0, len(ArrayZ)):
        pElement = ((SigmaIncoherentMoleculeElement[i] * MoleculeAtoms[i]) + (SigmaCoherentMoleculeElement[i] * MoleculeAtoms[i]))/SigmaTotalMolecule
        ProbabilityMoleculeElement = np.append(ProbabilityMoleculeElement, pElement)
    
    ProbabilityIncoherentElement = np.array([])
    ProbabilityCoherentElement = np.array([])
    for i in range(0, len(ArrayZ)):
        pElementInco = (SigmaIncoherentMoleculeElement[i] * MoleculeAtoms[i])/SigmaTotalIncoherentMolecule
        pElementCoh = (SigmaCoherentMoleculeElement[i] * MoleculeAtoms[i])/SigmaTotalCoherentMolecule
        ProbabilityIncoherentElement = np.append(ProbabilityIncoherentElement, pElementInco)
        ProbabilityCoherentElement = np.append(ProbabilityCoherentElement, pElementCoh)
    
    ProbabilityIncoherentMolecule = SigmaTotalIncoherentMolecule/SigmaTotalMolecule
    ProbabilityCoherentMolecule = SigmaTotalCoherentMolecule/SigmaTotalMolecule



    return [SigmaTotalMolecule, SigmaTotalIncoherentMolecule, SigmaTotalCoherentMolecule, \
            ProbabilityMoleculeElement, ProbabilityIncoherentElement, ProbabilityCoherentElement, \
            ProbabilityIncoherentMolecule, ProbabilityCoherentMolecule]


def PhotonScatteringMoleculeNumbers(n, ProbabilityIncoherentMolecule, ProbabilityCoherentMolecule):
    """
    Calculate the number of events n that are produced in the molecule by a distribution or another of 
    cross sections. In principle: incoherent and coherent scattering. That is, n_inco and n_coh.
    This function use the results for PhotonScatteringMoleculeProbabilities function
    """
    nMoleculeIncoherent = 0.0 #numero de fotones scattering incoherente
    nMoleculeCoherent = 0.0 #numero de fotones scattering coherente
    
    for i in range(0, len(n)):
        rdn_n = rd.uniform(0,1)
        
        if (rdn_n < ProbabilityIncoherentMolecule):
            #number for Incoherent scattering
            nMoleculeIncoherent = nMoleculeIncoherent + 1.0
        else: 
            #number for Coherent scattering
            nMoleculeCoherent = nMoleculeCoherent + 1.0   
            
    if ((nMoleculeIncoherent + nMoleculeCoherent) != len(n)):
        print('Warning in the number estimation !!!')
        
    return [int(nMoleculeIncoherent), int(nMoleculeCoherent)]


def PhotonScatteringElementNumbers(nMoleculeIncoherent, nMoleculeCoherent, ProbabilityIncoherentElement, ProbabilityCoherentElement):
    """
    Calculate the number of events n_inco_element and n_coh_element (per each element in the molecule). 
    This function use the results per molecule for PhotonScatteringMoleculeNumbers function
    """
    nMoleculeElementIncoherent = np.array([])
    nMoleculeElementCoherent = np.array([])
    
    #incoherent scattering
    for i in range(0, len(ProbabilityIncoherentElement)):
        p_incoh = ProbabilityIncoherentElement[i] * nMoleculeIncoherent
        nMoleculeElementIncoherent = np.append(nMoleculeElementIncoherent, int(round(p_incoh)))

    #coherent scattering
    for i in range(0, len(ProbabilityCoherentElement)):
        p_coh = ProbabilityCoherentElement[i] * nMoleculeCoherent
        nMoleculeElementCoherent = np.append(nMoleculeElementCoherent, int(round(p_coh)))
        
    return [nMoleculeElementIncoherent.astype(int), nMoleculeElementCoherent.astype(int)]



##############################################################################################
########################### DOSE AND FLUENCE #################################################
    
def Molar_mass(MoleculeElements, MoleculeAtoms):
    """
    Estimate molar mass of a molecule
    """
    Navogadro = 6.022140857*10**(23)
    M_element = np.array([])
    
    for i in range(0, len(MoleculeElements)):
        Element = MoleculeElements[i]
        Atoms = MoleculeAtoms[i]
        
        if (Element == 'Hydrogen') or (Element == 'hydrogen'):
            molar_mass = 1.00794 #[g/mol]
        if (Element == 'Carbon') or (Element == 'carbon'):
            molar_mass = 12.0107 #[g/mol]
        if (Element == 'Nitrogen') or (Element == 'nitrogen'):
            molar_mass = 14.0067 #[g/mol]
        if (Element == 'Oxygen') or (Element == 'oxygen'):
            molar_mass = 15.999 #[g/mol]
        if (Element == 'Fluorine') or (Element == 'fluorine'):
            molar_mass = 18.998403 #[g/mol]
        if (Element == 'Sulfur') or (Element == 'sulfur'):
            molar_mass = 32.065 #[g/mol]
        if (Element == 'Gold') or (Element == 'gold'):
            molar_mass = 196.96657 #[g/mol]
        
        M_element = np.append(M_element, Atoms*molar_mass)
        
    M_molec = sum(M_element)*(1/Navogadro) #[g/molecules]
    
    return M_molec #[g/molecules]

    
def Sigma_Mass_Omega(Sigma, Ageometric, Effi, M_molec):
    """
    Estimate the generic sigma prime (axuliar function in the Dose and Fluence evaluation)
    """
    sigma = Sigma * Ageometric * Effi / M_molec
    
    return sigma

def SigmaTotalMatter(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeAtoms):
    """
    Estimate the total cross section for the material of interest (axuliar function in the Dose and Fluence evaluation)
    """
    SigmaIncoherentMoleculeElement = np.array([])
    SigmaCoherentMoleculeElement = np.array([])
    for ind, Z in enumerate(ArrayZ):
        #Scattering incoherente
        DiffSigmaArrayIncoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = IncoherentCrossSection(InitEnergy, thetaRad[i], ArraySinterpolantSmall[ind], ArraySinterpolantBig[ind])
            DiffSigmaArrayIncoherentElement = np.append(DiffSigmaArrayIncoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        [SigmaIncoherentElement, SigmaIncoherentErrorElement] = IntegratedCrossSection(thetaRad, DiffSigmaArrayIncoherentElement)
        SigmaIncoherentMoleculeElement = np.append(SigmaIncoherentMoleculeElement, SigmaIncoherentElement)
        #SigmaIncoherentMoleculeElement es el array con la sección eficaz incoherente integrada para cada elemento
    
        #Scattering coherente
        DiffSigmaArrayCoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = CoherentCrossSection(InitEnergy, thetaRad[i], ArrayZ[ind], ArrayFinterpolant[ind])
            DiffSigmaArrayCoherentElement = np.append(DiffSigmaArrayCoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        [SigmaCoherentElement, SigmaCoherentError] = IntegratedCrossSection(thetaRad, DiffSigmaArrayCoherentElement)
        SigmaCoherentMoleculeElement = np.append(SigmaCoherentMoleculeElement, SigmaCoherentElement)
        #SigmaCoherentMoleculeElement es el array con la sección eficaz coherente integrada para cada elemento


    #secciones eficaces totales de la molecula incoherente y coherente
    SigmaTotalIncoherentMolecule = (np.sum(SigmaIncoherentMoleculeElement * MoleculeAtoms))
    SigmaTotalCoherentMolecule = (np.sum(SigmaCoherentMoleculeElement * MoleculeAtoms))
    SigmaTotalMatter = SigmaTotalIncoherentMolecule + SigmaTotalCoherentMolecule
    
    return SigmaTotalMatter


def Fluence_bkg_water(L, d_prime, d, rho_matter, rho_bkg, BkgElements, BkgAtoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter, Effi_matter, Sigma_bkg, Ageometric_bkg, Effi_bkg):
    """
    Estimate molar mass of molecule, then estimate sigma prime and then estimated the Fluence in a system: bkg-matter-bkg
    units: cm and g. All using previous auxiliar functions. With bkg only water-equivalent around the material of interest
    """
    M_molec_matter = Molar_mass(MoleculeElements, MoleculeAtoms) #[g/molecules]
    M_molec_bkg = Molar_mass(BkgElements, BkgAtoms) #[g/molecules]

    Sigma_Mass_Omega_matter = Sigma_Mass_Omega(Sigma_matter, Ageometric_matter, Effi_matter, M_molec_matter) 
    Sigma_Mass_Omega_bkg = Sigma_Mass_Omega(Sigma_bkg, Ageometric_bkg, Effi_bkg, M_molec_bkg) 

    #that is the minimun fluence
    #fluence = 25.0 * (2.0 * L * Sigma_Mass_Omega_bkg * rho_bkg) / (d**2 * d_prime**2 * (Sigma_Mass_Omega_matter * rho_matter - Sigma_Mass_Omega_bkg * rho_bkg)**2.0)
    # this is the complete equation without approximations
    fluence = 25.0 * (2.0 * Sigma_Mass_Omega_bkg * rho_bkg * (L - d) + Sigma_Mass_Omega_matter * rho_matter * d)/(d**2 * d_prime**2 * (Sigma_Mass_Omega_matter * rho_matter - Sigma_Mass_Omega_bkg * rho_bkg)**2)
    #print('minimun fluence: ', format(fluence*10**(-8), ".5E"), 'particles/mum2$')
    
    return fluence #[particles/cm2]


def Fluence_bkg_air_water(L, d_prime, d, a, rho_matter, rho_bkg_air, rho_bkg_water, Bkg_air_Elements, Bkg_water_Elements, \
                          Bkg_air_Atoms, Bkg_water_Atoms, MoleculeElements, MoleculeAtoms, Sigma_matter, Ageometric_matter, Effi_matter, \
                          Sigma_air_bkg, Sigma_water_bkg, Ageometric_air_bkg, Ageometric_water_bkg, Effi_air_bkg, Effi_water_bkg):
    """
    Estimate molar mass of molecule, then estimate sigma prime and then estimated the Fluence in a system: bkg-matter-bkg
    units: cm and g. All using auxiliar functions. With water-equivalent around the material of interest and then air around all. 
    It means, the realistic situation with cell (water-equivalent and material of interest) + bkg (air). 
    Where the parameter "a" is the size of air.
    """
    M_molec_matter = Molar_mass(MoleculeElements, MoleculeAtoms) #[g/molecules]
    M_molec_air_bkg = Molar_mass(Bkg_air_Elements, Bkg_air_Atoms) #[g/molecules]
    M_molec_water_bkg = Molar_mass(Bkg_water_Elements, Bkg_water_Atoms) #[g/molecules]
        
    Sigma_Mass_Omega_matter = Sigma_Mass_Omega(Sigma_matter, Ageometric_matter, Effi_matter, M_molec_matter) 
    Sigma_Mass_Omega_air_bkg = Sigma_Mass_Omega(Sigma_air_bkg, Ageometric_air_bkg, Effi_air_bkg, M_molec_air_bkg) 
    Sigma_Mass_Omega_water_bkg = Sigma_Mass_Omega(Sigma_water_bkg, Ageometric_water_bkg, Effi_water_bkg, M_molec_water_bkg) 
    
    # this is the complete equation without approximations
    fluence_numerador = 25.0 * (Sigma_Mass_Omega_water_bkg * rho_bkg_water * (2.0*L - d) + Sigma_Mass_Omega_air_bkg * rho_bkg_air * 4.0*a + Sigma_Mass_Omega_water_bkg * rho_bkg_water * d)
    fluence_denominador = d**2 * d_prime**2 * (Sigma_Mass_Omega_matter * rho_matter - Sigma_Mass_Omega_water_bkg * rho_bkg_water)**2
    fluence = fluence_numerador/fluence_denominador
    #print('minimun fluence: ', format(fluence*10**(-8), ".5E"), 'particles/mum2$')
    
    return fluence #[particles/cm2]



def PhotoabsorptionSigma(MoleculeElements, MoleculeAtoms, InitEnergy):
    """
    Estimate the photoabsorption sigma cross section for molecule with the NIST tables for each energy
    """
    Photoabsorption_sigma_molecule = np.array([])
    for i in range(0, len(MoleculeElements)):
        Element = MoleculeElements[i]
        Atoms = MoleculeAtoms[i]
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
    
        Photoabsorption_sigma_element = sigma_interpolant(InitEnergy)
        Photoabsorption_sigma_molecule = np.append(Photoabsorption_sigma_molecule, Photoabsorption_sigma_element*Atoms)
    
    Photoabsorption_sigma = sum(Photoabsorption_sigma_molecule)*10**(-24) #[cm2/atom]
    
    M_molec = Molar_mass(MoleculeElements, MoleculeAtoms) #[g/molecules]
    Sigma_photo = Sigma_Mass_Omega(Photoabsorption_sigma, 1.0, 1.0, M_molec)
    
    return Sigma_photo #[cm2/g·atom]

def PhotonScatteringMoleculeCrossSectionsCorrected(InitEnergy, thetaRad, ArrayZ, ArraySinterpolantSmall, ArraySinterpolantBig, ArrayFinterpolant, n, MoleculeElements, MoleculeAtoms):
    """
    Estimate the cross sections per element into molecule and the total cross section, both for the incoherent process 
    and for the coherent process with the energy and molecular mass corrections of the integral dose integral
    """
    SigmaIncoherentMoleculeElement = np.array([])
    SigmaCoherentMoleculeElement = np.array([])
    for ind, Z in enumerate(ArrayZ):
        #Scattering incoherente
        DiffSigmaArrayIncoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = IncoherentCrossSection(InitEnergy, thetaRad[i], ArraySinterpolantSmall[ind], ArraySinterpolantBig[ind])
            DiffSigmaArrayIncoherentElement = np.append(DiffSigmaArrayIncoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        SigmaIncoherentElement = IntegratedCrossSectionEnergyCorrection(thetaRad, DiffSigmaArrayIncoherentElement, InitEnergy, MoleculeElements, MoleculeAtoms)
        SigmaIncoherentMoleculeElement = np.append(SigmaIncoherentMoleculeElement, SigmaIncoherentElement)
        #SigmaIncoherentMoleculeElement es el array con la sección eficaz incoherente integrada para cada elemento
    
        #Scattering coherente
        DiffSigmaArrayCoherentElement = np.array([])
        for i in range(0, len(thetaRad)):
            value = CoherentCrossSection(InitEnergy, thetaRad[i], ArrayZ[ind], ArrayFinterpolant[ind])
            DiffSigmaArrayCoherentElement = np.append(DiffSigmaArrayCoherentElement, value*10**24)
        #integramos la sección eficaz diferencial     
        SigmaCoherentElement= IntegratedCrossSectionEnergyCorrection(thetaRad, DiffSigmaArrayCoherentElement, InitEnergy, MoleculeElements, MoleculeAtoms)
        SigmaCoherentMoleculeElement = np.append(SigmaCoherentMoleculeElement, SigmaCoherentElement)
        #SigmaCoherentMoleculeElement es el array con la sección eficaz coherente integrada para cada elemento

    #secciones eficaces totales de la molecula incoherente y coherente
    SigmaTotalIncoherentMolecule = (np.sum(SigmaIncoherentMoleculeElement * MoleculeAtoms))
    SigmaTotalCoherentMolecule = (np.sum(SigmaCoherentMoleculeElement * MoleculeAtoms))


    return [SigmaTotalIncoherentMolecule, SigmaTotalCoherentMolecule]
            
       
def IntegratedCrossSectionEnergyCorrection(theta, dSigmaArray, InitEnergy, MoleculeElements, MoleculeAtoms):
    """
    Integrated a differential cross section without phi dependence with de energy dependence 
    (see the equation in Villanueva's paper) corrected with the energy and molecular mass
    """
    Sigma = np.array([])
    SigmaError = np.array([])
    #integramos en cos(theta)
    Xp = np.cos(theta)
    Energy_difference = (1 - 1/(1 + (InitEnergy/0.510998928)*(1 - np.cos(theta))))
    Yaxis = dSigmaArray*Energy_difference
    Interp = interpolate.interp1d(Xp, Yaxis, kind='linear')
    a, err = quad(Interp, -0.9999, 1.0)
    Sigma = np.append(Sigma, a*2*np.pi)
    SigmaError = np.append(SigmaError, err)
    
    M_molec = Molar_mass(MoleculeElements, MoleculeAtoms) #[g/molecules]
    Sigma_corrected = Sigma_Mass_Omega(Sigma, 1.0, 1.0, M_molec)

    return Sigma_corrected #[cm2/g·atom]

def Dose(fluence, InitEnergy, Sigma_photo, Sigma_incoh):
    """
    Estimated the minimum dose according with Villanueva's paper. It is a local dose in [Gy]
    """
    D = fluence * InitEnergy * (Sigma_photo + Sigma_incoh) #fluence and sigmas in cm
    DoseGy = D*1.60213*10**(-13)*1000 #change units MeV to Jules and g to Kg
    #print('Dose in Gy: ', format(DoseGy, ".5E"), 'Gy')
    
    return DoseGy

def FluenceDosePlotterFeatureSize(InitEnergy, Fluence, DoseGy, FeatureSize, Sample):
    """
    The plotter for Fluence and Dose. Exist a file called PintaFluenceDosis.py to do that. 
    This function is not useful
    """
    fig=plt.figure(figsize=(20,10))
    fig.subplots_adjust(top=0.90, bottom=0.10, hspace=0.2, wspace=0.2)
    
    ax=fig.add_subplot(2,2,1)
    plt.title(str(Sample) + r' sample, initial energy E = %.3f MeV' %(InitEnergy), fontsize=15)
    ax.set_ylabel(r'Fluence [ph/$\mu$m$^{2}$]', fontsize=15)  
    ax.set_xlabel(r'Feature size [nm]', fontsize=15)
    plt.plot(FeatureSize, Fluence, color='b', linestyle='-', label='Ideal detector')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.grid(True)
    plt.legend()
    
    ax1=fig.add_subplot(2,2,3)
    plt.title(str(Sample) + r' sample, initial energy E = %.3f MeV' %(InitEnergy), fontsize=15)
    ax1.set_ylabel(r'Dose [Gy]', fontsize=15)  
    ax1.set_xlabel(r'Feature size [nm]', fontsize=15)
    plt.plot(FeatureSize, DoseGy, color='b', linestyle='-', label='Ideal detector')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.legend()
    plt.grid(True)
    
    plt.show()
    
    return fig
    
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    print('\n Ejecutando como programa principal \n')