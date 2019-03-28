#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:48:56 2018

@author: macbookpro

DETECTOR GEOMETRY, BACKGROUNDS, INTERACTIONS ON THE DETECTOR, GEOMETRIC ACEPTANCE AND EFFICIENCIES 
"""
import math
import numpy as np
import random as rd
import sys
import matplotlib.pyplot as plt

'''
------------------------------------------ Detector, Cell and Air origin coordinates --------------------------------
'''
#by defect,detector origin coordinates (cm units):
X_det = 0.0
Y_det = 0.0
Z_det = 0.0
#by defect, cell origin coordinates (cm units):
X_cell = 0.0
Y_cell = 0.0
Z_cell = 0.0

'''
--------------------------------------- Auxiliar Functions ---------------------------
'''

#------------------- probability law  ---------------------
def ProbabilityInteract(Path, Lambda):
    return (1 - np.exp(-Path/Lambda))


#----------------- pull the upper photons ----------------
def UpperPhotons(XR):
    UpperPhotons_XR = np.array([])
    for i in range(0,len(XR)):
        if (XR[i] > 0.0) & (XR[i] < math.pi):
            UpperPhotons_XR = np.append(UpperPhotons_XR, XR[i])
    return UpperPhotons_XR


#----------------- Correction for Detector's corners, Airblock's corners, Waterblock's corners -----------
def Corners(d, large, hole, theta1):
    thetacritic1 = np.arctan((hole + d)/(large/2.0))
    thetacritic2 = np.pi - np.arctan((hole + d)/(large/2.0))
    
    if (theta1 < thetacritic1):
        corner = (large/2.0)*(np.tan(theta1)) - hole
    if (theta1 > thetacritic2):
        corner = (large/2.0)*(np.tan(np.pi - theta1)) - hole
    else:
        corner = d

    return corner #the real distance of the corner


'''
------------------------------------------ Detector Enviroment -----------------------------------------------------
'''    
#-------------------probability to colision in the water block before the detector----------------
def WaterBlock(a_water, b_water, IntroducedPhotonsWater_XR):
    #calculate the interaction path of the photons into the water block
    r_water = np.array([]) #the path of photons in the water
    for i in range(0,len(IntroducedPhotonsWater_XR)):
        corner_water = Corners(a_water/2.0, b_water/2.0, 0.0, IntroducedPhotonsWater_XR[i])
        #Calculation of the delta-r path into the water block aka r_water
        if (IntroducedPhotonsWater_XR[i] < np.pi/2.0):
            r_water = np.append(r_water, corner_water/math.sin(IntroducedPhotonsWater_XR[i]))
        else:
            r_water = np.append(r_water, corner_water/math.sin(np.pi - IntroducedPhotonsWater_XR[i]))
            
    #interaction into the water: Probability_interaction = 1 - exp(-r_water/MeanFreePath)
    Lambda_water = 5.66671 #cm updated by nist xcom
    Pinteract_water = np.array([]) #probability of interact in the water
    print('\n')
    print('Interaction with the water block...')
    for i in range(0,len(IntroducedPhotonsWater_XR)):
        P_water = ProbabilityInteract(r_water[i], Lambda_water)
        Pinteract_water = np.append(Pinteract_water, P_water)
       
    #estimate the interactions in the water block and take the survival photons
    print('\n')
    print('Calculate the interactions in the water block...')
    Upper_NoWaterInteract_XR = np.array([]) #survival photons after interact
    for i in range(0,len(Pinteract_water)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_water[i] < rdn_water):
            Upper_NoWaterInteract_XR = np.append(Upper_NoWaterInteract_XR, IntroducedPhotonsWater_XR[i])
        
    #return the upper photons and the survival photons after te water block
    return Upper_NoWaterInteract_XR

#-------------------probability to colision in the air block before the detector----------------
def AirBlock(a_air, b_air, IntroducedPhotonsAir_XR):
    #calculate the interaction path of the photons into the air block
    r_air = np.array([]) #the path of photons in the air
    for i in range(0,len(IntroducedPhotonsAir_XR)):
        corner_air = Corners(a_air/2.0, b_air/2.0, 0.0, IntroducedPhotonsAir_XR[i])
        #Calculation of the delta-r path into the air block aka r_air
        if (IntroducedPhotonsAir_XR[i] < np.pi/2.0):
            r_air = np.append(r_air, corner_air/math.sin(IntroducedPhotonsAir_XR[i]))
        else:
            r_air = np.append(r_air, corner_air/math.sin(np.pi - IntroducedPhotonsAir_XR[i]))
        
    #interaction into the Air: Probability_interaction = 1 - exp(-r_air/MeanFreePath)
    Lambda_air = 4124.77 #cm updated by nist xcom
    Pinteract_air = np.array([]) #probability of interact in the air
    print('\n')
    print('Interaction with the air block...')
    for i in range(0,len(IntroducedPhotonsAir_XR)):
        P_air = ProbabilityInteract(r_air[i], Lambda_air)
        Pinteract_air = np.append(Pinteract_air, P_air)
       
    #estimate the interactions in the Air Block and take the survival photons
    print('\n')
    print('Calculate the interactions in the air block...')
    Upper_NoAirInteract_XR = np.array([]) #survival photons after interact
    for i in range(0,len(Pinteract_air)):
        rdn_air = rd.uniform(0,1)
        if (Pinteract_air[i] < rdn_air):
            Upper_NoAirInteract_XR = np.append(Upper_NoAirInteract_XR, IntroducedPhotonsAir_XR[i])
            
    #return the upper photons and the survival photons after te Air Block
    return Upper_NoAirInteract_XR

#---------------------------Detector box description and interaction ----------------------------
def DetectorBox(a_det, b_det, c_det, IntroducedPhotonsDetector_XR, Lambda_det): #dimensions (cm) of detector box a, b, c.
    #----------------------------Correct with X_det... and X_cell----------
    #----bla bla bla

    #forbiden angle zone detection with the before geometry (radians unit)
    theta_forbiden_A1 = np.arctan(c_det/(b_det/2.0))
    theta_forbiden_A2 = math.pi/2.0 + np.arctan((b_det/2.0)/c_det)
    
    #the permited photons that maybe will be detected (all are contained in the detection zone)
    permitted_XR = np.array([]) #los que entran en la region de interaccion con la diagonales y pueden ser detectados
    for i in range(0,len(IntroducedPhotonsDetector_XR)):
        if (IntroducedPhotonsDetector_XR[i] >= theta_forbiden_A1) & (IntroducedPhotonsDetector_XR[i] <= theta_forbiden_A2):
            permitted_XR = np.append(permitted_XR, IntroducedPhotonsDetector_XR[i])

    #Calculation of the parameters for each permitted photon
    r_det = np.array([])
    for i in range(0,len(permitted_XR)):
        #calculation of the corner with the previus function
        corner_det = Corners(a_det, b_det, c_det, permitted_XR[i])
        #Calculation of the delta-r path into the detector aka r_det
        if (permitted_XR[i] < np.pi/2.0):
            r_det = np.append(r_det, corner_det/math.sin(permitted_XR[i]))
        else:
            r_det = np.append(r_det, corner_det/math.sin(np.pi - permitted_XR[i]))
            
    #interaction into the gas/matter of the detector: Probability_interaction = 1 - exp(-r_det/FreePath)
    Pinteract_det = np.array([]) #probability of interactions with the detector
    MeasuredPhoton_XR = np.array([]) #measured photons with the detector (the final survival photons)
    print('\n')
    print('Measuring the photon events with the detector...')

    for i in range(0,len(permitted_XR)):
        if isinstance(Lambda_det, float) == True:
            P_det = ProbabilityInteract(r_det[i], Lambda_det)
        else:
            P_det = ProbabilityInteract(r_det[i], Lambda_det[i])
        Pinteract_det = np.append(Pinteract_det, P_det)

    for i in range(0,len(Pinteract_det)):
        rdn_det = rd.uniform(0,1)
        if (Pinteract_det[i] > rdn_det):
            MeasuredPhoton_XR = np.append(MeasuredPhoton_XR, permitted_XR[i])
        
        sys.stdout.write("\r{0}%".format(round(((float(i)/len(Pinteract_det))*100),1)))
        sys.stdout.flush()
            
    #return the final photons, the detected photons with its energy
    return [permitted_XR, MeasuredPhoton_XR]

'''
#----------------------------------------------- Backgrounds -------------------------------------
'''

def GenerateBackgrounds(l1_air, l1_water, l_cell, l2_water, l2_air, n):
    """
    probability to colision in the air line after the beam pipe, then in the water, then in the cell
    and again into water and finally air
    """
    
    '''
    -----------------------------------------
    ----- with reweighted probabilities -----
    -----------------------------------------
    '''
    
    AirLine_1_Bkg = np.array([]) #fotones que han interctuado en la linea de aire 1 y serán bkg
    WaterLine_1_Bkg = np.array([]) #fotones que han interctuado en la linea de agua 1 y serán bkg
    CellPhotonsArray = np.array([]) #los fotones que han interactuado en la célula que pasarán al MC
    WaterLine_2_Bkg = np.array([]) #fotones que han interctuado en la linea de agua 2 y serán bkg
    AirLine_2_Bkg = np.array([]) #fotones que han interctuado en la linea de aire 2 y serán bkg
    
    #Mean free paths:
    Lambda_air = 4124.77 #cm updated by nist xcom
    Lambda_water = 5.66671 #cm updated by nist xcom
    Lambda_cell = 3.3698418318647483 #cm evaluated by Dani for DNA
    #Lambda_cell = 4.11492449 #cm evaluated by Dani for Biomolecule
    #Lambda_cell = 4.684764 #cm evaluated by Dani for PMMA
    #Lambda_cell = 5.66671 #cm when the cell is water
    
    #Probabilities:
    P_air_1_old = ProbabilityInteract(l1_air, Lambda_air)
    P_water_1_old = ProbabilityInteract(l1_water, Lambda_water)
    P_cell_old = ProbabilityInteract(l_cell, Lambda_cell)
    P_air_2_old = ProbabilityInteract(l2_air, Lambda_air)
    P_water_2_old = ProbabilityInteract(l2_water, Lambda_water)
    
    #Weighted the probabilties
    Ptotal = P_air_1_old + P_water_1_old + P_cell_old + P_air_2_old + P_water_2_old
    
    P_air_1_new = P_air_1_old/Ptotal
    P_water_1_new = P_water_1_old/Ptotal
    P_cell_new = P_cell_old/Ptotal
    P_air_2_new = P_air_2_old/Ptotal
    P_water_2_new = P_water_2_old/Ptotal
    
    con = 1
    if (con == 1):
        print('aire1:', P_air_1_new)
        print('aire2:', P_air_2_new)
        print('agua1:', P_water_1_new)
        print('agua2:', P_water_2_new)
        print('cell:', P_cell_new)
        print('ptotal:', Ptotal)
        con = 2
    
    
    #-------------------------------------------------------------------------------------------
    #El beam pipe entra en aire lo primero:
    
    Pinteract_l1_air = np.array([]) #probability of interact in the line 1 air
    print('\n')
    print('Interaction with the air line 1, calculating backgrounds...')
    for i in range(0, n):
        Pinteract_l1_air = np.append(Pinteract_l1_air, P_air_1_new)
          
    NoInteract_l1_air = np.array([]) #fotones supervivientes que continuarán al agua
    for i in range(0,len(Pinteract_l1_air)):
        rdn_air = rd.uniform(0,1)
        if (Pinteract_l1_air[i] < rdn_air):
            NoInteract_l1_air = np.append(NoInteract_l1_air, 0.0)
        else:
            AirLine_1_Bkg = np.append(AirLine_1_Bkg, 0.0)
    
    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(AirLine_1_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        AirLine_1_Bkg[i] = rdn_air_angle

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del aire pasan al agua:
    
    Pinteract_l1_water = np.array([]) #probability of interact in the line 1 water
    print('\n')
    print('Interaction with the water line 1, calculating backgrounds...')
    for i in range(0,len(NoInteract_l1_air)):
        Pinteract_l1_water = np.append(Pinteract_l1_water, P_water_1_new)
     
    NoInteract_l1_water = np.array([]) #fotones supervivientes que continuarán a la celula
    for i in range(0,len(Pinteract_l1_water)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l1_water[i] < rdn_water):
            NoInteract_l1_water = np.append(NoInteract_l1_water, NoInteract_l1_air[i])
        else:
            WaterLine_1_Bkg = np.append(WaterLine_1_Bkg, NoInteract_l1_air[i])
                
    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(WaterLine_1_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        WaterLine_1_Bkg[i] = rdn_air_angle 

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del agua pasan a la célula:
    

    Pinteract_l_cell = np.array([]) #probability of interact in the cell
    print('\n')
    print('Interaction with the cell, calculating backgrounds...')
    for i in range(0,len(NoInteract_l1_water)):
        Pinteract_l_cell = np.append(Pinteract_l_cell, P_cell_new)
       
    NoInteract_l_cell = np.array([]) #fotones supervivientes que continuarán al agua siguiente
    for i in range(0,len(Pinteract_l_cell)):
        rdn_cell = rd.uniform(0,1)
        if (Pinteract_l_cell[i] < rdn_cell):
            NoInteract_l_cell = np.append(NoInteract_l_cell, NoInteract_l1_water[i])
        else:
            CellPhotonsArray = np.append(CellPhotonsArray, NoInteract_l1_water[i])
    
    #-------------------------------------------------------------------------------------------
    #Los cell photon array son los fotones que pasarán al MC simulation y son nuestra señal
    #-------------------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes de la célula pasan al agua 
    
    Lambda_water = 5.66671 #cm updated by nist xcom
    Pinteract_l2_water = np.array([]) #probability of interact in the line 2 water
    print('\n')
    print('Interaction with the water line 2, calculating backgrounds...')
    for i in range(0,len(NoInteract_l_cell)):
        Pinteract_l2_water = np.append(Pinteract_l2_water, P_water_2_new)

    NoInteract_l2_water = np.array([]) #fotones supervivientes de la célula que continuarán al aire siguiente
    for i in range(0,len(Pinteract_l2_water)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l2_water[i] < rdn_water):
            NoInteract_l2_water = np.append(NoInteract_l2_water, NoInteract_l_cell[i])
        else:
            WaterLine_2_Bkg = np.append(WaterLine_2_Bkg, NoInteract_l_cell[i])

    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(WaterLine_2_Bkg)):
        rdn_water_angle = rd.uniform(0,1)*np.pi*2
        WaterLine_2_Bkg[i] = rdn_water_angle

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del agua pasan al aire
    
    Lambda_air = 4124.77 #cm updated by nist xcom
    Pinteract_l2_air = np.array([]) #probability of interact in the line 2 air
    print('\n')
    print('Interaction with the air line 2, calculating backgrounds...')
    
    for i in range(0,len(NoInteract_l2_water)):
        Pinteract_l2_air = np.append(Pinteract_l2_air, P_air_2_new)

    NoInteract_l2_air = np.array([]) #fotones supervivientes de la célula que continuarán al aire siguiente
    for i in range(0,len(Pinteract_l2_air)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l2_air[i] < rdn_water):
            NoInteract_l2_air = np.append(NoInteract_l2_air, NoInteract_l2_water[i])
        else:
            AirLine_2_Bkg = np.append(AirLine_2_Bkg, NoInteract_l2_water[i])

    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(AirLine_2_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        AirLine_2_Bkg[i] = rdn_air_angle

    return [AirLine_1_Bkg, WaterLine_1_Bkg, CellPhotonsArray, WaterLine_2_Bkg, AirLine_2_Bkg]



'''
    AirLine_1_Bkg = np.array([]) #fotones que han interctuado en la linea de aire 1 y serán bkg
    WaterLine_1_Bkg = np.array([]) #fotones que han interctuado en la linea de agua 1 y serán bkg
    CellPhotonsArray = np.array([]) #los fotones que han interactuado en la célula que pasarán al MC
    WaterLine_2_Bkg = np.array([]) #fotones que han interctuado en la linea de agua 2 y serán bkg
    AirLine_2_Bkg = np.array([]) #fotones que han interctuado en la linea de aire 2 y serán bkg
    
    #-------------------------------------------------------------------------------------------
    #El beam pipe entra en aire lo primero:
    
    Lambda_air = 4124.77 #cm updated by nist xcom
    Pinteract_l1_air = np.array([]) #probability of interact in the line 1 air
    print('\n')
    print('Interaction with the air line 1, calculating backgrounds...')
    for i in range(0, n):
        P_air = ProbabilityInteract(l1_air, Lambda_air)
        Pinteract_l1_air = np.append(Pinteract_l1_air, P_air)
          
    NoInteract_l1_air = np.array([]) #fotones supervivientes que continuarán al agua
    for i in range(0,len(Pinteract_l1_air)):
        rdn_air = rd.uniform(0,1)
        if (Pinteract_l1_air[i] < rdn_air):
            NoInteract_l1_air = np.append(NoInteract_l1_air, 0.0)
        else:
            AirLine_1_Bkg = np.append(AirLine_1_Bkg, 0.0)
    
    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(AirLine_1_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        AirLine_1_Bkg[i] = rdn_air_angle

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del aire pasan al agua:
    
    Lambda_water = 5.66671 #cm updated by nist xcom
    Pinteract_l1_water = np.array([]) #probability of interact in the line 1 water
    print('\n')
    print('Interaction with the water line 1, calculating backgrounds...')
    for i in range(0,len(NoInteract_l1_air)):
        P_water = ProbabilityInteract(l1_water, Lambda_water)
        Pinteract_l1_water = np.append(Pinteract_l1_water, P_water)
     
    NoInteract_l1_water = np.array([]) #fotones supervivientes que continuarán a la celula
    for i in range(0,len(Pinteract_l1_water)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l1_water[i] < rdn_water):
            NoInteract_l1_water = np.append(NoInteract_l1_water, NoInteract_l1_air[i])
        else:
            WaterLine_1_Bkg = np.append(WaterLine_1_Bkg, NoInteract_l1_air[i])
                
    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(WaterLine_1_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        WaterLine_1_Bkg[i] = rdn_air_angle 

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del agua pasan a la célula:
    
    Lambda_cell = 4.243885 #cm evaluated by Dani
    Pinteract_l_cell = np.array([]) #probability of interact in the cell
    print('\n')
    print('Interaction with the cell, calculating backgrounds...')
    for i in range(0,len(NoInteract_l1_water)):
        P_cell = ProbabilityInteract(l_cell, Lambda_cell)
        Pinteract_l_cell = np.append(Pinteract_l_cell, P_cell)
       
    NoInteract_l_cell = np.array([]) #fotones supervivientes que continuarán al agua siguiente
    for i in range(0,len(Pinteract_l_cell)):
        rdn_cell = rd.uniform(0,1)
        if (Pinteract_l_cell[i] < rdn_cell):
            NoInteract_l_cell = np.append(NoInteract_l_cell, NoInteract_l1_water[i])
        else:
            CellPhotonsArray = np.append(CellPhotonsArray, NoInteract_l1_water[i])
    
    #-------------------------------------------------------------------------------------------
    #Los cell photon array son los fotones que pasarán al MC simulation y son nuestra señal
    #-------------------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes de la célula pasan al agua 
    
    Lambda_water = 5.66671 #cm updated by nist xcom
    Pinteract_l2_water = np.array([]) #probability of interact in the line 2 water
    print('\n')
    print('Interaction with the water line 2, calculating backgrounds...')
    for i in range(0,len(NoInteract_l_cell)):
        P_water = ProbabilityInteract(l2_water, Lambda_water)
        Pinteract_l2_water = np.append(Pinteract_l2_water, P_water)

    NoInteract_l2_water = np.array([]) #fotones supervivientes de la célula que continuarán al aire siguiente
    for i in range(0,len(Pinteract_l2_water)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l2_water[i] < rdn_water):
            NoInteract_l2_water = np.append(NoInteract_l2_water, NoInteract_l_cell[i])
        else:
            WaterLine_2_Bkg = np.append(WaterLine_2_Bkg, NoInteract_l_cell[i])

    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(WaterLine_2_Bkg)):
        rdn_water_angle = rd.uniform(0,1)*np.pi*2
        WaterLine_2_Bkg[i] = rdn_water_angle

    #-------------------------------------------------------------------------------------------
    #Ahora los supervivientes del agua pasan al aire
    
    Lambda_air = 4124.77 #cm updated by nist xcom
    Pinteract_l2_air = np.array([]) #probability of interact in the line 2 air
    print('\n')
    print('Interaction with the air line 2, calculating backgrounds...')
    
    for i in range(0,len(NoInteract_l2_water)):
        P_air = ProbabilityInteract(l2_air, Lambda_air)
        Pinteract_l2_air = np.append(Pinteract_l2_air, P_air)

    NoInteract_l2_air = np.array([]) #fotones supervivientes de la célula que continuarán al aire siguiente
    for i in range(0,len(Pinteract_l2_air)):
        rdn_water = rd.uniform(0,1)
        if (Pinteract_l2_air[i] < rdn_water):
            NoInteract_l2_air = np.append(NoInteract_l2_air, NoInteract_l2_water[i])
        else:
            AirLine_2_Bkg = np.append(AirLine_2_Bkg, NoInteract_l2_water[i])

    #ahora sobre los fotones que interaccionaron son bkg y estimamos los ángulos theta con el que salen que es un random
    for i in range(0, len(AirLine_2_Bkg)):
        rdn_air_angle = rd.uniform(0,1)*np.pi*2
        AirLine_2_Bkg[i] = rdn_air_angle

    return [AirLine_1_Bkg, WaterLine_1_Bkg, CellPhotonsArray, WaterLine_2_Bkg, AirLine_2_Bkg]
'''

def InteractDetectorBackgrounds(a_det, b_det, c_det, l1_air, l2_air, AirLine_1_Bkg, WaterLine_1_Bkg, WaterLine_2_Bkg, AirLine_2_Bkg, Lambda_det):
    """With the angles of each background independently we can evaluate the permitted angles and the
    interaction into de detector
    """
    
    #-------------------------------------Air line 1 interactions -----------------------------
    #permitted angles  for Air line 1
    permitted_AirLine_1_Bkg = np.array([])
    theta1_forbiden_l1_Air = (np.pi/2.0) -  np.arctan((b_det/2.0 + l1_air/2.0)/c_det) #desde el centro de la linea de aire
    theta2_forbiden_l1_Air = (np.pi/2.0) +  np.arctan((b_det/2.0 - l1_air/2.0)/c_det) #desde el centro de la linea de aire
    
    for i in range(0,len(AirLine_1_Bkg)):
        if ((AirLine_1_Bkg[i] >= theta1_forbiden_l1_Air) & (AirLine_1_Bkg[i] <= theta2_forbiden_l1_Air)):
            permitted_AirLine_1_Bkg = np.append(permitted_AirLine_1_Bkg, AirLine_1_Bkg[i])

    r_det_l1_air = np.array([])
    for i in range(0,len(permitted_AirLine_1_Bkg)):
        #Calculation of the delta-r path into the detector aka r_det for the air photons 
        if (permitted_AirLine_1_Bkg[i] < np.pi/2.0):
            r_det_l1_air = np.append(r_det_l1_air, a_det/math.sin(permitted_AirLine_1_Bkg[i]))
        else:
            r_det_l1_air = np.append(r_det_l1_air, a_det/math.sin(np.pi - permitted_AirLine_1_Bkg[i]))
        
    #Air bkgs into detector    
    Pinteract_det_l1_air = np.array([]) #probability of interactions with the detector
    MeasuredPhoton_l1_air = np.array([]) #measured photons with the detector (the final survival photons)
    print('\n')
    print('Measuring the air 1 background photon events with the detector...')

    for i in range(0,len(permitted_AirLine_1_Bkg)):
        P_det = ProbabilityInteract(r_det_l1_air[i], Lambda_det)
        Pinteract_det_l1_air = np.append(Pinteract_det_l1_air, P_det)

    for i in range(0,len(Pinteract_det_l1_air)):
        rdn_det = rd.uniform(0,1)
        if (Pinteract_det_l1_air[i] > rdn_det):
            MeasuredPhoton_l1_air = np.append(MeasuredPhoton_l1_air, permitted_AirLine_1_Bkg[i])


    #-----------------------------------Water line 1 interactions--------------------------------
    #permitted angles  for Water line 1
    permitted_WaterLine_1_Bkg = np.array([])
    theta1_forbiden_l1_Water = (np.pi/2.0) - np.arctan((b_det/2.0)/c_det) #se aproxima que sale desde la celula central
    theta2_forbiden_l1_Water = (np.pi/2.0) + np.arctan((b_det/2.0)/c_det) #se aproxima que sale desde la celula central
    
    for i in range(0,len(WaterLine_1_Bkg)):
        if ((WaterLine_1_Bkg[i] >= theta1_forbiden_l1_Water) & (WaterLine_1_Bkg[i] <= theta2_forbiden_l1_Water)):
            permitted_WaterLine_1_Bkg = np.append(permitted_WaterLine_1_Bkg, WaterLine_1_Bkg[i])

    r_det_l1_water = np.array([])
    for i in range(0,len(permitted_WaterLine_1_Bkg)):
        #Calculation of the delta-r path into the detector aka r_det for the water photons 
        if (permitted_WaterLine_1_Bkg[i] < np.pi/2.0):
            r_det_l1_water = np.append(r_det_l1_water, a_det/math.sin(permitted_WaterLine_1_Bkg[i]))
        else:
            r_det_l1_water = np.append(r_det_l1_water, a_det/math.sin(np.pi - permitted_WaterLine_1_Bkg[i]))
            
    #Water bkgs into detector    
    Pinteract_det_l1_water = np.array([]) #probability of interactions with the detector
    MeasuredPhoton_l1_water = np.array([]) #measured photons with the detector (the final survival photons)
    print('\n')
    print('Measuring the water 1 background photon events with the detector...')

    for i in range(0,len(permitted_WaterLine_1_Bkg)):
        P_det = ProbabilityInteract(r_det_l1_water[i], Lambda_det)
        Pinteract_det_l1_water = np.append(Pinteract_det_l1_water, P_det)

    for i in range(0,len(Pinteract_det_l1_water)):
        rdn_det = rd.uniform(0,1)
        if (Pinteract_det_l1_water[i] > rdn_det):
            MeasuredPhoton_l1_water = np.append(MeasuredPhoton_l1_water, permitted_WaterLine_1_Bkg[i])
    
    
    #-----------------------------------Water line 2 interactions--------------------------------
    #permitted angles  for Water line 2
    permitted_WaterLine_2_Bkg = np.array([])
    theta1_forbiden_l2_Water = (np.pi/2.0) - np.arctan((b_det/2.0)/c_det)
    theta2_forbiden_l2_Water = (np.pi/2.0) + np.arctan((b_det/2.0)/c_det)
    
    for i in range(0,len(WaterLine_2_Bkg)):
        if ((WaterLine_2_Bkg[i] >= theta1_forbiden_l2_Water) & (WaterLine_2_Bkg[i] <= theta2_forbiden_l2_Water)):
            permitted_WaterLine_2_Bkg = np.append(permitted_WaterLine_2_Bkg, WaterLine_2_Bkg[i])

    r_det_l2_water = np.array([])
    for i in range(0,len(permitted_WaterLine_2_Bkg)):
        #Calculation of the delta-r path into the detector aka r_det for the water photons
        if (permitted_WaterLine_2_Bkg[i] < np.pi/2.0):
            r_det_l2_water = np.append(r_det_l2_water, a_det/math.sin(permitted_WaterLine_2_Bkg[i]))
        else:
            r_det_l2_water = np.append(r_det_l2_water, a_det/math.sin(np.pi - permitted_WaterLine_2_Bkg[i]))

    #Water bkgs into detector    
    Pinteract_det_l2_water = np.array([]) #probability of interactions with the detector
    MeasuredPhoton_l2_water = np.array([]) #measured photons with the detector (the final survival photons)
    print('\n')
    print('Measuring the water 2 background photon events with the detector...')

    for i in range(0,len(permitted_WaterLine_2_Bkg)):
        P_det = ProbabilityInteract(r_det_l2_water[i], Lambda_det)
        Pinteract_det_l2_water = np.append(Pinteract_det_l2_water, P_det)

    for i in range(0,len(Pinteract_det_l2_water)):
        rdn_det = rd.uniform(0,1)
        if (Pinteract_det_l2_water[i] > rdn_det):
            MeasuredPhoton_l2_water = np.append(MeasuredPhoton_l2_water, permitted_WaterLine_2_Bkg[i])

    
    #-------------------------------------Air line 2 interactions -----------------------------
    #permitted angles  for Air line 2
    permitted_AirLine_2_Bkg = np.array([])
    theta1_forbiden_l2_Air = (np.pi/2.0) -  np.arctan((b_det/2.0 - l2_air/2.0)/c_det) #salen desde el medio de la linea de aire
    theta2_forbiden_l2_Air = (np.pi/2.0) +  np.arctan((b_det/2.0 + l2_air/2.0)/c_det) #salen desde el medio de la linea de aire
    
    for i in range(0,len(AirLine_2_Bkg)):
        if ((AirLine_2_Bkg[i] >= theta1_forbiden_l2_Air) & (AirLine_2_Bkg[i] <= theta2_forbiden_l2_Air)):
            permitted_AirLine_2_Bkg = np.append(permitted_AirLine_2_Bkg, AirLine_2_Bkg[i])

    r_det_l2_air = np.array([])
    for i in range(0,len(permitted_AirLine_2_Bkg)):
        #Calculation of the delta-r path into the detector aka r_det for the air photons 
        if (permitted_AirLine_2_Bkg[i] < np.pi/2.0):
            r_det_l2_air = np.append(r_det_l2_air, a_det/math.sin(permitted_AirLine_2_Bkg[i]))
        else:
            r_det_l2_air = np.append(r_det_l2_air, a_det/math.sin(np.pi - permitted_AirLine_2_Bkg[i]))
        
    #Air bkgs into detector    
    Pinteract_det_l2_air = np.array([]) #probability of interactions with the detector
    MeasuredPhoton_l2_air = np.array([]) #measured photons with the detector (the final survival photons)
    print('\n')
    print('Measuring the air 2 background photon events with the detector...')

    for i in range(0,len(permitted_AirLine_2_Bkg)):
        P_det = ProbabilityInteract(r_det_l2_air[i], Lambda_det)
        Pinteract_det_l2_air = np.append(Pinteract_det_l2_air, P_det)

    for i in range(0,len(Pinteract_det_l2_air)):
        rdn_det = rd.uniform(0,1)
        if (Pinteract_det_l2_air[i] > rdn_det):
            MeasuredPhoton_l2_air = np.append(MeasuredPhoton_l2_air, permitted_AirLine_2_Bkg[i])
    
    
    return [MeasuredPhoton_l1_air, MeasuredPhoton_l1_water, MeasuredPhoton_l2_water, MeasuredPhoton_l2_air]
  
     
'''
#----------------------------------------------- Efficiencies -------------------------------------
'''

def Efficiencies(UpperPhotons_XR, Upper_NoWaterInteract_XR, Upper_NoAirInteract_XR, permitted_XR, MeasuredPhoton_XR):
    
    #Efficiency in the water block, survival photons vs the upper photons
    try:
        Eff_water = len(Upper_NoWaterInteract_XR)/len(UpperPhotons_XR)
        print('\n')
        print('The efficiency in the water is ', Eff_water*100.0, '%')
    except:
        Eff_water = 0.0
        print('\n')
        print('The efficiency in the water is ', Eff_water*100.0, '%')
        
    #Efficiency in the Air block, survival photons vs survival photons after water
    try:
        Eff_air = len(Upper_NoAirInteract_XR)/len(Upper_NoWaterInteract_XR)
        print('\n')
        print('The efficiency in the air is ', Eff_air*100.0, '%')
    except:
        Eff_air = 0.0
        print('\n')
        print('The efficiency in the air is ', Eff_air*100.0, '%')
           
    #Efficiency of the photons contained in the detection area between theta forbiden angles (geometric aceptance)
    try:
        A_geo = len(permitted_XR)/len(Upper_NoAirInteract_XR)
        print('\n')
        print('The geometric aceptance is ', A_geo*100.0, '%')
    except:
        A_geo = 0.0
        print('\n')
        print('The geometric aceptance is ', A_geo*100.0, '%')
    
    #Efficiency of the photons measured in the detector over permitted (detector efficiency)
    try:
        Eff_det = len(MeasuredPhoton_XR)/len(permitted_XR)
        print('\n')
        print('The efficiency in the detector is ', Eff_det*100.0, '%')
    except:
        Eff_det = 0.0
        print('\n')
        print('The efficiency in the detector is ', Eff_det*100.0, '%')
    
    #Efficiency of the photons measured in the detector over all generated before pass trough the Air block (general efficiency)
    try:
        Eff_gen = len(MeasuredPhoton_XR)/len(UpperPhotons_XR)
        print('\n')
        print('The general efficiency is ', Eff_gen*100.0, '%')
    except:
        Eff_gen = 0.0
        print('\n')
        print('The general efficiency is ', Eff_gen*100.0, '%')

    return [Eff_water, Eff_air, A_geo, Eff_det, Eff_gen]


'''
#------------------------------------------------- Plotter ----------------------------------------
'''
def PlotterEfficiencyDistributions(n, InitEnergy, UpperPhotons_XR, Upper_NoAirInteract_XR, permitted_XR, MeasuredPhoton_XR): 
    # Distribution of the photon diferents efficiencies on each part
    fig=plt.figure(figsize=(20,10))
    fig.subplots_adjust(top=0.90, bottom=0.10, hspace=0.2, wspace=0.2)

    ax=fig.add_subplot(1,1,1)
    plt.title(r'%i X-ray photon distribution in the detector E = %f MeV' %(len(n), InitEnergy), fontsize=15)
    plt.hist([UpperPhotons_XR, Upper_NoAirInteract_XR, permitted_XR, MeasuredPhoton_XR], bins = 25, histtype='bar', linewidth=3, 
             color=['blue', 'red', 'green', 'purple'], label= [r'upper # = %i' %len(UpperPhotons_XR), r'no air interact # = %i' %len(Upper_NoAirInteract_XR),
                   r'permitted # = %i' %len(permitted_XR), r'measured # = %i' %len(MeasuredPhoton_XR)])
    ax.set_xlabel(r'photons along the detector size', fontsize=15)
    ax.set_ylabel('photon distributions', fontsize=15)
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return fig


#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    print('\n')
    print('Ejecutando como programa principal \n')
  