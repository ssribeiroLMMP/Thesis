#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 08:55:13 2019

@author: sergio

Properties:
    Final:
    Momentum and Mass Conservation
    Body Forces Term: Gravity
    Concentration advection and diffusion
    
    Partial:
    Newtonian Fluids
    Density Constant with time

"""
#%%#################    Libraries and Pre-Defined Functions  
# Third Party Packages
import numpy as np
import shelve
import timeit
import sys
import os
sys.path.append(os.path.abspath('.'))
from ProblemInputs import *
from PreProcessing.meshConversion import *
from Solver.Equations import *
from PostProcessing.Saving import *

def main(inputs):
    
    # Prepare for saving
    meshpath, paraFullPath, imFullPath,parapath, imagespath,geopath,  \
    postpath = assemblePaths(inputs)
    paraFiles = savingPreparation(paraFullPath,inputs.ParaViewFilenames)
    
    #%%#################    Mesh Applicaticvvon  
    # Type 1 - from XML File            -- OK
    # Type 2 - Converts from .msh       -- to do
    # Type 3 - Manual FEniCs Creation   -- to do
    meshObj, boundaries, markers, Subdomains = applyMesh(1,geopath,meshpath,inputs.meshFile)
    
    # Get Mesh Dimension: 1, 2 or 3
    Dim = meshObj.geometric_dimension()
    
    #%%#################    Time Loop - If Transient
    t = inputs.t0
    saveDt = min([inputs.savedt, min(inputs.plotTimeList)])
    #solutions = []
    
    # Mixture Properties
    D = Constant(inputs.D)
    
    # Start timer
    start = timeit.default_timer()
    
    # Space Functions: Flow and Scalar Field
    W = flowSpaceCreation(inputs,meshObj)
    C = fieldSpaceCreation(inputs,meshObj)
    w0 = Function(W)
    rho_cem_t = project(inputs.shrinkage_rhoMax,C)
    
    # Initial Conditions
    c0 = initialMixture(C,inputs)
        
    # Timestep
    dt = inputs.dt
    
    # Initialize loop variable values
    lastStep = False
    lastTry = False
    outputMassFlowrate = 0
    pInlet = inputs.pInlet
    listId = 0
    cbarU = 0
    cbarP = 0
    TOC = 0

    # Create CSV File
    createCSVOutput(inputs.outputFlowrate,inputs.fieldnamesFlow)
    createCSVOutput(inputs.outputPressure,inputs.fieldnamesPre)
    initializeVelocityProfile(inputs)
    initializePressureProfile(inputs)
    initializePressurePerDepth(inputs)

    TOC = inputs.Zmin
    rhoMix = inputs.CInitialMixture*inputs.rho_values[0] + \
            (1-inputs.CInitialMixture)*inputs.rho_values[1]
    
    while t <= inputs.tEnd:
        rho_cem_t0 = rho_cem_t
        # Initialize results Vector
        results = []
        
        # Assign Fluids Properties
        (u0, p0) = w0.leaf_node().split()
        
        #assignFluidProperties(inputs,C,c,u,t)
        rho,mu,rho_cem_t = assignFluidProperties(inputs,C,c0,u0,t)

        # Set Initial Density Field
        if t<=inputs.dt:
            rho0=rho
    
        # Solve Equations
        # try:
        begin('Flow - Time:{:.3f}s and dt:{:.5f}s'.format(t,dt))
        if t>0:
            pInlet, TOC, rhoMix = calculateNewInletPressure(TOC,outputMassFlowrate,rho,t,dt,boundaries,Subdomains,inputs)
            
        (w,no_iterations,converged) = transientFlow(t,W,w0,dt,rho,rho0,mu,inputs,meshObj,boundaries,Subdomains,pInlet)
        end()
        # except: 
        #     no_iterations = inputs.maxIter
        #     converged = False
        #     begin('Reducing timestep')
        #     end()
        #     end()
        
        if converged:
            (u1, p1) = w.leaf_node().split()
            outputMassFlowrate = calculateOutletFlowrate(u1,inputs,boundaries,Subdomains)
            
            begin('Concentration')
            c1 = transienFieldTransport(C,c0,dt,u1,D,rho_cem_t,rho_cem_t0,mu,inputs,meshObj,boundaries,Subdomains)
            # c1 = c0
            end()
            
            rho0 = rho
            
            # Save Paraview Files
            if t==0 or t >= saveDt:
                begin('---------------- Saving ----------------')
                # Append and save Results
                results.append(u1)
                results.append(p1)
                results.append(c1)
                results.append(rho)
                results.append(mu)
                saveResults(results,paraFiles,inputs.ParaViewFilenames,inputs.ParaViewTitles)
                # Pressure per Depth
                if not(t==0):
                    savePressurePerDepthDuringRun(inputs,p1,t)
                
                # Save Flowrate
                results.append(outputMassFlowrate)
                # Save TOC
                results.append(TOC)
                # Save Fluid Properties
                results.append(rhoMix)
                # results.append(rho)
                # results.append(mu)
                cbarU, cbarP, listId = logResults(t,results,listId,inputs,meshObj,cbarU=cbarU,cbarP=cbarP)
                if listId < len(inputs.plotTimeList):
                    saveDt = min([t + inputs.savedt, inputs.plotTimeList[listId]])
                end()
                # Store Initial Solution in Time t=t
                # solutions.append({'t':t, 'variables':results})
                    
        	   # Update current time #ERROR ON VERSION 1.0.4
               # Log into CSVFiles
                       
            
            w0.assign(w)
            c0.assign(c1)
            t += dt
        
        if t > inputs.tEnd and not(lastStep):
            t = inputs.tEnd
            lastStep = True
            # Get last timestep velocity profile: Only turn on in MeshTest cases
            if inputs.caseId.find('MeshTest'):
                saveVelocityProfileDuringRun(inputs,u1)

        ### Set next timestep
        #Dynamic Timestep
        # dt = dynamicTimestep(t,inputs.dtMax,inputs.dtMin,inputs.tChange)
        # Auto-adjustable timestep
        dt = autoTimestep(no_iterations,dt,inputs)
        if dt <= inputs.dtMin and lastTry:
            begin('Did not converge: Minimum timestep achieved')
            end()
            break
        elif dt <= inputs.dtMin:
            lastTry = True
          
    #####################  Post Processing
    print('Finished')
    # End Time
    stop = timeit.default_timer()
    total_time = stop - start
    
    # Output running time in a nice format.
    mins, secs = divmod(total_time, 60)
    hours, mins = divmod(mins, 60)
    
    begin("Total running time: %dh:%dmin:%ds \n" % (hours, mins, secs))
    end()
    
     
if __name__ == '__main__':
    inputs = Inputs()
    mainPaths = Paths()
    
    # Check directories for saving
    savingCheckings(inputs,mainPaths)

    # Prepare for saving
    main(inputs)
