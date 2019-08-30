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
    
def main(inputs, geopath, meshpath):
    
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
    saveDt = inputs.savedt
    #solutions = []
    
    # Fluid Properties
    rho = Constant(inputs.rho0)
    mu = Constant(inputs.mu0)
    
    # Start timer
    start = timeit.default_timer()
    
    # Time Loop
    w0 = []
    while t <= inputs.tEnd:
           
        # Set Dynamic Timestep
        dt = dynamicTimestep(t,inputs.dtMax,inputs.dtMin,inputs.tChange)
        
        # Initialize results Vector
        results = []
        
        begin('Flow - Time:{:.3f}s'.format(t))
        
        # Solve Equations
        w0 = transientFlow(w0,dt,rho,mu,inputs,meshObj,boundaries,Subdomains)
        end()
        
        (u1, p1) = w0.leaf_node().split()
        
        
        # Save Paraview Files
        if t==0 or t >= saveDt:
            begin('---------------- Saving ----------------')
            # Append and save Results
            results.append(u1)
            results.append(p1)
            saveResults(results,paraFiles,inputs.ParaViewFilenames,inputs.ParaViewTitles)
            saveDt = t + inputs.savedt
            end()
            # Store Initial Solution in Time t=t
        #        solutions.append({'t':t, 'variables':results})
                
        # Update current time
        t += dt    
            
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
    
    #%%################     Save Simulation Properties and Results  : TO DO
    #file = open(fileDir+'/Results_'+tag+'.txt','w') 
    #file.write('----'+tag +'\n') # Title
    #
    #file.write('-- Properties of '+'\n') # Properties
    #file.write("Total running time: %dh:%dmin:%ds" % (hours, mins, secs)+'\n') 
    #file.write('rhoH20:{:.0f}'.format(rho0)+' | rhoCement:{:.0f}'.format(rho1)+'\n') 
    #file.write('muH20:{:.0f}'.format(mu0)+' | muCement:{:.0f}'.format(mu1)+'\n')
    #file.write('Newtonian Fluids'+'\n')
    #file.write('Density not a function of time'+'\n') 
    #
    #file.write('-- Results of '+ '\n') # Results
    #file.write('Mass Out:') 
    #
    #file.close() 
     
if __name__ == '__main__':
    inputs = Inputs()
    mainPaths = Paths()
    
    # Check directories for saving
#    meshFile = savingCheckings(inputs.meshFile,inputs.caseId,mainPaths)

    # Prepare for saving
    main(inputs,mainPaths.geopath,mainPaths.meshpath)