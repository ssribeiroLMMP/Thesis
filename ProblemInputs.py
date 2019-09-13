#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 08:21:14 2019

@author: Sergio Ribeiro
Description: Problem Definitions: Input Dictionary
"""
from dolfin import *
import math

# Sigmoid Function S(x) = (dtMax-dtMin)/(1+exp(Inclination*(-t+t0)))+dtMin
def dynamicTimestep(t,dtMax,dtMin,t0): 
    Inclination = 3
    return (dtMax-dtMin)/(1+math.exp(Inclination*(-t+t0)))+dtMin

def autoTimestep(no_iterations,dt,inputs,MinIterations=4, MaxIterations=5,increment=1.5):
    # Check if 
    if no_iterations < MinIterations:
        dt = min(increment*dt,inputs.dtMax)
    elif no_iterations > MaxIterations:
        dt = max((1/increment)*dt,inputs.dtMin)
    
    return dt

class Inputs():
    def __init__(self):
        #%%############ Case Definition    ##############################
        self.caseId = 'WellSimulator_FiltrationTest' ## If name already exists in folder ./PostProcessing/Cases, 
                         ## old data will be overwritten.
        
        # Output Variables
        self.ParaViewFilenames = []; self.ParaViewTitles = []
        self.ParaViewFilenames.append("velocity"); self.ParaViewTitles.append('Velocity (m/s)')
        self.ParaViewFilenames.append("pressure"); self.ParaViewTitles.append('Pressure (Pa)')
        self.ParaViewFilenames.append("concentration"); self.ParaViewTitles.append('Mass Fraction (Fluid Tags)')
        
        #%%############ Gravitationa Field ##############################
        # Gravity Acceleration (m/s²)
        self.g =  9.81
        
        #%%############ Fluids' Properties ##############################
        # Tags
        self.Fluid0 = 0 # Cement 
        self.Fluid1 = 1 - self.Fluid0 # Water
                
        # Density (kg/m³)
        self.rho_values = [1000, 1000]
        
        # Initial Interface position
        self.InterfaceX0 = 0.05
        # self.InterfaceY0 = 0.01
        
        # Diffusity of Between species (m²/s)
        self.D = 10
        
        # Rheology
        # Newtonian Viscosity (Pa.s)
        self.mu_values = [0.01 , 0.01]
        
        #%%############ Time Parameters #################################
        # Start Time
        self.t0 = 0 # s
        
        # Simulation Time
        self.tEnd = 10 # s
        
        # Variable time step
        self.dtMin = 0.0000001    # s
        self.dtMax = 1  # s
        self.tChange = 0   # point in time of sigmoid inflection occurs (s)
        self.dt = 0.01
#        self.dt = dynamicTimestep(self.t0,self.dtMax,self.dtMin,self.tChange)    # s
        
        #%%############ Logging Options   ###############################
        # Result Saving time step
        self.savedt = 0.01 # s
        
        #%%############ Problem Geometry   ##############################
        ## Mesh File
        self.meshFile = 'WellSimulator'
        ## Mesh Elements
        # Velocity
        self.velocityElementfamily = 'Lagrange'
        self.velocityElementOrder = 2
        # Pressure
        self.pressureElementfamily = 'Lagrange'
        self.pressureElementOrder = 1
        # Advected Scalar Field
        self.scalarFieldElementfamily = 'Lagrange'
        self.scalarFieldElementOrder = 1
        
        #%%############ Boundary Conditions
        
        ## No slip Boundaries
        self.noSlipBCs = []
        self.noSlipBCs.append('InnerPipe')
        self.noSlipBCs.append('OuterWall')
        self.noSlipBCs.append('BottomWall')
#        self.noSlipBCs.append('TopWall')
#        self.noSlipBCs.append('BottomWall')
        #noSlipBoundaries.append('InnerWalls')
        
        ## Pressure Inputs
        self.pressureBCs = {}
        self.pressureBCs.update({'Inlet' : 1000}) # Pa
#        self.pressureBCs.update({'Outlet' : 0}) # Pa
        
        ## Advected Scalar Field Inputs
        self.CInitialMixture = 0.5      # Mass fraction of Fluid0. Fluid1 = 1-Flui0
        self.scalarFieldBCs = {}
        self.scalarFieldBCs.update({'Inlet' : Constant(self.CInitialMixture)}) # CMix
        self.scalarFieldBCs.update({'Outlet': Constant(self.Fluid1)}) # C1
        
        ## Velocity Inputs
        self.velocityBCs = {}
        self.velocityBCs.update({'Oulet' : Constant(1e-5)}) # m/s
        
        #%%############ Solver parameters ###############################
        # Absolute Tolerance    
        self.absTol = 1e-12
        
        # Relative Tolerance
        self.relTol = 1e-10
        
        # Maximum Iterations
        self.maxIter = 15
        
        # Linear Solver
        self.linearSolver = 'mumps'
            
        # Relaxation Factors
        self.alpha = 1
        self.alphaC = 1
            
        #%% Possible Solvers
        # Solver method  |  Description    
        # -----------------------------------------------------------
        # bicgstab       |  Biconjugate gradient stabilized method                      
        # cg             |  Conjugate gradient method                                   
        # default        |  default linear solver                                       
        # gmres          |  Generalized minimal residual method                         
        # minres         |  Minimal residual method                                     
        # mumps          |  MUMPS (MUltifrontal Massively Parallel Sparse direct Solver)
        # petsc          |  PETSc built in LU solver                                    
        # richardson     |  Richardson method                                           
        # superlu        |  SuperLU                                                     
        # tfqmr          |  Transpose-free quasi-minimal residual method                
        # umfpack        |  UMFPACK (Unsymmetric MultiFrontal sparse LU factorization) 
        #
        # LU method  |  Description                                                 
        # --------------------------------------------------------------------------
        # default    |  default LU solver                                           
        # mumps      |  MUMPS (MUltifrontal Massively Parallel Sparse direct Solver)
        # petsc      |  PETSc built in LU solver                                    
        # superlu    |  SuperLU                                                     
        # umfpack    |  UMFPACK (Unsymmetric MultiFrontal sparse LU factorization) 
        #
        # Krylov method  |  Description                                 
        # --------------------------------------------------------------
        # bicgstab       |  Biconjugate gradient stabilized method      
        # cg             |  Conjugate gradient method                   
        # default        |  default Krylov method                       
        # gmres          |  Generalized minimal residual method         
        # minres         |  Minimal residual method                     
        # richardson     |  Richardson method                           
        # tfqmr          |  Transpose-free quasi-minimal residual method
        # 
        # 