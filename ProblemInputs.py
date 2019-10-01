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

def autoTimestep(no_iterations,dt,inputs,limitIterations=4,increment=2):
    # Check if 
    if no_iterations < limitIterations:
        dt = min(increment*dt,inputs.dtMax)
    else:
        dt = max((1/increment)*dt,inputs.dtMin)
    
    return dt

class Inputs():
    def __init__(self):
        #%%############ Case Definition    ##############################
        self.caseId = 'TransPoiseulle_nonNewtonian_ShortTest' ## If name already exists in folder ./PostProcessing/Cases, 
                         ## old data will be overwritten.
        
        # Output Variables
        self.ParaViewFilenames = []; self.ParaViewTitles = []
        self.ParaViewFilenames.append("velocity"); self.ParaViewTitles.append('Velocity (m/s)')
        self.ParaViewFilenames.append("pressure"); self.ParaViewTitles.append('Pressure (Pa)')
        self.ParaViewFilenames.append("concentration"); self.ParaViewTitles.append('Mass Fraction (Fluid Tags)')
        
        #%%############ Gravitationa Field ##############################
        # Gravity Acceleration (m/s²) on axis X
        self.g = 0
        
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
        self.D = 1e-1
        
        # Rheology
        # Newtonian Viscosity (Pa.s)
        self.mu_values = [0.01 , 0.01]  #
        # Modified SMD Model Variables
        self.tau0 = 19.19           # Dinamic Yield Stress               
        self.etaInf = 0.295         # Equilibrium Viscosity(Newtonian Plato: High shear rates)
        self.eta0 = 1e3             # Newtonian Plato: Low shear rates
        self.K = 1.43               # Consistency Index
        self.n = 0.572              # Power-law Index
        self.ts = 60000             # Caracteristic curing time


        #%%############ Time Parameters #################################
        # Start Time
        self.t0 = 0 # s
        
        # Simulation Time
        self.tEnd = 200 # s
        
        # Variable time step
        self.dtMin = 0.005    # s
        self.dtMax = 10  # s
        self.tChange = 0   # point in time of sigmoid inflection occurs (s)
        self.dt = 0.1
#        self.dt = dynamicTimestep(self.t0,self.dtMax,self.gging Options   ###############################
        # Result Saving time step
        self.savedt = 5*self.dt # s
        
        #%%############ Problem Geometry   ##############################
        ## Mesh File
        self.meshFile = 'MacroParallelPlates'#'WellSimulator'#
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
        self.noSlipBCs.append('TopWall') # Both
        # self.noSlipBCs.append('InnerPipe') # WelSimulator only
        # self.noSlipBCs.append('OuterWall') # WelSimulator only
        self.noSlipBCs.append('BottomWall')  # Both
        #noSlipBoundaries.append('InnerWalls')
        
        ## Pressure Inputs
        self.pressureBCs = {}
        self.pInlet = 0.3164557 #self.rho_values[0]*2*self.g
        self.pressureBCs.update({'Inlet' : self.pInlet}) # Pa
        # self.pOutlet = 0.2*(self.pInlet + self.rho_values[0]*self.g*1)
        self.pressureBCs.update({'Outlet' : 0}) # Pa
        
        ## Advected Scalar Field Inputs
        self.CInitialMixture = 0.5      # Mass fraction of Fluid0. Fluid1 = 1-Flui0
        self.scalarFieldBCs = {}
        self.scalarFieldBCs.update({'Inlet' : Constant(self.CInitialMixture)}) # CMix
        self.scalarFieldBCs.update({'Outlet' : Constant(self.CInitialMixture)}) # CMix: REGULAR FLOW
        # self.scalarFieldBCs.update({'Outlet': Constant(self.Fluid1)}) # C1: FILTRATION
        
        ## Velocity Inputs
        t=0
        self.velocityBCs = {}
        # self.VxInlet = '0.0025+0*t'
        # self.VrOutlet = '0.0000043+0*t'#'2*exp(1-(t/200))/300'#
        # self.velocityBCs.update({'Inlet' : Expression((self.VxInlet,'0.0'),t=t,degree=1)}) # m/s
        # self.velocityBCs.update({'Outlet' : Expression((self.VrOutlet,'0.0'),t=t,degree=1)}) # m/s
        
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
        self.alpha = 0.9
        self.alphaC = 0.9
            
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