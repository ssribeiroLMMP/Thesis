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

def autoTimestep(no_iterations,dt,inputs,limitIterations=4,increment=1.1):
    # Check if 
    if no_iterations < limitIterations:
        dt = min(increment*dt,inputs.dtMax)
    else:
        dt = max((1/increment)*dt,inputs.dtMin)
    
    return dt

class Inputs():
    def __init__(self):
        #%%############ Case Definition    ##############################
        self.caseId = 'vargesPR_HeleShawCell_sigma0_01_miStar_5e0_CSF_CurvedInterface_NavierSlip'
#        self.caseId = 'vargesPR_HeleShawCell_etaStar5_3_RectTriangMesh' ## If name already exists in folder ./PostProcessing/Cases, 
                         ## old data will be overwritten.
        
        # Output Variables
        self.ParaViewFilenames = []; self.ParaViewTitles = []
        self.ParaViewFilenames.append("velocity"); self.ParaViewTitles.append('Velocity (m/s)')
        self.ParaViewFilenames.append("pressure"); self.ParaViewTitles.append('Pressure (Pa)')
        self.ParaViewFilenames.append("concentration"); self.ParaViewTitles.append('Mass Fraction (Fluid Tags)')
        
        #%%############ Gravitationa Field ##############################
        # Gravity Acceleration (m/s²)
        self.g = 9.81
        
        #%%############ Fluids' Properties ##############################
        # Tags
        self.FluidTags = [0,1]
                
        # Density (kg/m³)
        self.rho_values = [1191, 867.6]
        
        # Initial Interface position
        self.InterfaceR = 0.05
        self.InterfaceX0 = 0.02 + self.InterfaceR
        self.InterfaceY0 = 0.071
        
        # Diffusity of Between species (m²/s)
        self.D = 1e-6
        
        # Interfacial Tension (N/m)
        self.sigma = 0.01
        # slip Coefficient (m)
        self.dSlip = 1e-6 # (1 um)
                
        
        # Rheology
        # Newtonian Viscosity (Pa.s)
        self.mu_values = [0.0251,0.134]
        self.mu_star = self.mu_values[0]/self.mu_values[1]
        
        #%%############ Time Parameters #################################
        # Start Time
        self.t0 = 0 # s
        # Simulation Time
        self.tEnd = 10 # s
        
        # Variable time step
        self.dtMin = 1e-3   # s
        self.dtMax = 1e-1  # s
        self.tChange = 0   # point in time of sigmoid inflection occurs (s)
        self.dt = 5e-2
#        self.dt = dynamicTimestep(self.t0,self.dtMax,self.dtMin,self.tChange)    # s
        
        #%%############ Logging Options   ###############################
        # Result Saving time step
        self.savedt = 1e-2 # s
        
        #%%############ Problem Geometry   ##############################
        ## Mesh File
        self.meshFile = 'OriginalHeleShaw'
        self.CellThickness = 0.7e-3
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
               
        ## Velocity Boundaries
        # Navier-Slip
        self.navierslip = {}
        # self.navierslip.update({'BottomWall' : self.dSlip})
        # self.navierslip.update({'TopWall'    : self.dSlip})
        
        
        # No-Slip
        self.noSlipBCs = []
        # self.noSlipBCs.append('TopWall')
        # self.noSlipBCs.append('BottomWall')
#        self.noSlipBCs.append('RightWall')
#        self.noSlipBCs.append('LeftWall')
        #noSlipBoundaries.append('InnerWalls')
        # Inputs
        self.velocityBCs = {}                #  Vy = 0 m/s
        self.velocityBCs.update({'TopWall'    : {2 : 0.0}}) 
        self.velocityBCs.update({'BottomWall' : {2 : 0.0}})
        # self.VxInlet = 1e-2
        # self.velocityBCs.update({'Inlet' : Constant((self.VxInlet,0.0))}) # m/s
        
        
        ## Pressure Inputs
        self.pressureBCs = {}
        self.PInlet = 5
        self.POutlet = 0
        self.pressureBCs.update({'Inlet' : self.PInlet}) # Pa
        self.pressureBCs.update({'Outlet' : self.POutlet}) # Pa
        
        ## Advected Scalar Field Inputs
        self.CInlet = self.FluidTags[0]
        self.scalarFieldBCs = {}
        self.scalarFieldBCs.update({'Inlet' : Constant(self.CInlet)}) # % Only fluid 0 enters
        #self.scalarFieldBCs.update({'TopWall': self.TTopWall}) # T
        
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
        self.alpha = 0.85
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
