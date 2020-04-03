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
    elif no_iterations > limitIterations + 2:
        dt = max((1/increment)*dt,inputs.dtMin)
    else:
        dt = dt
    
    return dt

def dynamicSaveDt(dt):
    return 5*dt

class Inputs():
    def __init__(self):
        #%%############ Case Definition    ##############################
        self.caseId = 'TransWellSimulator_BaseCase_14000_7' ## If name already exists in folder ./PostProcessing/Cases, 
                         ## old data will be overwritten.
        
        # Output Variables
        self.ParaViewFilenames = []; self.ParaViewTitles = []
        self.ParaViewFilenames.append("velocity"); self.ParaViewTitles.append('Velocity (m/s)')
        self.ParaViewFilenames.append("pressure"); self.ParaViewTitles.append('Pressure (Pa)')
        self.ParaViewFilenames.append("concentration"); self.ParaViewTitles.append('Mass Fraction (Fluid Tags)')
        self.outputFlowrate = './PostProcessing/Cases/'+self.caseId+'/flowrateInput.csv'
        self.outputPressure = './PostProcessing/Cases/'+self.caseId+'/pressureOutput.csv'
        self.outputTOC = './PostProcessing/Cases/'+self.caseId+'/TOC.csv'
        
    
        #%%############ Time Parameters #################################
        # Start Time
        self.t0 = 0 # s
        
        # Simulation Time
        self.tEnd = 14000 # s

        # Plot Time List
        self.plotTimeList = [1, 200, 2500, 4500, 7000, 8250, 9000, 10500, 12000, self.tEnd]#
        self.plotDepthList = [8, 7.6, 7.2, 6.8, 6.4, 6]
        # self.plotTimeList = [1, self.tEnd]
        self.fieldnamesFlow = ['Time(s)','outletFlowRate(Kg/s)']
        self.dZPlot = 0.01
        self.fieldnamesPre = ['Time(s)','rhoInlet(Kg/m3)']

        # Variable time step
        self.dtMin = 0.005    # s
        self.dtMax = 20  # s
        self.tChange = 0   # point in time of sigmoid inflection occurs (s)
        self.dt = 0.01
#        self.dt = dynamicTimestep(self.t0,self.dtMax,self.gging Options   ###############################
        # Result Saving time step
        self.savedt = 100 # s

        #%%############ Gravitationa Field ##############################
        # Gravity Acceleration (m/s²) on axis X
        self.g = 9.81
        
        #%%############ Fluids' Properties ##############################
        # Tags
        self.Fluid0 = 0 # Cement 
        self.Fluid1 = 1 - self.Fluid0 # Water
        self.CInitialMixture = 0.5      # Mass fraction of Fluid0. Fluid1 = 1-Flui0
                       
        # Initial Interface position: Two-Phase Flow
        # self.InterfaceX0 = 0.05
        # self.InterfaceY0 = 0.01
        
        # Diffusity of Between species (m²/s)
        self.D = 1e-3
        
        # Rheology
        # Newtonian Viscosity
        self.mu_water = 0.01    # Pa.s
        self.mu_cem = 0.01       # Pa.s
        self.mu_values = [self.mu_cem , self.mu_water]  
        
        ## Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t)) 
        # Input Variables
        self.tau0 = 19.019         # Dinamic Yield Stress               
        self.etaInf = 0.295        # Equilibrium Viscosity(Newtonian Plato: Lowgh shear rates)
        self.eta0 = 1e4            # Viscosity Value for Low shear rates
        self.K = 1.43              # Consistency Index
        self.n = 0.572             # Power-law Index
        self.ts = 6000             # Caracteristic viscosity buildup time
        self.eps = 1e-10
        self.tauY_t = Expression('tau0*exp(t/ts)',degree=1,\
                                  tau0 = self.tau0, \
                                  ts = self.ts, t = 0)
        # Modified SMD Model Expression
        self.rheologicalModel = Expression('(1 - exp(-eta0*(gammaDot)/tauY_t))* \
                                            (tauY_t/(gammaDot+eps) + \
                                            K*(pow(gammaDot+eps,nPow)/(gammaDot+eps))) \
                                            + etaInf', degree=2, \
                                            etaInf=self.etaInf, eta0=self.eta0, \
                                            K = self.K, nPow = self.n, ts = self.ts, \
                                            eps = self.eps, gammaDot = Constant(1/self.eps),\
                                            tauY_t = Constant(self.eps))

        # Density (kg/m³)
        # Experimental Values
        self.rho_water = 1000           # kg/m³
        self.rho_bulk0 = 1737.48        # kg/m³ = 14.5ppg
        self.rho_bulkInf = 1200      # kg/m³    
        self.rho_cem0 = (self.rho_bulk0 - (1-self.CInitialMixture)*self.rho_water)/(self.CInitialMixture)
        
        # Initial Density Values per Component
        self.rho_values = [self.rho_cem0 , self.rho_water] # kg/m³

        # Shrinkage Equation: rhoMax - ((rhoMax-rhoMin)/(1+math.exp(Inclination*(-t+t0)))+rhoMin) +rhoMin
        # if Inclination is zero, shrinkage is neglected
        self.shrinkage_inclination = 0.0004  # kg/m³ / s
        self.shrinkage_rhoMin = self.rho_bulkInf    # kg/m³
        self.shrinkage_t0 = 1.2*self.ts            # s  

        # Shrinkage Model Expression
        self.shrinkageModel = Expression('(rhoMax-((rhoMax-rhoMin)/  \
                                         (1 + exp(Inclination*(-t+t0))) \
                                         +rhoMin) + rhoMin)*cFrac + rho_water*(1-cFrac)', degree=2, \
                                         cFrac = Constant(self.eps), \
                                         rhoMax = self.rho_values[self.Fluid0], \
                                         rhoMin = self.shrinkage_rhoMin, \
                                         Inclination = self.shrinkage_inclination, \
                                         t=self.t0, t0 = self.shrinkage_t0, \
                                         rho_water = self.rho_values[self.Fluid1])
                
        #%%############ Problem Geometry   ##############################
        ## Mesh File
        self.meshFile = 'WellSimulator'#'MacroParallelPlates'#'WellSimulator'#
        # Geometric Values
        self.Zmin = 6
        self.ZFL = 7
        self.Zmax = 8
        self.ROut = 0.22
        self.HFluidLoss = .1
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
        # self.noSlipBCs.append('TopWall') # Both
        self.noSlipBCs.append('InnerPipe') # WelSimulator only
        self.noSlipBCs.append('OuterWall') # WelSimulator only
        self.noSlipBCs.append('BottomWall')  # Both
        #noSlipBoundaries.append('InnerWalls')
        
        ## Advected Scalar Field Inputs
        self.CInitialMixture = 0.5      # Mass fraction of Fluid0. Fluid1 = 1-Flui0
        self.scalarFieldBCs = {}
        self.COutlet = 0.8
        self.scalarFieldBCs.update({'Inlet' : Constant(self.CInitialMixture)}) # CMix
        # self.scalarFieldBCs.update({'Outlet' : Constant(self.CInitialMixture)}) # CMix: REGULAR FLOW
        self.scalarFieldBCs.update({'Outlet': Constant(self.COutlet)}) # C1: FILTRATION
        
        ## Velocity Inputs
        t=self.dtMin
        self.AOut = 2*pi*self.ROut*self.HFluidLoss
        self.rhoOut = self.rho_values[self.Fluid1] # Water Only
        self.muOut = self.mu_values[self.Fluid1] # Water Only
        self.velocityBCs = {}
        # self.VrOutlet = '0.00043 + 0*t*A*rho'
        # self.VrOutlet = 't <= 100 ? (1/(rho*A))*0.00163 : (1/(rho*A))*0.061/((pow(t,0.78)))' # BaseCase 0
        # self.VrInlet = Expression('VrInlet',VrInlet=0.0,degree=1)
        self.VrInlet = Expression('VrInlet',VrInlet=0.0,degree=2)
        self.VzOutlet = Expression('VzOutlet',VzOutlet=0.0,degree=2)
        self.VrOutlet = Expression('t <= 100 ? (1/(rho*A))*0.00163 : (1/(rho*A))*0.055 /((pow(t,0.78)))',\
                                   A=self.AOut,rho=self.rhoOut,t=t,degree=2) # All BaseCases
        # self.VrOutlet = 't <= 100 ? (1/(rho*A))*0.0163 : (1/(rho*A))*0.55 /((pow(t,0.78)))'
        # self.VrOutlet = 't <= 100 ? (1/(rho*A))*0.000163 : (1/(rho*A))*0.0055 /((pow(t,0.78)))'#'2*exp(1-(t/200))/300'#'2*exp(1-(t/200))/300'#
        # self.velocityBCs.update({'Inlet' : {2: self.VrInlet}}) # m/s
        # self.VOutlet_exp = Expression('VrOutlet',VrOutlet=self.VrOutlet,degree=2)
        self.velocityBCs.update({'Inlet' : {1: self.VrInlet}}) # m/s
        # self.velocityBCs.update({'Outlet' : {1: self.VzOutlet}}) # m/s
        self.velocityBCs.update({'Outlet' : {0: self.VzOutlet}}) # m/s
        self.velocityBCs.update({'Outlet' : {1: self.VrOutlet}}) # m/s
        
        ## Pressure Inputs
        self.pressureBCs = {}
        # Mixture Initial Density
        self.rho_mix0 = (self.rho_values[0]*self.CInitialMixture + self.rho_values[1]*(1-self.CInitialMixture))
        # Inlet Pressure
        self.pInlet = self.rho_mix0*self.Zmin*self.g #0.3164557 #self.rho_values[0]*2*self.g
        self.pressureBCs.update({'Inlet' : self.pInlet}) # Pa
        
        

        #%%############ Solver parameters ###############################
        # Absolute Tolerance    
        self.absTol = 1e-12
        
        # Relative Tolerance
        self.relTol = 1e-20
        
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