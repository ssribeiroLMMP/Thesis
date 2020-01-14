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

def Calculations(self):
    # #%%############ Fluids' Properties ##############################
    ## Tags
    self.Fluid1 = 1 - self.Fluid0 # Water
  
    # Rheology
    # Newtonian Viscosity
    self.mu_values = [self.mu_cem , self.mu_water]  
    
    # Density (kg/m³)
    # Experimental Values 
    self.rho_cem0 = (self.rho_bulk0 - (1-self.CInitialMixture)*self.rho_water)/(self.CInitialMixture)
    
    # Initial Density Values per Component
    self.rho_values = [self.rho_cem0 , self.rho_water] # kg/m³

    # Shrinkage Equation: rhoMax - ((rhoMax-rhoMin)/(1+math.exp(Inclination*(-t+t0)))+rhoMin) +rhoMin
    # if Inclination is zero, shrinkage is neglected
    self.shrinkage_rhoMin = self.rho_bulkInf    # kg/m³
    self.shrinkage_t0 = self.shrinkage_t0_ts_delay*self.ts            # s  
    # Model Equation
    self.shrinkageModel = Expression(shrinkageEquation,\
                            rhoMax=self.rho_values[self.Fluid0],rhoMin=self.shrinkage_rhoMin,Inclination = self.shrinkage_inclination, \
                            t=self.t0, t0 = self.shrinkage_t0, degree=2)
        
    # ## Pressure Inputs
    self.rhoInitialInlet = self.rho_values[0]*(self.CInitialMixture) + self.rho_values[1]*(1-self.CInitialMixture)
    self.pInlet = self.rhoInitialInlet*self.Zmin*self.g #0.3164557 #self.rho_values[0]*2*Inputs.g
    self.pressureBCs.update({'Inlet' : self.pInlet}) # Pa
    # self.pOutlet = 0.6*(self.pInlet + self.rho_values[0]*self.g*1)
    # self.pressureBCs.update({'Outlet' : self.pOutlet}) # Pa
    
    ## Advected Scalar Field Inputs
    self.scalarFieldBCs = {}
    self.scalarFieldBCs.update({'Inlet' : Constant(self.CInitialMixture)}) # CMix
    self.scalarFieldBCs.update({'Outlet': Constant(self.COutlet)}) # C1: FILTRATION
    
    ## Velocity Inputs
    t=self.dtMin
    self.AOut = 2*pi*self.ROut*self.HFluidLoss
    self.rhoOut = (1-self.COutlet)*self.rho_values[self.Fluid0]
    self.velocityBCs = {}
    # self.velocityBCs.update({'Inlet' : Expression((self.VxInlet,'0.0'),t=t,degree=1)}) # m/s
    self.velocityBCs.update({'Outlet' : Expression(('0.0',self.VrOutlet), A=self.AOut,rho=self.rhoOut,t=t,degree=2)}) # m/s
    
    return self