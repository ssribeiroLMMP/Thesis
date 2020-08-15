#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 08:21:14 2019

@author: Sergio Ribeiro
Description: Module that processes all Boundary and Initial Conditions
"""
from dolfin import *
import numpy as np

#%%############	   Coordinates at certain Subdomain
def coordinatesAt(boundaries,SubdomainVal):
    SubdomainVertices = SubsetIterator(boundaries,SubdomainVal)
    x = []
    y = []
    for f in SubdomainVertices:
        for v in vertices(f):
            x.append(v.point().x())
            y.append(v.point().y())

    x.sort()
    y.sort()

    return x,y

#RZ
def calculateCyllinderOuterArea(xOut,yOut):
    # Outlet cross-section Area: 2*pi*ROut(max(Zfl) - min(Zfl))
    outletArea = 2*np.pi*max(xOut)*(max(yOut)-min(yOut))
    return outletArea

#RZ
def calculateAnnulusCrossArea(xIn,yIn):
    # Inlet cross-section area: pi*(Rout²-Rin²)
    inletArea = np.pi*(max(xIn)**2 - min(xIn)**2)
    return inletArea
    
#%%############     Initial Conditions
# Flow
# None Given

#%%############	    Boundary Conditions
def noSlip(Dim):
# Impermeable and No-slip Condition
    noSlipU = [] # m/s
    for i in range(0,Dim):
        noSlipU.append(0.0)
    
    return Constant(noSlipU)

# Flow Boundary Conditions
def flowBC(t,rho,U,inputs,meshId,boundariesId,subdomainsDict):
    Dim = meshId.geometric_dimension()
    
    noSlipU = noSlip(Dim)
    
    # Outlet Vertices Coordinates
    xOut,yOut = coordinatesAt(boundariesId,subdomainsDict['Outlet'])
    
    # Outlet cross-section Area
    outletArea = calculateCyllinderOuterArea(xOut,yOut)

    # OPnly Fluid 0 leaves through outlet
    rhoOut = inputs.rho_values[inputs.Fluid1]

    # Vertices Inlet Coordinates
    xIn,yIn = coordinatesAt(boundariesId,subdomainsDict['Inlet'])

    # Concentration at the inlet
    # TODO: Add temporal variant cInlet
    cInlet = inputs.CInitialMixture

    #Mixture density at the inlet

    # Loop over vertices and sum the rho_cem_inlet
    cumsum = 0
    n = 0
    for i in range(0,len(yIn)):
        # Local Cement Density by Shrinkage Model
        rho_Cem_i = rho(xIn[i],yIn[i])
        cumsum = cumsum + rho_Cem_i
        n += 1
    
    # avg Inlet Cement Density
    rhoMix = cumsum/n  
    
    # Initialize Boundary Condition
    bc = []
    
    ## Dirichlet Conditions
    # No-slip Condition
    for i in range(0,len(inputs.noSlipBCs)):
        bc.append(DirichletBC(U,noSlipU,boundariesId,subdomainsDict[inputs.noSlipBCs[i]]))
    
    
    # Velocity Conditions  #ERROR ON VERSION 1.0.4
    for DomainKey,value in inputs.velocityBCs.items():
        # valueExp.A = outletArea
        # valueExp.rho = rhoOut
        for key2, valueExp in value.items():
            if t < 4:
                valueExp.t = 4
            else:
                valueExp.t = t
            dim = key2
            if DomainKey == 'Inlet' and key2 == 1:
                valueExp.rho = rhoMix
            bc.append(DirichletBC(U.sub(dim),valueExp,boundariesId,subdomainsDict[DomainKey]))
    
    return bc

# Advected Field Boundary Conditions
def fieldTransportBC(C,inputs,meshId,boundariesId,subdomainsDict):
    Dim = meshId.geometric_dimension()
    
    # Initialize Boundary Condition
    bc = []

    # Scalar Field Conditions
    # Outlet Vertices Coordinates
    xOut,yOut = coordinatesAt(boundariesId,subdomainsDict['Outlet'])
    outletArea = calculateCyllinderOuterArea(xOut,yOut)

    for DomainKey,valueExp in inputs.scalarFieldBCs.items():
        valueExp.rho = inputs.rho_values[1]
        valueExp.Area = outletArea
        bc.append(DirichletBC(C,valueExp,boundariesId,subdomainsDict[DomainKey]))
    return bc
