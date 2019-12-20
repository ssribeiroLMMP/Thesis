#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 08:21:14 2019

@author: Sergio Ribeiro
Description: Module that processes all Boundary and Initial Conditions
"""
from dolfin import *
import numpy as np

#%%############	   Get new Velocity Boundary Conditions for Velocity
def calculareNewOutletFlowrate(initialP, velocityBCs,inputs):
    dPdL = (inputs. POutlet - initialP)/inputs.PorousRadius
    VrOutlet = '-(FormationK*dPdL/mu)'
    velocityBCs.update({'Outlet' : Expression(('0.0',VrOutlet), mu = inputs.mu_values[inputs.COutlet],dPdL=dPdL,\
                                                                FormationK = inputs.FormationK, degree=2)}) # m/s
    return velocityBCs

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

def calculateCyllinderOuterArea(xOut,yOut):
    # Outlet cross-section Area
    outletArea = 2*np.pi*(max(xOut)-min(xOut))*max(yOut)
    return outletArea


def calculateAnnulusCrossArea(xIn,yIn):
    # Inlet cross-section area
    inletArea = np.pi*(max(yIn)**2 - min(yIn)**2)
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
def flowBC(t,U,inputs,meshId,boundariesId,subdomainsDict):
    Dim = meshId.geometric_dimension()
    
    noSlipU = noSlip(Dim)
    
    # Outlet Vertices Coordinates
    xOut,yOut = coordinatesAt(boundariesId,subdomainsDict['Outlet'])
    
    # Outlet cross-section Area
    outletArea = 2*np.pi*(max(xOut)-min(xOut))*max(yOut)

    # OPnly Fluid 0 leaves through outlet
    rhoOut = inputs.rho_values[inputs.Fluid1]

    # Initialize Boundary Condition
    bc = []
    
    ## Dirichlet Conditions
    # No-slip Condition
    for i in range(0,len(inputs.noSlipBCs)):
        bc.append(DirichletBC(U,noSlipU,boundariesId,subdomainsDict[inputs.noSlipBCs[i]]))
    
    # Update Velocity Boundary Conditions
    inputs.velocityBCs = calculareNewOutletFlowrate(initialP, inputs.velocityBCs,inputs)

    # Velocity Conditions  #ERROR ON VERSION 1.0.4
    for DomainKey,valueExp in inputs.velocityBCs.items():
        if t < 4:
            valueExp.t = 4
        else:
            valueExp.t = t
        valueExp.A = outletArea
        valueExp.rho = rhoOut
        bc.append(DirichletBC(U,valueExp,boundariesId,subdomainsDict[DomainKey]))   
    
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
