#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 08:21:14 2019

@author: Sergio Ribeiro
Description: Module that processes all Boundary and Initial Conditions
"""
from dolfin import *

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
def flowBC(U,inputs,meshId,boundariesId,subdomainsDict):
    Dim = meshId.geometric_dimension()
    
    noSlipU = noSlip(Dim)
    
    # Initialize Boundary Condition
    bc = []
    
    ## Dirichlet Conditions
    # No-slip Condition
    for i in range(0,len(inputs.noSlipBCs)):
        bc.append(DirichletBC(U,noSlipU,boundariesId,subdomainsDict[inputs.noSlipBCs[i]]))
    
    return bc
    # Velocity Condition
    for key,value in inputs.velocityBCs.items:
        Vi = interpolate(U,value)
        bc.append(DirichletBC(U,Vi,boundariesId,subdomainsDict[key]))
    
    return bc
# Concentration Boundary Conditions
#def concentrationBC(C,meshId,boundariesId,Subdomains):
#    Dim = meshId.geometric_dimension()
#    
#    # Initialize Boundary Condition
#    bc = []
#    # Dirichlet Conditions
#    bc.append(DirichletBC(C,InletConcentration,boundariesId,Subdomains['Inlet']))
#    bc.append(DirichletBC(C,OutletConcentration,boundariesId,Subdomains['Outlet']))
#    
#    return bc
