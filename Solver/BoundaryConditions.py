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
        
    # Velocity Conditions  #ERROR ON VERSION 1.0.4
    for key,value in inputs.velocityBCs.items():
<<<<<<< HEAD
        for key2, value2 in value.items():
            dim = key2
            velValue = value2
            bc.append(DirichletBC(U.sub(dim-1),velValue,boundariesId,subdomainsDict[key]))
=======
        dim = value[0]
        velValue = value[1]
        bc.append(DirichletBC(U.sub(dim-1),velValue,boundariesId,subdomainsDict[key]))
>>>>>>> f6494501cd5278610775d64b6616f070559a267b
    
    return bc

# Advected Field Boundary Conditions
def fieldTransportBC(C,inputs,meshId,boundariesId,subdomainsDict):
    Dim = meshId.geometric_dimension()
    
    # Initialize Boundary Condition
    bc = []
    # Scalar Field Conditions
    for key,value in inputs.scalarFieldBCs.items():
        bc.append(DirichletBC(C,value,boundariesId,subdomainsDict[key]))
    return bc
