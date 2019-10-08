#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 08:47:41 2019

@author: Sergio Ribeiro
@description: Module that defines problem Equations
"""
from dolfin import *
import sys
import os
sys.path.append(os.path.abspath('..'))
from Solver.BoundaryConditions import *

# Body Forces Term: Gravity
def fb(inputs):
    # Body Forces: Gravity
    return Constant((inputs.g, 0.0))

# Rheological Model Function
def smdM(C,inputs,u,t):

    # Determine gammaDot from deformation tensor D
    D = sym(grad(u))
    gammaDot = project(sqrt(2*tr(dot(D,D))),C)

    # Time dependent Yield Stress - Curing Process: tauY(t) 
    tauY_t = inputs.tau0*exp(t/inputs.ts)

    eps = 1e-10

    # Modified smd Model (Souza Mendes e Dutra (2004)) + Cure(tauY(t))  
    smdEquation = Expression('K*(pow((abs(gammaDot)+eps),nPow)/(gammaDot+eps))',\
                            K = inputs.K,nPow = inputs.n, \
                            gammaDot = gammaDot,eps = eps, degree=2)

    # smdEquation = Expression('(1 - exp(-eta0*(gammaDot)/tauY_t))* \
    #                          (tauY_t/(gammaDot+eps)+ K*(pow((abs(gammaDot)+eps),nPow)/(gammaDot+eps))) + etaInf',\
    #                         etaInf=inputs.etaInf,eta0=inputs.eta0,K = inputs.K,nPow = inputs.n,ts = inputs.ts, \
    #                         gammaDot = gammaDot, t=t, tauY_t = tauY_t,eps = eps, degree=2)
    return project(smdEquation,C)

# Calculates Fluid properties by Mesh Cell
def assignFluidProperties(inputs,c0,C=0,u=0,t=0):
    # Constant Viscosity of Each Specie
    if C == 0:
        mu = inputs.mu_values[1]*c0 + inputs.mu_values[0]*(1-c0)
    # Cement is modeled by Modified SMD Non-Newtonian Model + Cure(tauY(t))
    else: 
        mu = smdM(C,inputs,u,t)*c0 + inputs.mu_values[0]*(1-c0)

    rho = inputs.rho_values[1]*c0 + inputs.rho_values[0]*(1-c0)
#    rho = Constant(inputs.rho_values[0])
    return rho, mu


def meshMeasures(meshObj,boundaries):
    # Define any measure associated with domain and subdomains
    dx = Measure('dx', domain=meshObj)
    ds = Measure('ds', domain=meshObj, subdomain_data=boundaries)
    
    # Vectors Normal to the Mesh
    n = FacetNormal(meshObj) # Normal vector to mesh
    return dx, ds, n

def flowSpaceCreation(inputs,meshObj):
    # Get Element Shape: Triangle, etc...
    elementShape = meshObj.ufl_cell()
    
    # Set Mesh Elements
    Uel = VectorElement(inputs.velocityElementfamily, elementShape, inputs.velocityElementOrder) # Velocity vector field
    Pel = FiniteElement(inputs.pressureElementfamily, elementShape, inputs.pressureElementOrder) # Pressure field
    UPel = MixedElement([Uel,Pel])
    
    # Mixed Function Space: Pressure and Velocity
    W = FunctionSpace(meshObj,UPel)
    
    return W

def steadyStateFlow(rho,mu,inputs,meshObj,boundaries,Subdomains):
    W = flowSpaceCreation(inputs,meshObj)    
    
    ## Result Functions
    w = Function(W)
    # Split into Velocity and Pressure
    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = W.split()
    
    # Initial Condition Function
    w0 = Function(W)
    (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
    
    # Load Important Measures: Omega, deltaOmega, Normal Vector
    dx, ds, n = meshMeasures(meshObj,boundaries)
    
    # Time step Constant
    Dt = Constant(inputs.dt)
   
    alpha = Constant(inputs.alpha)
    #%%##############################   Equations
    # Linear Momentum Conservation
         # Transient Term: dot((u - u0) / Dt, v)*dx()
         # Inertia Term             # Surface Forces Term                  # Pressure Force
    a1 = alpha*(inner(grad(u)*u , v) + (mu/rho)*inner(grad(u), grad(v)) - div(v)*p/rho)*dx() + \
    (1-alpha)*(inner(grad(u0)*u0,v) + (mu/rho)*inner(grad(u0),grad(v))- div(v)*p0/rho)*dx()  # Relaxation
                      
            # Pressure Force: Natural Boundary Conditions             # Body Forces Term: Gravity         
    L1 = 0
    for key, value in inputs.pressureBCs.items():
        Pi = Constant(value)
        L1 = L1 + (Pi/rho)*dot(v,n)*ds(Subdomains[key])
    L1 = - L1
    
    # Add Mass Conservation
    a2 = (q*div(u))*dx() 
    L2 = 0
    
    # Weak Complete Form
    F = a1 + a2 - (L1 + L2)
        
    # Jacobian Matrix
    J = derivative(F,w,dw)
    
    bcU = flowBC(U,inputs,meshObj,boundaries,Subdomains)
    
    
    #%%##############################   Numerical Solver Properties
    # Problem and Solver definitions
    problemU = NonlinearVariationalProblem(F,w,bcU,J)
    solverU = NonlinearVariationalSolver(problemU)
    # Solver Parameters
    prmU = solverU.parameters
    #info(prmU,True)  #get full info on the parameters
    prmU['nonlinear_solver'] = 'newton'
    prmU['newton_solver']['absolute_tolerance'] = inputs.absTol
    prmU['newton_solver']['relative_tolerance'] = inputs.relTol
    prmU['newton_solver']['maximum_iterations'] = inputs.maxIter
    prmU['newton_solver']['linear_solver'] = inputs.linearSolver
    
    # Solve Problem
    solverU.solve()
    
    # Append Flow Problem
    return w

#%% Transient Coupled scheeme for Flow 
def transientFlow(t,W,w0,dt,rho,mu,inputs,meshObj,boundaries,Subdomains):    
    #####  Functions and Constants
        ## Trial and Test function(s)
    dw = TrialFunction(W)
    (v, q) = TestFunctions(W)
    w = Function(W)
    
    # Split into Velocity and Pressure
    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = W.split()
    
    # Initial Conditions or previous timestep
    (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
    
    # Calculate Important Measures: Omega, deltaOmega, Normal Vector
    dx, ds, n = meshMeasures(meshObj,boundaries)
    
    # Time step Constant
    Dt = Constant(dt)
   
    alpha = Constant(inputs.alpha)
    ##########   Equations
    # Linear Momentum Conservation
          
           # Transient Term            # Inertia Term             # Surface Forces Term           # Pressure Force
    a1 = inner((u-u0)/Dt,v)*dx() + alpha*(inner(grad(u)*u , v) + (mu/rho)*inner(grad(u), grad(v)) - div(v)*p /rho)*dx() \
                            + (1-alpha)*(inner(grad(u0)*u0,v) + (mu/rho)*inner(grad(u0),grad(v)) - div(v)*p/rho)*dx()    # Relaxation
                      
    L1 = 0
    for key, value in inputs.pressureBCs.items():
        Pi = Constant(value)
               # Pressure Force: Natural Boundary Conditions
        L1 = L1 + (Pi/rho)*dot(v,n)*ds(Subdomains[key])
    # Body Forces Term: Gravity 
    L1 = - L1 + inner(fb(inputs),v)*dx()
    
    # Add Mass Conservation
    a2 = (q*div(u))*dx() 
    L2 = 0
    
    # Weak Complete Form
    F = a1 + a2 - (L1 + L2)
        
    # Jacobian Matrix
    J = derivative(F,w,dw)
    
    # Apply Flow Boundary Conditions
    bcU = flowBC(t,U,inputs,meshObj,boundaries,Subdomains)
        
    ##########   Numerical Solver Properties
    # Problem and Solver definitions
    problemU = NonlinearVariationalProblem(F,w,bcU,J)
    solverU = NonlinearVariationalSolver(problemU)
    # Solver Parameters
    prmU = solverU.parameters
    #info(prmU,True)  #get full info on the parameters
    prmU['nonlinear_solver'] = 'newton'
    prmU['newton_solver']['absolute_tolerance'] = inputs.absTol
    prmU['newton_solver']['relative_tolerance'] = inputs.relTol
    prmU['newton_solver']['maximum_iterations'] = inputs.maxIter
    prmU['newton_solver']['linear_solver'] = inputs.linearSolver
    
    # Solve Problem
    (no_iterations,converged) = solverU.solve()
    
    # Append Flow Problem
    return w,no_iterations,converged

#%% Field Advection Functions

def fieldSpaceCreation(inputs,meshObj):
    # Get Element Shape: Triangle, etc...
    elementShape = meshObj.ufl_cell()
    
    # Set Mesh Elements
    Cel = FiniteElement(inputs.scalarFieldElementfamily, elementShape, inputs.scalarFieldElementOrder) # Scalar Field
  
    # Mixed Function Space: Pressure and Velocity
    C = FunctionSpace(meshObj,Cel)
    
    return C

## Fluid Mixture
def initialConditionField(C,inputs):
    init = Expression('C0','C0',C0=inputs.CMixture,degree=2)
    c0 = Function(C)
    c0.assign(project(init,C))
    return c0

## Fluid Interface 
def initialMixture(C,inputs):
    c0 = Function(C)
    c0.assign(project(Constant(inputs.CInitialMixture),C))
    return c0

## Fluid Interface 
def initialInterface(C,inputs):
    smoothstep = Expression('(CMax-CMin)/(1+exp(IntIncl*(-x[0]+x0)))+CMin',IntIncl = 20000,x0=inputs.InterfaceX0,CMax=inputs.Fluid1,CMin=inputs.Fluid0,degree=2)
    c0 = Function(C)
    c0.assign(project(smoothstep,C))
    return c0

def transienFieldTransport(C,c0,dt,u1,D,rho,mu,inputs,meshObj,boundaries,Subdomains):
    ## Trial and Test function(s)
    c = TrialFunction(C)
    l = TestFunction(C)
    
    ## Result Functions
    c1 = Function(C)
    
    # Load Important Measures: Omega, deltaOmega, Normal Vector
    dx, ds, n = meshMeasures(meshObj,boundaries)
    
    # Time step Constant
    Dt = Constant(dt)
    alphaC = Constant(inputs.alphaC)
    
    # Concentration Equation
          # Transient Term   #                 Advection Term                         # Diffusion Term                            
    F = rho*inner((c - c0)/Dt,l)*dx() + alphaC*(rho*inner(u1,(grad(c ))*l) + D*dot(grad(c ), grad(l)))*dx() \
                                  + (1-alphaC)*(rho*inner(u1,(grad(c0))*l) + D*dot(grad(c0), grad(l)))*dx() # Relaxation
    a, L = lhs(F), rhs(F)

    # Boundary Conditions    
    bcC = fieldTransportBC(C,inputs,meshObj,boundaries,Subdomains)
    
    # Solve Problem
    solve(a == L, c1, bcC)
    
    return c1