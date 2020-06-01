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

# Shrinkage Model Function
def shrinkage(inputs,C,t):
    shrinkageModel = inputs.shrinkageModel
    # Density of assignment
    shrinkageModel.t = t
    rho_cem_t = project(inputs.shrinkageModel,C)
    
    return rho_cem_t   
# def shrinkage(inputs,C,c,t):
#     shrinkageModel = inputs.shrinkageModel
#     # Density of assignment
#     shrinkageModel.t = t
#     shrinkageModel.cFrac = c
    
#     return project(inputs.shrinkageModel,C)   

# Rheological Model Function
def smdM(inputs,C,u,t):
    
    # Determine gammaDot from deformation tensor D
    D = sym(grad(u))
    gammaDot = project(sqrt(2*tr(dot(D,D))),C)
    gammaDotArray = gammaDot.vector().get_local()
    gammaDot.vector().set_local(abs(gammaDotArray))
    # Time dependent Yield Stress - Curing Process: tauY(t) 
    inputs.tauY_t.t = t
    tauY_t = project(inputs.tauY_t,C)

    # Model variable Inputs (gammaDot, tauY(t), t)
    inputs.rheologicalModel.t=t
    inputs.rheologicalModel.gammaDot=gammaDot
    inputs.rheologicalModel.tauY_t = tauY_t

    return project(inputs.rheologicalModel,C)

def calculateNewInletPressure(TOC,massFlowrate,rho,t,dt,boundaries,Subdomains,inputs):
    # Vertices Inlet Coordinates
    xIn,yIn = coordinatesAt(boundaries,Subdomains['Inlet'])

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

    # Cement                            # Water
    # rhoMix = rho_cem_inlet*cInlet + inputs.rho_values[inputs.Fluid1]*(1-cInlet)
    
    # rhoMix = cInlet*inputs.rho_values[inputs.Fluid0] + (1-cInlet)*inputs.rho_values[inputs.Fluid1]
    
    # Mass Variation in the last timestep
    deltaM = massFlowrate*dt

    # Volume variation due to the decay of TopOfCement
    deltaV = deltaM/rhoMix

    # Calculate Inlet Area    
    inletArea = calculateAnnulusCrossArea(xIn,yIn)

    # Variation of TopOfCement
    deltaTOC = deltaV/inletArea

    # New Inlet Pressure
    TOC = TOC-deltaTOC
    pInlet = rhoMix*inputs.g*(TOC)
    
    return pInlet, TOC, rhoMix

def calculateOutletFlowrate(u1,inputs,boundaries,Subdomains):
    # Outlet Vertices Coordinates
    xOut,yOut = coordinatesAt(boundaries,Subdomains['Outlet'])
    
    outletArea = calculateCyllinderOuterArea(xOut,yOut)
    
    # Initialize flowrate
    cumsum = 0
    n = 0
    # Loop over vertices and sum the normal velocity
    for i in range(0,len(yOut)):
        # TODO: add variable velocity with outlet area normal vector
        cumsum = cumsum + u1(xOut[i],yOut[i])[1]
        n += 1
    
    # vAvg
    avgVelocity = cumsum/n

    # mDot = Q*rho : Q = vAvg.A => mDot = vAvg*A*rho
    #TODO: Insert variable outlet rho
    massFlowrate = avgVelocity*outletArea*inputs.rho_values[inputs.Fluid1]

    return massFlowrate

# Body Forces Term: Gravity
def fb(inputs):
    # Body Forces: Gravity
    return Constant((inputs.g, 0.0)) 

# # Calculates Fluid properties by Mesh Cell
# def assignFluidProperties(inputs,c0,C=0,u=0,t=0):
#     # Constant Viscosity of Each Specie
#     if C == 0:
#         mu = inputs.mu_values[inputs.Fluid0]*c0 + inputs.mu_values[inputs.Fluid1]*(1-c0)
#     # Cement is modeled by Modified SMD Non-Newtonian Model + Cure(tauY(t))
#     else: 
#         mu = smdM(C,inputs,u,t)*c0 + inputs.mu_values[inputs.Fluid1]*(1-c0)

#     rho = inputs.rho_values[inputs.Fluid0]*c0 + inputs.rho_values[inputs.Fluid1]*(1-c0)
# #    rho = Constant(inputs.rho_values[0])
#     return rho, mu

# Calculates Fluid properties by Mesh Cell
def assignFluidProperties(inputs,C,c,u,t):
    
    # Density Shrinkage
    rho_cem_t = shrinkage(inputs,C,t) 
    rho = project(Expression('rho_1*c + rho_2*(1-c)',degree=1,\
                            c=c,rho_1=rho_cem_t,rho_2=inputs.rho_water),C)
    # rho = shrinkage(inputs,C,c,t)
    
    # Rheological Model
    mu = smdM(inputs,C,u,t)
    
    return rho, mu, rho_cem_t

def meshMeasures(meshObj,boundaries):
    # Define any measure associated wshrinkageith domain and subdomains
    dx = Measure('dx', domain=meshObj)
    ds = Measure('ds', domain=meshObj, subdomain_data=boundaries)
    
    # Vectors Normal to the Mesh
    n = FacetNormal(meshObj) # Normal vector to mesh
    return dx, ds, n

# Function Space Creation - Mixed space: Flow varibles pressure and velocity
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

# def steadyStateFlow(rho,mu,inputs,meshObj,boundaries,Subdomains):
#     W = flowSpaceCreation(inputs,meshObj)    
    
#     ## Result Functions
#     w = Function(W)
#     # Split into Velocity and Pressure
#     (u, p) = (as_vector((w[0], w[1])), w[2])
#     (U, P) = W.split()
    
#     # Initial Condition Function
#     w0 = Function(W)
#     (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
    
#     # Load Important Measures: Omega, deltaOmega, Normal Vector
#     dx, ds, n = meshMeasures(meshObj,boundaries)
    
#     # Time step Constant
#     Dt = Constant(inputs.dt)
   
#     alpha = Constant(inputs.alpha)
#     #%%##############################   Equations
#     # Linear Momentum Conservation
#          # Transient Term: dot((u - u0) / Dt, v)*dx()
#          # Inertia Term             # Surface Forces Term                  # Pressure Force
#     a1 = alpha*(inner(grad(u)*u , v) + (mu/rho)*inner(grad(u), grad(v)) - div(v)*p/rho)*dx() + \
#     (1-alpha)*(inner(grad(u0)*u0,v) + (mu/rho)*inner(grad(u0),grad(v))- div(v)*p0/rho)*dx()  # Relaxation
                      
#             # Pressure Force: Natural Boundary Conditions             # Body Forces Term: Gravity         
#     L1 = 0
#     for key, value in inputs.pressureBCs.items():
#         Pi = Constant(value)
#         L1 = L1 + (Pi/rho)*dot(v,n)*ds(Subdomains[key])
#     L1 = - L1
    
#     # Add Mass Conservation
#     a2 = (q*div(u))*dx() 
#     L2 = 0
    
#     # Weak Complete Form
#     F = a1 + a2 - (L1 + L2)
        
#     # Jacobian Matrix
#     J = derivative(F,w,dw)
    
#     bcU = flowBC(U,inputs,meshObj,boundaries,Subdomains)
    
    
#     #%%##############################   Numerical Solver Properties
#     # Problem and Solver definitions
#     problemU = NonlinearVariationalProblem(F,w,bcU,J)
#     solverU = NonlinearVariationalSolver(problemU)
#     # Solver Parameters
#     prmU = solverU.parameters
#     #info(prmU,True)  #get full info on the parameters
#     prmU['nonlinear_solver'] = 'newton'
#     prmU['newton_solver']['absolute_tolerance'] = inputs.absTol
#     prmU['newton_solver']['relative_tolerance'] = inputs.relTol
#     prmU['newton_solver']['maximum_iterations'] = inputs.maxIter
#     prmU['newton_solver']['linear_solver'] = inputs.linearSolver
    
#     # Solve Problem
#     solverU.solve()
    
#     # Append Flow Problem
#     return w

#%% Deformation Tensor 
def DD(u):
    return sym(nabla_grad(u))

# Define stress tensor
def TT(u, p, mu):
    return 2*mu*DD(u) - p*Identity(len(u))

#%% Transient Coupled scheeme for Flow 
def transientFlow(t,W,w0,dt,rho,rho0,mu,inputs,meshObj,boundaries,Subdomains,Pin=0):    
    #####  Functions and Constants
        ## Trial and Test function(s)
    dw = TrialFunction(W)
    (v, q) = TestFunctions(W)
    w = Function(W)
    
    # Split into Velocity and Pressure
    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = W.split()
    
    # Initial Conditions or previous timestep
    (u0, p0) = w0.leaf_node().split()
    
    # Calculate Important Measures: Omega, deltaOmega, Normal Vector
    dx, ds, n = meshMeasures(meshObj,boundaries)
    
    # Time step Constant
    Dt = Constant(dt)
    
    alpha = Constant(inputs.alpha)
    ##########   Equations
    # Linear Momentum Conservation
    #        # Transient Term            # Inertia Term             # Surface Forces Term           # Pressure Force
    # a1 = inner((u-u0)/Dt,v)*dx() + alpha*(inner(grad(u)*u , v) + (mu/rho)*inner(grad(u), grad(v)) - div(v)*p /rho)*dx() + \
    #                            (1-alpha)*(inner(grad(u0)*u0,v) + (mu/rho)*inner(grad(u0),grad(v)) - div(v)*p /rho)*dx()    # Relaxation
            # Transient Term                # Inertia Term                           # Surface Forces Term           
    a1 = rho*dot((u-u0)/Dt,v)*dx() + alpha *(rho*dot(dot(u,nabla_grad(u)), v) + inner(TT(u,p,mu),DD(v)))*dx() + \
                                (1 - alpha)*(rho*dot(dot(u0,nabla_grad(u0)),v) + inner(TT(u,p,mu),DD(v)))*dx()  # Relaxation
                                
    L1 = 0
    for key, value in inputs.pressureBCs.items():
        if Pin>0 and key == 'Inlet':
            Pi = Constant(Pin)
        else:
            Pi = Constant(value)   
        # Pressure Force: Natural Boundary Conditions
        # L1 = L1 + (Pi/rho)*dot(v,n)*ds(Subdomains[key])
        L1 = L1 + (Pi)*dot(n,v)*ds(Subdomains[key])
    
    # Body Forces Term: Gravity 
    L1 = - L1 + inner(rho*fb(inputs),v)*dx()
    
    # Add Mass Conservation Equation
    # if t <= inputs.dt:
    #     a2 = ((dot(u,grad(rho))*q + rho*div(u)*q)*dx())
    # else:
    a2 = ((rho-rho0)/Dt)*q*dx() + (alpha)*((dot(u,grad(rho))*q + rho*div(u)*q)*dx()) + \
                                    (1-alpha)*((inner(u,grad(rho0))*q + rho0*div(u)*q)*dx())
                                
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

# Function Space Creation - Mixed space: Flow varibles pressure and velocity
def fieldSpaceCreation(inputs,meshObj):
    # Get Element Shape: Triangle, etc...
    elementShape = meshObj.ufl_cell()
    
    # Set Mesh Elements
    Cel = FiniteElement(inputs.scalarFieldElementfamily, elementShape, inputs.scalarFieldElementOrder) # Scalar Field
    # Mixed Function Space: Pressure and Velocity
    C = FunctionSpace(meshObj,Cel)
    
    return C

## Initial Conditions - 
# Expression Fluid Mixture
def initialConditionField(C,inputs):
    init = Expression('C0','C0',C0=inputs.CMixture,degree=2)
    c0 = Function(C)
    c0.assign(project(init,C))
    return c0

## Constant Fluid Mixture 
def initialMixture(C,inputs):
    c0 = Function(C)
    c0.assign(project(Constant(inputs.CInitialMixture),C))
    return c0

## Fluids smooth Interface
def initialInterface(C,inputs):
    smoothstep = Expression('(CMax-CMin)/(1+exp(IntIncl*(-x[0]+x0)))+CMin',IntIncl = 20000,x0=inputs.InterfaceX0,CMax=inputs.Fluid1,CMin=inputs.Fluid0,degree=2)
    c0 = Function(C)
    c0.assign(project(smoothstep,C))
    return c0

# Mass Transport of Component i Mass fraction
def transienFieldTransport(C,c_i0,dt,u1,D,rho_i_t,rho_i_t0,mu,inputs,meshObj,boundaries,Subdomains):
    ## Trial and Test function(s)
    c_i = TrialFunction(C)  # Mass fraction of i component
    l = TestFunction(C)   
    
    ## Result Functions
    c_i1 = Function(C)
    
    # Load Important Measures: Omega, deltaOmega, Normal Vector
    dx, ds, n = meshMeasures(meshObj,boundaries)
    
    # Time step Constant
    Dt = Constant(dt)
    alphaC = Constant(inputs.alphaC)
    
    # Concentration Equation                         
                # Transient Term   #                                                # Advection Term                                  # Diffusion Term                            
    F = inner((rho_i_t*c_i - rho_i_t0*c_i0)/Dt,l)*dx() + alphaC *(inner(u1,(grad(c_i*rho_i_t )) *l) + dot(u1, grad(c_i ))*l + c_i *div(u1)*l + (D)*dot(grad(c_i ),grad(l)))*dx() +\
                                                    (1 - alphaC)*(inner(u1,(grad(c_i0*rho_i_t0))*l) + dot(u1, grad(c_i0))*l + c_i0*div(u1)*l + (D)*dot(grad(c_i0),grad(l)))*dx() # Relaxation
    a, L = lhs(F), rhs(F)

    # Boundary Conditions    
    bcC = fieldTransportBC(C,inputs,meshObj,boundaries,Subdomains)
    
    # Solve Problem
    solve(a == L, c_i1, bcC)
    
    return c_i1

    #%% Transient Coupled scheeme for Flow 
# def transientImplicitFlow(t,W,C,w0,c0,dt,rho,mu,inputs,meshObj,boundaries,Subdomains):    
#     #####  Functions and Constants
#         ## Trial and Test function(s)
#     dw = TrialFunction(W)
#     (v, q) = TestFunctions(W)
#     w = Function(W)
#     c = TrialFunction(C)
#     l = TestFunction(C)
    
    
#     # Split into Velocity and Pressure
#     (u, p) = (as_vector((w[0], w[1])), w[2])
#     (U, P) = W.split()
#     c1 = Function(C)
    
#     # Initial Conditions or previous timestep
#     (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
    
#     # Calculate Important Measures: Omega, deltaOmega, Normal Vector
#     dx, ds, n = meshMeasures(meshObj,boundaries)
    
#     # Time step Constant
#     Dt = Constant(dt)
   
#     alpha = Constant(inputs.alpha)
#     alphaC = Constant(inputs.alphaC)
    
#     # Concentration Equation
#           # Transient Term   #                 Advection Term                         # Diffusion Term                            
#     F0 = inner((c - c0)/Dt,l)*dx() + alphaC*(inner(u,(grad(c ))*l) + (D/rho)*dot(grad(c ), grad(l)))*dx() +\
#                                     (1-alphaC)*(inner(u,(grad(c0))*l) + (D/rho)*dot(grad(c0), grad(l)))*dx() # Relaxation
#     a0, L0 = lhs(F0), rhs(F0)

#     # Boundary Conditions    
#     bcC = fieldTransportBC(C,inputs,meshObj,boundaries,Subdomains)

#     ##########   Equations
#     # Linear Momentum Conservation

#            # Transient Term            # Inertia Term             # Surface Forces Term           # Pressure Force
#     a1 = inner((u-u0)/Dt,v)*dx() \
#         + alpha*\
#             (inner(grad(u)*u , v) \
#             + ((inputs.mu_values[1]*c  + inputs.mu_values[0]*(1-c))/ \
#                (inputs.rho_values[1]*c  + inputs.rho_values[0]*(1-c)))*inner(grad(u), grad(v)) \
#             - div(v)*p / \
#                 (inputs.rho_values[1]*c  + inputs.rho_values[0]*(1-c)))*dx() \
#         + (1-alpha)* \
#             (inner(grad(u0)*u0 , v) \
#             + ((inputs.mu_values[1]*c0  + inputs.mu_values[0]*(1-c0))/ \
#                (inputs.rho_values[1]*c0  + inputs.rho_values[0]*(1-c0)))*inner(grad(u0), grad(v)) \
#             - div(v)*p / \
#                 (inputs.rho_values[1]*c0  + inputs.rho_values[0]*(1-c0)))*dx() \    
                      
#     L1 = 0
#     for key, value in inputs.pressureBCs.items():
#         Pi = Constant(value)
#                # Pressure Force: Natural Boundary Conditions
#         L1 = L1 + (Pi/(inputs.rho_values[1]*c0  + inputs.rho_values[0]*(1-c0)))*dot(v,n)*ds(Subdomains[key])
#     # Body Forces Term: Gravity 
#     L1 = - L1 + inner(fb(inputs),v)*dx()
    
#     # Add Mass Conservation
#     a2 = (q*div(u))*dx() 
#     L2 = 0
    
#     # Weak Complete Form
#     F = a0 + a1 + a2 - (L0 + L1 + L2)
        
#     # Jacobian Matrix
#     J = derivative(F,w,dw)
    
#     # Apply Flow Boundary Conditions
#     bcU = flowBC(t,U,inputs,meshObj,boundaries,Subdomains)
        
#     ##########   Numerical Solver Properties
#     # Problem and Solver definitions
#     problemU = NonlinearVariationalProblem(F,w,bcU,J)
#     solverU = NonlinearVariationalSolver(problemU)
#     # Solver Parameters
#     prmU = solverU.parameters
#     #info(prmU,True)  #get full info on the parameters
#     prmU['nonlinear_solver'] = 'newton'
#     prmU['newton_solver']['absolute_tolerance'] = inputs.absTol
#     prmU['newton_solver']['relative_tolerance'] = inputs.relTol
#     prmU['newton_solver']['maximum_iterations'] = inputs.maxIter
#     prmU['newton_solver']['linear_solver'] = inputs.linearSolver
    
#     # Solve Problem
#     (no_iterations,converged) = solverU.solve()
    
#     # Append Flow Problem
#     return w,no_iterations,converged