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
import numpy as np
sys.path.append(os.path.abspath('..'))
from Solver.BoundaryConditions import *

def IP(c,l,n,meshObj):
    one = Constant(1.0)
    h = CellDiameter(meshObj)
    h_avg = (h('+') + h('-'))/2.0  
    r = one('+')*h_avg*h_avg*inner(jump(grad(c),n), jump(grad(l),n))*dS
    return (r)

# Normalization of Mollified Color Function(c = rho - Brackbird 1992)
def normalization(cFLimits):
    diffC = max(cFLimits)-min(cFLimits)
    medC = np.average(cFLimits)
    return diffC, medC

# Dirac's Delta to find Interface
def interface(c,meshObj):
    N = VectorFunctionSpace(meshObj, "CG", 2, dim=2) 
    DD = FunctionSpace(meshObj,"CG",1)
    grad_c = project(grad(c),N)
    gradientMag = sqrt(dot(grad_c,grad_c))
    gradProj = project(gradientMag,DD)
    nGamma = grad_c/gradProj
    deltaDirac = gradientMag/gradProj.vector().max()
    #nGammaMag = project(nGamma,DD)
    #nGammaNorm = nGamma/nGammaMag
    k = div(nGamma)
    
    #deltaDirac = interpolate(Expression("c > 0.32 & c < 0.68  ? 1 : 0", c=c, degree=1),V0)
    return k,nGamma,deltaDirac

# Body Forces Term: Gravity
def fb(g):
    # Body Forces: Gravity
    return Constant((0.0, -g))


def assignFluidProperties(inputs,c0):
    mu = inputs.mu_values[1]*c0 + inputs.mu_values[0]*(1-c0)
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
def transientFlow(W,w0,c0,dt,rho,mu,inputs,meshObj,boundaries,Subdomains):    
    #####  Functions and Constants
        ## Trial and Test function(s)
    dw = TrialFunction(W)
    (v, q) = TestFunctions(W)
    w = Function(W)
    
    # Split into Velocity and Pressure
    (u, p) = split(w)
    (U, P) = W.split()
    
    # Initial Conditions or previous timestep
    (u0, p0) = split(w0)
    
    # Calculate Important Measures: Omega(dx), deltaOmega(ds), Normal Vector(n)
    dx, ds, n = meshMeasures(meshObj,boundaries)
    # Rot = Expression((('0','1'),('-1','0')),degree=1)
    
    # Interface Properties
    k,nGamma,dDirac = interface(rho,meshObj)
    # Tangent of the Interface
    # tGamma = Rot* nGamma

    #tGamma = I - outer(nGamma,nGamma)
    
    # Time step Constant
    Dt = Constant(dt)
   
    alpha = Constant(inputs.alpha)
    
    ###  Equations
    # Linear Momentum Conservation
          
           # Transient Term            # Inertia Term             # Surface Forces Term           # Pressure Force              
    a1 = inner((u-u0)/Dt,v)*dx() + alpha*(inner(grad(u)*u , v) + (mu/rho)*inner(grad(u), grad(v)) - div(v)*p /rho)*dx() + \
                               (1-alpha)*(inner(grad(u0)*u0,v) + (mu/rho)*inner(grad(u0),grad(v)) - div(v)*p/rho)*dx()    # Relaxation

                      
    # Body Forces Term - Example: Weight(Gravity)
    Fb = + dot(fb(inputs.g),v)

    # Surface Tension Term
    # CSF - Brackbill 1992
    # Normalization Factors diffC = [c] = c2-c1 ; medC = <c> = (c1+c2)/2
    diffrho,medrho = normalization(inputs.rho_values)
    
    # # Interface Properties
    # k,nGamma,dDirac = interface(c0,meshObj)
    
    # Surface Tension force (but written as a Body Force)
    rhoNorm = rho/diffrho
    grad_rhoNorm = grad(rho)/medrho
    FsV = (inputs.sigma)*(k)*dot(rhoNorm*grad_rhoNorm,v)*dDirac
    # uSlip = inputs.lSlip
    
    # Viscous Drag
    # Dvwalls = - (12/(inputs.CellThickness**2))*(mu/rho)*inner(u,v)*dx()
    
    # Right Hand Side of Mommentum Conservation Equation
      # Surface Tension   # Body Forces      # Viscous Drad(2D proxy due to top and bottom walls)   
    L1 = 1/rho*FsV*dx()          # Fb*dx() +          Dvwalls          

    ## Natural Boundary Conditions
    # Pressure BCs
    for key, value in inputs.pressureBCs.items():
        Pi = Constant(value)
               # Pressure Force: Natural Boundary Conditions
        L1 = L1 - (Pi/rho)*dot(v,n)*ds(Subdomains[key])
    
    # Velocity BCs (Navier-Slip)
    for key, value in inputs.navierslip.items():
        uslip = value*dot(grad(u),n)
        L1 = L1 + (mu/rho)*value*dot(dot(grad(uslip),n),v)*ds(Subdomains[key])
    
    ## Add Mass Conservation
    # Left Hand Side of Mass Conservation Equation
    a2 = (q*div(u))*dx() 
    # Right Hand Side of Mass Conservation Equation 
    L2 = 0
    
    ## Weak Complete Form
    F = a1 + a2 - (L1 + L2)
        
    ## Jacobian Matrix
    J = derivative(F,w,dw)
    
    # Apply Flow Boundary Conditions
    bcU = flowBC(U,inputs,meshObj,boundaries,Subdomains)

    ###   Numerical Solver Properties
    # General Parameters
    # parameters['dof_ordering_library'] = 'Boost'
    # parameters['form_compiler']['quadrature_degree'] = 2
    
    ## Problem and Solver definitions
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
    init = Expression('C0','C0',C0=inputs.TInlet,degree=2)
    c0 = Function(C)
    c0.assign(project(init,C))
    return c0

# Fluid Interface 
# def initialInterface(C,inputs):
#     CIn = inputs.FluidTags[1] - 1e-10
#     COut = inputs.FluidTags[0] + 1e-10
#     IntIncl = 200
#     inflection = inputs.InterfaceX0

#     smoothstep = Expression('(CMax-CMin)/(1+exp(IntIncl*(-x[0]+x0)))+CMin',CMin=CIn,CMax=COut, IntIncl=IntIncl,x0=inflection,degree=1)
#     # smoothstep = Expression('(CMin*exp(IntIncl*x0)+CMax*exp(IntIncl*x[0]))/(exp(IntIncl*x0)+exp(IntIncl*x[0]))',IntIncl = 300,x0=inputs.InterfaceX0,CMin=CIn,CMax=COut,degree=2)
#     # smoothstep = Expression('x[0] < x0 ? 0 : 1',x0=inputs.VxInlet,degree=1)
#     c0 = Function(C)
#     c0.assign(project(smoothstep,C))
#     return c0
def initialInterface(C,inputs):
    CIn = inputs.FluidTags[1] - 1e-10
    COut = inputs.FluidTags[0] + 1e-10
    IntIncl = 500
    x0 = inputs.InterfaceX0
    y0 = inputs.InterfaceY0
    R = inputs.InterfaceR
    # (Cmax-Cmin) / (1+exp(IntIncl*(-sqrt(  (x[0]-x0)^2  +  (x[1]-y0)^2)      +R)))+Cmin
    
    smoothstep = Expression('(CMax-CMin) / (1 + exp(IntIncl* (-pow(pow(x[0]-x0,2) + pow(x[1]-y0,2),0.5)+R)))+CMin',CMin=CIn,CMax=COut, IntIncl=IntIncl,x0=x0,y0=y0,R=R, degree=2)
    # smoothstep = Expression('(CMin*exp(IntIncl*x0)+CMax*exp(IntIncl*x[0]))/(exp(IntIncl*x0)+exp(IntIncl*x[0]))',IntIncl = 300,x0=inputs.InterfaceX0,CMin=CIn,CMax=COut,degree=2)
    # smoothstep = Expression('x[0] < x0 ? 0 : 1',x0=inputs.VxInlet,degree=1)
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
          # Transient Term   #                 Advection Term                         # Diffusion Term        # Penalty function: #\#+ IP(c,l,n,meshObj) \
    F = rho*inner((c - c0)/Dt,l)*dx() + alphaC*(rho*inner(u1,(grad(c ))*l) + D*dot(grad(c ), grad(l)))*dx()  \
        + (1-alphaC)*(rho*inner(u1,(grad(c0))*l) + D*dot(grad(c0), grad(l)))*dx() # Relaxation
    a, L = lhs(F), rhs(F)

    # Boundary Conditions    
    bcC = fieldTransportBC(C,inputs,meshObj,boundaries,Subdomains)
    
    # Solve Problem
    solve(a == L, c1, bcC)
    
    return c1

def simpleReinit(C, c1, inputs, IntIncl = 200, inflection = 0.5):
    CIn = inputs.FluidTags[1] - 1e-10
    COut = inputs.FluidTags[0] + 1e-10

    # sharpener = Expression('c < 0.5 ? 0 : 1',c=c1,degree=1)
    # sharpener = Expression('(CMin*exp(IntIncl*x0)+CMax*exp(IntIncl*x[0]))/(exp(IntIncl*x0)+exp(IntIncl*c))',IntIncl = IntIncl,c=c1, x0=inflection,CMin=CIn,CMax=COut,degree=2)
    sharpener = Expression('(CMax-CMin)/(1+exp(IntIncl*(-cAdv+x0)))+CMin',CMin=CIn,CMax=COut, IntIncl=IntIncl,cAdv=c1,x0=inflection,degree=1)
    c0 = project(sharpener,C)
    return c0
    
