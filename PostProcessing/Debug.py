def newPressureInput(boundaries,Subdomains,inputs):
    xOut,yOut = coordinatesAt(boundaries,Subdomains['Outlet'])

    outletArea = 2*np.pi*(max(xOut)-min(xOut))*max(yOut)
    cumsum = 0
    n = 0
    for i in range(0,len(y)):
        cumsum = cumsum + u1(xOut[i],yOut[i])[1]
        n += 1

    flowRate = outletArea*cumsum/n
    massFlowrate = inputs.rho_values[inputs.Fluid1]*flowRate
    cInlet = inputs.CInitialMixture

    rhoMix = (1-cInlet)*inputs.rho_values[0] + cInlet*inputs.rho_values[1]
    deltaM = massFlowrate*dt

    deltaV = deltaM/rhoMix

    xIn,yIn = coordinatesAt(boundaries,Subdomains['Inlet'])

    inletArea = np.pi*(max(yIn)**2 - min(yIn)**2)

    deltaTOC = deltaV/inletArea

    pInlet = inputs.pInlet - rhoMix*inputs.g*deltaTOC