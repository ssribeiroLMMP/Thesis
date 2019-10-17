#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 08:21:14 2019

@author: Sergio Ribeiro
Description: Set of functions to plot pressure, velocity, concentration and 
             other variables
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from dolfin import *


def plotResult(output,outputFile,outType = 0):
    fig, ax = plt.subplots()
    plt.clf
    im = ax.imshow(output,cmap = 'gray')
    plt.savefig('./'+outputFile+'.png',dpi=200)
    if outType ==0 :
        plt.clf
    else: 
        return fig, ax, im

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def prettyplot(fig,mesh,t,ui,pi,li,dicTitle, pnlevels=10,resultspath='',tag='',cbarU=0,cbarP=0,cbarDirection = 1):
    # Mesh Vertices' Coordinates
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    nVertices = len(x)
    
    mycmap1 = cm.get_cmap('jet')
    mycmap2 = cm.get_cmap('autumn')
    
    shape = (nVertices, 2)

    fontSize = 20

    if pi != 0:
        # Get P Values
        pValues = pi.compute_vertex_values(mesh)
        # Plot Pressure
        fig1 = plt.figure(num=fig, figsize=(30, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axp = fig1.add_subplot(211)
        pax = plot(pi, cmap = mycmap1) #title=dicTitle[1],
        axp.set_xticks([6,7,8])
        axp.set_xticklabels(['6','7','8'],rotation=90)
        axp.set_yticks([y.min(),y.max()])
        axp.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)
        # plt.axis('equal')
        minP = pValues.min()
        meanP = pValues.mean()
        maxP = pValues.max()
        pticks = [minP]
        pticklabels = ['{:.0f} Pa'.format(minP)]
        for l in range(1,pnlevels+1):
            levelPvalue = (l*(maxP-minP)/pnlevels) + minP
            pticks.append(levelPvalue)
            pticklabels.append('{:.0f} Pa'.format(levelPvalue))
        
        if cbarDirection == 1:
            cbarP = plt.colorbar(pax,orientation='vertical',cmap = mycmap1) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_yticklabels(pticklabels)
        else:
            cbarP = plt.colorbar(pax,orientation='horizontal',cmap = mycmap1) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_xticklabels(pticklabels,rotation=90)   
        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Pressure_t='+'{:.2f}'.format(t) +'.png')
        
        dpdx = ((pValues[len(pValues)-1]-pValues[1])/x.max())
        
        
    # Plot Level or concentration
    if li != 0:
        # Get L Values
        cValues = li.compute_vertex_values(mesh)
        # Plot Level
        fig2 = plt.figure(num=fig+1, figsize=(30, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axc = fig2.add_subplot(211)
        lax = plot(li, cmap = mycmap2) #title=dicTitle[2]
        axc.set_xticks([6,7,8])
        axc.set_xticklabels(['6','7','8'],rotation=90)
        axc.set_yticks([y.min(),y.max()])
        axc.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)

        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Concentration_t='+'{:.2f}'.format(t) +'.png')
        
    # Plot Velocities
    if ui != 0:
        # Get Cell Values
        uValues = ui.compute_vertex_values(mesh)
        # Plot Velocities
        fig3 = plt.figure(num=fig+2, figsize=(30, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axu = fig3.add_subplot(211)
        uax = plot(ui, cmap = mycmap1) #title=dicTitle[3]
        axu.set_xticks([6,7,8])
        axu.set_xticklabels(['6','7','8'],rotation=90)
        axu.set_yticks([y.min(),y.max()])
        axu.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)

        # Get Velocity Values    
        uXYValues = np.zeros(shape)    
        
        # Colect velocity data in Arrays
        for j in range(0,nVertices):
            uXYValues[j,0] = uValues[j]
            uXYValues[j,1] = uValues[j+nVertices]
        
        # Calculate Arrow Sizes
        C = np.hypot(uXYValues[:,0], uXYValues[:,1])
        # plt.axis('equal')
        minVel = '{:.2e} m/s'.format(C.min()) 
        meanVel = '{:.2e} m/s'.format(C.mean())
        maxVel = '{:.2e} m/s'.format(C.max())
        
        if cbarDirection == 1:
            cbarU = plt.colorbar(uax,orientation='vertical', cmap = mycmap1) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])    
            cbarU.ax.set_yticklabels([minVel, meanVel, maxVel])
        else:
            cbarU = plt.colorbar(uax,orientation='horizontal', cmap = mycmap1) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])    
            cbarU.ax.set_xticklabels([minVel, meanVel, maxVel],rotation=90)
        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Velocities_t='+'{:.2f}'.format(t) +'.png')
    
        
        
    return cbarU, cbarP, uXYValues, cValues, pValues, nVertices