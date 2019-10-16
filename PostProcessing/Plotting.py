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

def prettyplot(fig,mesh,t,ui,pi,li,dicTitle, pnlevels=10,resultspath='',tag='',cbarU=0,cbarP=0,cbarDirection = 1):
    # Mesh Vertices' Coordinates
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    nVertices = len(x)
    
    mycmap = cm.get_cmap('jet')
    
    shape = (nVertices, 2)

    
    if pi != 0:
        # Get P Values
        pValues = pi.compute_vertex_values(mesh)
        # Plot Pressure
        plt.figure(num=fig, figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')
        plt.clf()
        pax = plot(pi,title=dicTitle[1], cmap = mycmap)
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
            cbarP = plt.colorbar(pax,orientation='vertical',cmap = mycmap) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_yticklabels(pticklabels)
        else:
            cbarP = plt.colorbar(pax,orientation='horizontal',cmap = mycmap) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_xticklabels(pticklabels)   
        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Pressure_t='+'{:.2f}'.format(t) +'.png')
        
        dpdx = ((pValues[len(pValues)-1]-pValues[1])/x.max())
        
        
    # Plot Level or concentration
    if li != 0:
        # Get L Values
        cValues = li.compute_vertex_values(mesh)
        # Plot Level
        plt.figure(num=fig+1, figsize=(10, 20), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        lax = plot(li,title=dicTitle[2], cmap = mycmap)
    
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Concentration_t='+'{:.2f}'.format(t) +'.png')
        
    # Plot Velocities
    if ui != 0:
        # Get Cell Values
        uValues = ui.compute_vertex_values(mesh)
        # Plot Velocities
        plt.figure(num=fig+2, figsize=(10, 20), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        uax = plot(ui,title=dicTitle[3], cmap = mycmap)
    
        # Get Velocity Values    
        uXYValues = np.zeros(shape)    
        
        # Colect velocity data in Arrays
        for j in range(0,nVertices):
            uXYValues[j,0] = uValues[j]
            uXYValues[j,1] = uValues[j+nVertices]
        
        # Calculate Arrow Sizes
        C = np.hypot(uXYValues[:,0], uXYValues[:,1])
        # plt.axis('equal')
        minVel = '{:.5f} m/s'.format(C.min()) 
        meanVel = '{:.5f} m/s'.format(C.mean())
        maxVel = '{:.5f} m/s'.format(C.max())
        
        if cbarDirection == 1:
            cbarU = plt.colorbar(uax,orientation='vertical', cmap = mycmap) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])    
            cbarU.ax.set_yticklabels([minVel, meanVel, maxVel])
        else:
            cbarU = plt.colorbar(uax,orientation='horizontal', cmap = mycmap) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])    
            cbarU.ax.set_xticklabels([minVel, meanVel, maxVel])
        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Velocities_t='+'{:.2f}'.format(t) +'.png')
    
    
        
    return cbarU, cbarP, uXYValues, lValues, pValues, nVertices