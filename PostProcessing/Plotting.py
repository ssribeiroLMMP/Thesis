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
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
import matplotlib.ticker as ticker
from dolfin import *
import pandas as pd

class Style():
    def __init__(self):
        # Colormap
        self.colors = ['#800000',\
                        #'#A60A0F',\
                        '#ff1a1a',\
                        '#cc7a00',\
                        '#ff9900',\
                        '#ffcc80',\
                        '#ffcc00',\
                        '#9DD52D',\
                        '#4AD77E',\
                        '#4EDEDE',\
                        '#2880C3',\
                        '#99ebff']

        self.font = {'family': 'serif',
                'weight': 'normal',
                'size': 10}

def plotPressureProfileDF(inputFile):

    dataframe = pd.read_csv(inputFile) #Time(s),Depth(m),Pressure(Pa)
    fig, ax = plt.subplots()
    
    i=0

    for key, grp in dataframe.groupby(['Time']):
        timeLabel = '{:.0f}s'.format(key)
        ax = grp.plot(ax=ax, kind='line', x='Depth(m)', y='Pressure(Pa)', \
                      c = Style().colors[i], label=timeLabel)
        i+=1
    plt.legend(loc='best')
    plt.show()


def plotPressureProfile(inputFile,toutputFile):
    fig, ax = plt.subplots()
    
    dataframe = pd.read_csv(inputFile) #Time(s),P6(Pa),P6.5(Pa),P7(Pa),P7.5(Pa),P8(Pa)
    
    time = dataframe['Time(s)']

    # Z Depth series
    # 'P6(Pa)','P6.125(Pa)','P6.25(Pa)','P6.375(Pa)','P6.5(Pa)','P6.625(Pa)','P6.75(Pa)','P6.875(Pa)',\
    # 'P7(Pa)','P7.125(Pa)','P7.25(Pa)','P7.375(Pa)','P7.5(Pa)','P7.625(Pa)','P7.75(Pa)','P7.875(Pa)',\
    # 'P8(Pa)']

    x = pd.Series([6, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 6.9375,\
                   7, 7.0625,7.10, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875,  8])

    # Loop on time series
    for i in range(0,len(time)):
        y = dataframe.iloc[i,2:]
        timeLabel = '{:.0f}s'.format(round(time[i],0))
        plt.plot(x,y,c=colors[i],label = timeLabel)
        print(i)

    plt.xlim((6,8))
    plt.ylim((0,200000))
    plt.xlabel('Column Height (m)', fontdict = font)
    plt.ylabel('Pressure (Pa)', fontdict = font)
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    plt.grid(color='0.8', linestyle='-')
    plt.legend(loc = 'lower center',bbox_to_anchor=(0.5, 1.01),ncol=6,fontsize = 8)
    plt.savefig('./'+outputFile+'.png',dpi=300)

    return fig, ax

def plotResult(output,outputFile,outType = 0):
    fig, ax = plt.subplots()
    plt.clf
    im = ax.imshow(output,cmap = 'gray')
    plt.savefig('./'+outputFile+'.png',dpi=200)
    if outType ==0 :
        plt.clf
    else: 
        return fig, ax

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def prettyplot(fig,mesh,t,ui,pi,li,dicTitle, pnlevels=10,resultspath='',tag='',cbarU=0,cbarP=0,cbarDirection = 1):
    # Mesh Vertices' Coordinates
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]

    X, Y = np.meshgrid(x, y)
    
    nVertices = len(x)
    
    mycmap1 = cm.get_cmap('jet')
    mycmap2 = cm.get_cmap('autumn')
    
    shape = (nVertices, 2)

    fontSize = 20

    if pi != 0:
        # Get P Values
        pValues = pi.compute_vertex_values(mesh)
        # Plot Pressure
        fig1 = plt.figure(num=fig, figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axp = fig1.add_subplot(111)
        axp_divider = make_axes_locatable(axp)
        pax = plot(pi, cmap = mycmap1) #title=dicTitle[1],
        axp.set_aspect('auto')
        axp.set_xticks([x.min(),x.max()])
        axp.set_xticklabels(['$z_{Min}$','$z_{Max}$'],rotation=90)
        axp.set_yticks([y.min(),y.max()])
        axp.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)
        # plt.axis('equal')
        minP = pValues.min()
        meanP = pValues.mean()
        maxP = pValues.max()
        pticks = [minP]
        pticklabels = ['{:.0f} kPa'.format(minP/1000)]

        
        for l in range(1,pnlevels+1):
            levelPvalue = (l*(maxP-minP)/pnlevels) + minP
            pticks.append(levelPvalue)
            pticklabels.append('{:.0f} kPa'.format(levelPvalue/1000))
        
        caxp = axp_divider.append_axes("top", size="5%", pad="5%")

        if cbarDirection == 1:
            cbarP = plt.colorbar(pax,cax=caxp,orientation='vertical',cmap = mycmap1) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_yticklabels(pticklabels)
        else:
            cbarP = plt.colorbar(pax,cax=caxp,orientation='horizontal',cmap = mycmap1) #
            cbarP.set_ticks(pticks)
            cbarP.ax.set_xticklabels(pticklabels,rotation=90)
            caxp.xaxis.set_ticks_position("top")

        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Pressure_t='+'{:.2f}'.format(t) +'.png')
        
        dpdx = ((pValues[len(pValues)-1]-pValues[1])/x.max())
        
        
    # Plot Level or concentration
    if li != 0:
        # Get L Values
        cValues = li.compute_vertex_values(mesh)
        # Plot Level
        fig2 = plt.figure(num=fig+1, figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axc = fig2.add_subplot(111)
        axc_divider = make_axes_locatable(axc)
        lax = plot(li, cmap = mycmap2) #title=dicTitle[2]
        axc.set_aspect('auto')
        axc.set_xticks([x.min(),x.max()])
        axc.set_xticklabels(['$z_{Min}$','$z_{Max}$'],rotation=90)
        axc.set_yticks([y.min(),y.max()])
        axc.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)

        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Concentration_t='+'{:.2f}'.format(t) +'.png')
        
    # Plot Velocities
    if ui != 0:
        # Get Cell Values
        uValues = ui.compute_vertex_values(mesh)
        uXYValues = np.zeros(shape)
    
        # Colect velocity data in Arrays
        for j in range(0,nVertices):
            uXYValues[j,0] = uValues[j]
            uXYValues[j,1] = uValues[j+nVertices]

        # Plot Velocities
        fig3 = plt.figure(num=fig+2, figsize=(20, 10), dpi=180, facecolor='w', edgecolor='k')
        plt.clf()
        plt.rc('font', size=fontSize)
        axu = fig3.add_subplot(111)
        axu_divider = make_axes_locatable(axu)
        
        # Calculate Arrow Sizes
        # Get Velocity Values   
        # Initialize arrays
        # U = np.zeros_like(X); V = np.zeros_like(X); Amod = np.zeros_like(X) 
        # for i in range(0,X.shape[0]):
        #     for j in range(0,X.shape[1]):
        #         U[i][j] = ui(X[i][j],Y[i][j])[0]
        #         V[i][j] = ui(X[i][j],Y[i][j])[1]
        #         Amod[i][j] = sqrt(U[i][j]*U[i][j]+V[i][j]*V[i][j])
        
        uax = plot(ui, cmap = mycmap1) #title=dicTitle[3]
        # plt.quiver(X,Y,U,V,Amod,alpha=0.8)
        axu.set_aspect('auto')
        axu.set_xticks([x.min(),x.max()])
        axu.set_xticklabels(['$z_{Min}$','$z_{Max}$'],rotation=90)
        axu.set_yticks([y.min(),y.max()])
        axu.set_yticklabels(['$R_{In}$','$R_{Out}$'],rotation=90)

        # Calculate Arrow Sizes
        C = np.hypot(uXYValues[:,0], uXYValues[:,1])
        # plt.axis('equal')
        minVel = '{:f}'.format(C.min()) 
        meanVel = '{:f}'.format(C.mean())
        maxVel = '{:f}'.format(C.max())
        
        caxu = axu_divider.append_axes("top", size="5%", pad="10%")

        if cbarDirection == 1:
            cbarU = plt.colorbar(uax,cax=caxu,orientation='vertical', cmap = mycmap1,format=ticker.FuncFormatter(fmt)) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])    
            cbarU.ax.set_yticklabels([minVel, meanVel, maxVel])
        else:
            cbarU = plt.colorbar(uax,cax=caxu,orientation='horizontal', cmap = mycmap1,format=ticker.FuncFormatter(fmt)) #,
            cbarU.set_ticks([C.min(), C.mean(), C.max()])   
            caxu.xaxis.set_ticks_position("top")
        
        # Save Figure as .PNG file
        if resultspath != '':
            plt.savefig(resultspath+tag+'Velocities_t='+'{:.2f}'.format(t) +'.png')
    
        
        
    return cbarU, cbarP, uXYValues, cValues, pValues, nVertices