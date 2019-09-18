#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Set of functions to convert mesh files to run on FEniCs
"""
Created on Fri Jul 26 13:04:40 2019
Version: 1.0
@author: Sergio Ribeiro
"""
import os
from dolfin import *


### Function Definition - Converts .msh to .xml
# no formats should be passed
def msh2xml(inPath,inFile,outPath,outFile):
    
    cmd = 'dolfin-convert '+inPath+inFile+'.msh '+outPath+outFile+'.xml'
    os.system(cmd)
    
### Function Definition - Converts .msh to .xml
# Example of msh file header to Physical Names
#$PhysicalNames
#6
#1 1 "Inlet"
#1 2 "Outlet"
#1 3 "InnerPipe"
#1 4 "OuterWall"
#1 5 "BottomWall"
#2 6 "Fluid"
#$EndPhysicalNames  
#
# no formats should be passed
def readDomains(inPath,inFile):
    # Read .msh File
    fid = open(inPath+inFile+'.msh', 'r')
    
    # Initialize variables
    found = 0
    finished = 0
    physicalNames = {}
    
    
    # Loop througn .msh lines
    for line in fid:
        if '$EndPhysicalNames' in line:
            finished == 1
            break
        elif '$PhysicalNames' in line:
            found = 1
        elif found==1 and finished == 0:
            word=line.split()
            if len(word)==3:
                physicalNames[word[2][1:len(word[2])-1]] = int(word[1])
     
    return physicalNames
 
def applyMesh(Type,geopath, meshpath,meshFile):
    # Load Subdomains
    Subdomains = readDomains(geopath,meshFile)
    
    if Type == 1:
        # Option 1 - Gmsh generation
        meshObj = Mesh(meshpath+meshFile +'.xml')
        
        # Initialize boundaries (inlet, outlet, etc...)
        boundaries = MeshFunction('size_t',meshObj,meshpath+ meshFile + "_facet_region.xml")
        
        # Initialize subdomain (fluid)
        markers = MeshFunction('size_t',meshObj,meshpath+ meshFile + '_physical_region.xml')

        
    return meshObj, boundaries, markers, Subdomains            
   
    
File = 'WellSimulator'
msh2xml('./Geometry/',File,'./Mesh/',File)
    
