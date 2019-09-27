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
import meshio

def msh2hdf5(inPath,inFile,outPath,outFile, option=1):
    if option ==1:
        # Read .msh file
        msh = meshio.read(inPath+inFile)
        
        # Write h5+xdmf files
        for key in msh.cell_data:
            # Write Physical Lines
            if key == "line":
                meshio.write(Wri_path+"mf.xdmf",meshio.Mesh(points=msh.points,cells={"line": msh.cells["line"]},cell_data={"line": {"/subdomains": msh.cell_data["line"]["gmsh:physical"]}}))
            # Write Cells
            else:
                meshio.write(outPath+outFile+".xdmf",meshio.Mesh(points=msh.points,cells={key: msh.cells[key]}))
    else:
        cmd = 'meshio-convert '+inPath+inFile+'.msh '+outPath+outFile+'.xdmf'
        os.system(cmd)
        meshObj,boundaries,markers = applyMesh(1,inPath, outPath,outFile)
        write_hdf5(meshpath+File,meshObj,boundaries,markers)
        
#    mesh_file = HDF5File(MPI.comm_world, File+'.h5','w')
#    mesh_file.write(meshObj, '/mesh')
#    mesh_file.write(markers, '/subdomains')
#    mesh_file.write(boundaries, '/boundaries')

# Read Mesh from .xml 3 Files
def read_xml(File):
    # Read Mesh Object
    meshObj = Mesh(File+'.xml')        
            
    # Initialize boundaries (boundaries facets)
    boundaries = MeshFunction('size_t',meshObj, File + "_facet_region.xml")
    
    # Initialize subdomain (inlet,outlet, etc...)
    markers = MeshFunction('size_t',meshObj, File + '_physical_region.xml')
    
    return meshObj, markers, boundaries

# Write Mesh into .h5 and .xdmf pair of Files
def write_hdf5(File):
    meshObj, markers, boundaries = read_xml(File)
    hdfw = HDF5File(self.mesh.mpi_comm(), File+'.h5', "w")
    hdfw.write(meshObj, "/mesh")
    hdfw.write(markers, "/markers")
    hdfw.write(boundaries, '/boundaries')
    hdfw.close()

# Read Mesh from .h5 and .xdmf pair of Files
def read_hdf5(File):
    hdfr = HDF5File(mesh.mpi_comm(), File, 'r')
    
    hdf.read(meshObj, '/mesh', False)
    markers = MeshFunction('size_t', meshObj, meshObj.topology().dim)
    boundaries = MeshFunction('size_t', meshObj, meshObj.topology().dim -1) 

    hdf.read(markers, '/markers')
    hdf.read(boundaries, '/boundaries')
    return meshObj, markers, boundaries

    
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
 
def applyMesh(Type,geopath,meshpath,meshFile):
    # Load Subdomains
    if Type < 3:
        Subdomains = readDomains(geopath,meshFile)
    
    # Read from XML Data
    if Type == 1: 
        # Option 1 - Gmsh generation
        meshObj, markers, boundaries = read_xml(meshpath+ meshFile)
    
    elif Type == 2:
        # Read Parallelly from H5 Data
        meshObj, markers, boundaries = read_hdf5(File)
    
    elif Type == 3:
        lx = geopath[0]
        ly = geopath[1]
        nx = meshpath[0]
        ny = meshpath[1]
        method = meshFile
        
        # Method inside must be defined for all SubdDomain Classes
        class TopWall(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1],ly)
        
        class BottomWall(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1],0)
            
        class Inlet(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],0)
            
        class Outlet(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],lx)
        
        meshObj = RectangleMesh(Point(0,0),Point(lx,ly),nx,ny,method)
        markers = MeshFunction("size_t", meshObj, meshObj.topology().dim())
        boundaries = MeshFunction("size_t", meshObj, meshObj.topology().dim() - 1)
        
        Subdomains = {}
        boundaries.set_all(0)
        
        Subdomains['Inlet'] = 1
        inlet = Inlet(); inlet.mark(boundaries,Subdomains['Inlet']); 
        Subdomains['Outlet'] = 2
        outlet = Outlet(); outlet.mark(boundaries, Subdomains['Outlet']); 
        Subdomains['TopWall'] = 3 
        topWall = TopWall(); topWall.mark(boundaries, Subdomains['TopWall']); 
        Subdomains['BottomWall'] = 4
        bottomWall = BottomWall(); bottomWall.mark(boundaries, Subdomains['BottomWall']);
        Subdomains['Fluid'] = 5
        markers.set_all(Subdomains['Fluid'])
        
        
    return meshObj, boundaries, markers, Subdomains            
     
    
meshFile = 'OriginalHeleShaw'
geopath = './Geometry/'
meshpath = './Mesh/'
#msh2xml(geopath,meshFile,meshpath,meshFile)
##msh2xdmf(geopath,meshFile,meshpath,meshFile)
#
#
#lx = 0.770;
#ly = 0.15;
#nx = int(0.42*400);
#ny = int(0.42*150);
#num_el = nx*ny*4
#meshObj, boundaries, markers, Subdomains = applyMesh(1,geopath,meshpath,meshFile)
#plot(meshObj)
#meshObj, boundaries, markers, Subdomains = applyMesh (3,[lx,ly],[nx,ny],'crossed')


