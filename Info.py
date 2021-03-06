#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 16:00:58 2019

@author: adminlinux
"""
import datetime
from ProblemInputs import *
from Solver.main import *
from PostProcessing.Saving import *
from PostProcessing.Plotting import *

def infoCommand():    
    return 'python3 ./Info.py>./PostProcessing/Cases/'+inputs.caseId+'/Execution.txt' 

def solveCommand():
    return 'python3 ./Solver/main.py>>./PostProcessing/Cases/'+inputs.caseId+'/Execution.txt'

   
inputs = Inputs()
mainPaths = Paths()

# Store Simulation Execution Info
print('Simulation Execution: '+str(datetime.datetime.now()))
print('Case: '+inputs.caseId)
print('Source: '+os.path.basename(os.getcwd()))
    
inputsDict = inputs.__dict__
for key, value in inputsDict.items():
#    String = 'type('+key+')==str or type('+key+')==float or type('+key+')==int'
#    if eval(String) and not(key.startswith('_') or key.startswith('stored_')or key.startswith('RTLD_')):
    print('{0}={1}'.format(key,value))

print('\n')