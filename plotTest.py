#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PostProcessing.Plotting import *
from ProblemInputs import *

def plotTest(inputs):
    dir = './PostProcessing/Cases/'+inputs.caseId
    inputFile = dir+'/pressureProfile.csv'
    outputFile = dir+'/ pressureProfile'
    # plotPressureProfile(inputFile,outputFile)
    plotPressureProfileDF(inputFile,outputFile)

if __name__ == '__main__':
    inputs = Inputs()

    # Prepare for saving
    plotTest(inputs)
    