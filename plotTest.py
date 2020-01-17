#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PostProcessing.Plotting import *


def plotPressureProfile(CaseId):
    inputFile = '/home/sergio/Documents/Thesis/GitHub/Master/Thesis/PostProcessing/Cases/'+CaseId+'/pressureProfile.csv'
    outputFile = '/home/sergio/Documents/Thesis/GitHub/Master/Thesis/PostProcessing/Cases/'+CaseId+'/Images/PressureProfile'

    # plotPressureProfile(inputFile,outputFile)

    plotPressureProfileDF(inputFile,outputFile)

# Calls Plot Functions
if __name__ == "__main__":
    CaseId = sys.argv[1]
    plotPressureProfile(CaseId)