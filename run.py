#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:07:36 2019

@author: sergio
@description: calls main.py execution
"""
from Info import *

createDirectories(inputs)

# Overwrite the previous Execution.txt with new Header
os.system(infoCommand())

# # Append Execut                   vi  on.txt with Simulation Residues
os.system(solveCommand())
# main(inputs)
os.system(plotCommand())