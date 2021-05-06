#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 14:37:11 2021

@author: abqmr
"""

import fileinput
import os
import sys


path ='/home/nathan/.local/share/kicad/5.99/footprints/HalbachGradCoils.pretty/'
footprint='4LayerGradCoils.kicad_mod'

filein=path+footprint
print(os.listdir(path))

f = open(filein,'r')
filedata = f.read()
f.close()

newdata = filedata.replace("In1.Cu","In2.Cu")

fileout=filein
f = open(fileout,'w')
f.write(newdata)
f.close()
