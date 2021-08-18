# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 16:09:50 2019

@author: Yu-Chuan Chen
"""

import pandas as pd
import numpy as np
"""import all files in data folder"""
from os import walk
from os.path import join

dataPath = ['data/group/purkinje2stage_para']
neuronPath = []

for path in dataPath:
    for root, dires, files in walk(path):
        for f in files:
            fullpath = join(root, f)
            neuronPath.append(fullpath)

dataCollect = []
           
for index, fileName in enumerate(neuronPath):
    
    data = pd.DataFrame(np.load(str(fileName)), columns = ['start', 'end'])
    dataCollect.append (data)
    print(fileName)
    nameOfdata = fileName.split('.')
    nameNet = nameOfdata[0].split('\\')
    newName = 'CSV data/Purkinje_para_0604/' + nameNet[1] + '.csv'
    
    """ data output CSV file"""
    data.to_csv(newName, index=False)
    
    
    
    
    
    