# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:13:56 2019

@author: Yu-Chuan Chen
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
neuronType = 'Tm'
labname = 'LeeCH'
neuron = ['Tm1', 'Tm2', 'Tm9', 'Tm20']
dendrite = 'Dendrite'
Segment_name = ['Br', 'Te']

dataPath = '/Yu-Chuan/ForYuChuan/python program/dedritic growth model/meanOfACF/'
neuronPath = []
meanOfACF = []
for name in neuron:
    neuron_filesName = []
    data_kb_kt = []
    for Prop in Segment_name:
        fileName = dataPath + 'MeanOfACF_' + neuronType + '_' + labname + '_' + name  + '_' +\
                dendrite + '_' + Prop + '.csv'
        data_kb_kt.append(np.genfromtxt(fileName, delimiter=',')) 
        neuron_filesName.append(fileName)
            
    neuronPath.append(neuron_filesName)
    meanOfACF.append(data_kb_kt)
    
#%%
for neuron_n in meanOfACF:
    ##plot kb
    plt.plot(neuron_n[0])
plt.legend(neuron)   
plt.xlabel("lags")
plt.ylabel("autocorrelation function")
plt.xlim([-0.5, 60])
plt.show()

for neuron_n in meanOfACF:
    ##plot kt
    plt.plot(neuron_n[1])
plt.legend(neuron)   
plt.xlabel("lags")
plt.xlim([-0.5, 120])
plt.ylabel("autocorrelation function")
plt.show()