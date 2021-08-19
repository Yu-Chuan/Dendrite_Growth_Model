# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:13:56 2019

@author: Yu-Chuan Chen
"""
#%% library
import numpy as np
import os
import matplotlib.pyplot as plt

#%% create folder
dataPath = '/Yu-Chuan/ForYuChuan/python program/dedritic growth model/meanOfACF/'
neuronType = 'Tm'
labname = 'LeeCH'
neuron = ['Tm1', 'Tm2', 'Tm9', 'Tm20']
dendrite = 'Dendrite'
Segment_name = ['Br', 'Te']
folder = 'figure/'
f_dir = dataPath + folder

def folderCreate(f_dir):
    try:
        os.makedirs(f_dir)
    except FileExistsError:
        print("The directory has been created on %s" % f_dir)      
    except OSError:
        print ("Creation of the directory %s failed" % f_dir)  
    else:
        print ("Successfully created the directory %s" % f_dir)
    return

folderCreate(f_dir)

#%% Running program
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
    
#%% Figure output
##plot kb
savename = 'meanOfACF_' + labname + '_' + neuronType + '_' + Segment_name[0] + '.png'
            
plt.figure(figsize=(10,6), linewidth = 1.5)
for neuron_n in meanOfACF:
    plt.plot(neuron_n[0])
plt.legend(neuron)   
plt.xlabel("lags", fontsize = 18)
plt.ylabel("ACF", fontsize = 18)
plt.xlim([-0.5, 60])
plt.savefig(f_dir + savename)

##plot kt
savename = 'meanOfACF_' + labname + '_' + neuronType + '_' + Segment_name[1] + '.png'
plt.figure(figsize=(10,6), linewidth = 1.5)
for neuron_n in meanOfACF:
    plt.plot(neuron_n[1])
plt.legend(neuron)   
plt.xlabel("lags", fontsize = 18)
plt.xlim([-0.5, 120])
plt.ylabel("ACF", fontsize = 18)
plt.savefig(f_dir + savename)