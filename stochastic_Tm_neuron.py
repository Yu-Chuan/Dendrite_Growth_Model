# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 09:44:03 2019

@author: Yu-Chuan Chen (b308)
"""
#%% library
import numpy as np
import matplotlib.pyplot as plt
import neuronclass as nc
import time
""" import all files in data folder"""
import os 
import fnmatch
import pandas as pd

#%% create save folder
path = 'D:/Yu-Chuan/ForYuChuan/python program/dedritic growth model/'
folder = 'meanTrajectory/' 
f_dir = path + folder 

try:
    os.makedirs(f_dir)
except FileExistsError:
    print("The directory has been created on %s" % f_dir)      
except OSError:
    print ("Creation of the directory %s failed" % f_dir)  
else:
    print ("Successfully created the directory %s" % f_dir)

#%% 
global n0, v0, kb, kt, diffk, sumk
v0 = 0
neuronPath = []  
allClassOfBr = []
allClassOfTe = []

neuron_type = ['Tm1', 'Tm2', 'Tm9', 'Tm20']
for T, key in enumerate(neuron_type):
    class_name = 'Neuron_two_' + key + '_movemean'
    if key == 'Tm1':
        neuron2stage = nc.Neuron_two_Tm1_movemean
    elif key == 'Tm2':
        neuron2stage = nc.Neuron_two_Tm2_movemean
    elif key == 'Tm9':
        neuron2stage = nc.Neuron_two_Tm9_movemean
    elif key == 'Tm20':
        neuron2stage = nc.Neuron_two_Tm20_movemean
    
    n0 = 10 # first assuption
    neuron = neuron2stage(n0)
    neuron.nrtraj()  
    
    neuronType = neuron.neuronType
    labname = neuron.labname
    neuronName = neuron.neuronName
    dendrite = neuron.dendrite

    workingFolder = 'Working_' + str(neuronType).capitalize()
    dataPath = ['D:/Yu-Chuan/' + workingFolder+ '/figure/NewKbAndKt/numberChange/AtRisk/Excel/' + \
                str(labname)]
      
    dataName = '*' + str(neuronType) + '_' + str(labname) + '_' + str(neuronName) + '_' + dendrite + '*.csv'

    for path in dataPath:
        for root, dires, files in os.walk(path):
            for f in files: 
                if fnmatch.fnmatch( f, dataName):
                    fullpath = os.path.join(root, f)
                    neuronPath.append(fullpath)
               
    #%% experiment data (pandas datafram)
    numberChange_exp = pd.read_csv(neuronPath[T])
    numberChange_exp.columns = ['number']

    n0 = int(numberChange_exp['number'].head(1))
    neuron.n0 = n0
    
    ### Steps setting
    nclass = 1 # to get the group of neuron 
    nsample = 500
    num_collect_data = 300
    ### inital condition
    R95=np.zeros(0)
    Rg=np.zeros(0)
    densityR95=np.zeros(0)
    densityRg=np.zeros(0)
    rterm=np.zeros(0)
    plotobj=plt.figure()
    se=np.zeros((0,2))
    se_b=np.zeros((0,2))
    se_t=np.zeros((0,2))
    groupOfneuron = []
    segment_group =[]
    segment_nsample =[]
    meanOfmultiprocess = []
    varianceOfmultiprocess = []
    NumOfrunning = []
    n_traj_collect = []
    r_traj_collect = []
    Br_collect = []
    Te_collect = []
    groupOfBr = []
    groupOfTe = []
    all_group_Br = pd.DataFrame()
    all_group_Te = pd.DataFrame()

    #time start
    start = time.time()
    for group in range(1):
    
        n_collect_1sum = np.zeros(num_collect_data)
        n_collect_2sum = np.zeros(num_collect_data)
    
        se=np.zeros((0,2))
        se_b=np.zeros((0,2))
        se_t=np.zeros((0,2))
    
        start1 = time.time()  
        print('group = '+ str(group))
        
        for isample in range(nsample):      
            #print('isample = '+ str(isample))
            
            rmax = neuron.tmax
            r_collect=np.linspace(0,rmax,num_collect_data) 
            n_collect=np.zeros(num_collect_data)    

            itraj=0
            r_collect[0]=neuron.r_traj[0]
            
            for i in range(1,num_collect_data):
                while itraj<len(neuron.r_traj) and neuron.r_traj[itraj]<r_collect[i] :
                    itraj+=1
                
                n_collect[i]=neuron.n_traj[itraj-1]
            
            f=plt.step(neuron.r_traj, neuron.n_traj, '-')
            plt.setp(f, 'color', '#cccccc', 'linewidth', 0.5)
        
            #collect trajectory
            n_traj_collect.append(neuron.n_traj)
            r_traj_collect.append(neuron.r_traj)
        
            NumOfrunning.append(neuron.n_traj.shape[0])
        
            n_collect_1sum += n_collect
            n_collect_2sum += np.multiply(n_collect,n_collect)

            se=np.append(se,neuron.startend,axis=0)
            se_b=np.append(se_b,neuron.startend_b,axis=0)
            se_t=np.append(se_t,neuron.startend_t,axis=0)
        
            # collect the br and te events
            Br_collect.append(neuron.startend_b[:, 1])
            Te_collect.append(neuron.startend_t[:, 1])
            
        #end isample
        segment_group.append(segment_nsample)
        groupOfBr.append(se_b[:, 1])
        groupOfTe.append(se_t[:, 1])
        groupOfBr_df = pd.DataFrame(se_b[:, 1], columns = ['group_' + str(group)])
        groupOfTe_df = pd.DataFrame(se_t[:, 1], columns = ['group_' + str(group)])
        all_group_Br = pd.concat([all_group_Br, groupOfBr_df], axis = 1)
        all_group_Te = pd.concat([all_group_Te, groupOfTe_df], axis = 1)
        
        # collect the mean and variance of one process
        n_collect_mean = n_collect_1sum/nsample
        n_collect_std = np.sqrt(n_collect_2sum/nsample-np.multiply(n_collect_mean, n_collect_mean))
        
        meanOfmultiprocess.append(n_collect_mean )
        varianceOfmultiprocess.append(n_collect_std)
        
        #Calculate the time of one process
        end1 = time.time()
        #print(end1-start1)  
    # end nclass 
    allClassOfBr.append(all_group_Br)
    allClassOfTe.append(all_group_Te)
    #Calculate the time of multiprocess   
    end = time.time()
    #print(end - start)
    
    n_collect_mean = n_collect_1sum/nsample
    n_collect_std=np.sqrt(n_collect_2sum/nsample-np.multiply(n_collect_mean, n_collect_mean))

    f1=plt.plot(r_collect,n_collect_mean,'-')
    plt.setp(f1, 'color', '#000000', 'linewidth', 2.0)
    f2=plt.plot(r_collect,n_collect_mean+n_collect_std, '-.',r_collect, n_collect_mean-n_collect_std, '-.')
    plt.setp(f2, 'color', '#000000', 'linewidth', 1.0)
    plt.xlabel("length of dendrite segments")
    plt.ylabel("number of dendrites")
    plt.xlim(0, 100)
    savename = neuronType + '_' + labname + '_' + str(key) + '.png'
    plt.savefig(f_dir + savename, dpi= 1500)