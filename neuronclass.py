#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:16:40 2019

@author: cherri, Yu-Chuan Chen

"""
#%% Package import
import numpy as np
import random as rand
import pandas as pd
"""import all files in data folder"""
import os
import fnmatch

#%% class defination
class Neuron:
    "this is a class for one neuron randomly growing"
    def __init__(self,n0):
        self.n0 = n0
        self.rmax = 0.0
        self.tmax = 0.0
        self.neuronType = 'neuronType'
        self.labname = 'labname'
        self.neuronName = 'neuronName'
        self.dendrite = 'dendrite'
        
    def statistics(self):
        self.dendritelength=self.startend[:,1]-self.startend[:,0]
        self.totallength=np.sum(self.dendritelength)
        self.R_terminating=self.r_terminating**0.5
        self.R95=np.percentile(self.r_terminating**0.5,95)
        self.Rg=np.sqrt(np.sum(self.n_traj/2*(np.append(np.diff(self.r_traj**2),0))/self.totallength))
        self.densityR95=self.totallength/np.pi/self.R95**2
        self.densityRg=self.totallength/np.pi/self.Rg**2
        return

class Neuron_one(Neuron):
    type='one-stage growing neoron'
    
    def nrtraj(self,kb,kt):
        diffk=kb-kt
        rmax=1/abs(diffk)*10
        r=0.
        n=self.n0
        self.n_traj=np.array([n])
        self.r_traj=np.array([r])
        self.dr=np.zeros(0) #[0]
        self.r_terminating=np.zeros(0) #[]
        #n_terminating=0
        n0=self.n0
        self.startend=np.zeros((n0,2)) #keep the starting and ending r for each dendrite
        self.startend_b=np.zeros((0,2))
        self.startend_t=np.zeros((0,2))
        surviving=np.arange(n0)
        nindexing=n0
        #Gillesipe alogrithm
        while n>0 and r<rmax:
            prop=[n*kb, n*kt]
            a=sum(prop)
            deltar=-1./a*np.log(rand.random())
            r=r+deltar
            self.dr=np.append(self.dr,deltar)#dr.append(deltar)
            if rand.random()<prop[0]/a:
                i_branch=rand.randrange(n)
                n=n+1 # branching event
                nindexing+=2
                branch=surviving[i_branch]
                surviving=np.delete(surviving,i_branch)
                self.startend[branch,1]=r
                self.startend_b=np.append(self.startend_b, [self.startend[branch]],axis=0)
                self.startend=np.append(self.startend,[[r,0]],axis=0)
                self.startend=np.append(self.startend,[[r,0]],axis=0)
                surviving=np.append(surviving, [nindexing-2, nindexing-1])
            else:
                i_term=rand.randrange(n) #index of the terminating dendrite
                n=n-1 # terminating event
                self.r_terminating=np.append(self.r_terminating,r)#r_terminating.append(r)
                term=surviving[i_term]
                surviving=np.delete(surviving,i_term)
                self.startend[term,1]=r
                self.startend_t=np.append(self.startend_t, [self.startend[term]],axis=0)
                    #n_terminating+=1
#            if n!=0:
            self.n_traj=np.append(self.n_traj,n)
            self.r_traj=np.append(self.r_traj,r)
        #end while loop here
        return
     
 
class Neuron_two(Neuron):
    type='two-stage growing neoron'
    
    def nrtraj(self,kb_in,kt_in,t):
        diffk=kb_in[-1]-kt_in[-1]
        rmax=1/abs(diffk)*10+t
        r=0.
        n=self.n0
        self.n_traj=np.array([n])
        self.r_traj=np.array([r])
        self.dr=np.zeros(0) #[0]
        self.r_terminating=np.zeros(0) #[]
        #n_terminating=0
        n0=self.n0
        self.startend=np.zeros((n0,2)) #keep the starting and ending r for each dendrite
        self.startend_b=np.zeros((0,2))
        self.startend_t=np.zeros((0,2))
        surviving=np.arange(n0)
        nindexing=n0
        while n>0 and r<rmax:
            if r<t:
                kb=kb_in[0]
                kt=kt_in[0]
            else:
                kb=kb_in[1]
                kt=kt_in[1]
                
            prop=[n*kb, n*kt]
            a=sum(prop)
            deltar=-1./a*np.log(rand.random())
            r=r+deltar
            self.dr=np.append(self.dr,deltar)#dr.append(deltar)
            if rand.random()<prop[0]/a:
                i_branch=rand.randrange(n)
                n=n+1 # branching event
                nindexing+=2
                branch=surviving[i_branch]
                surviving=np.delete(surviving,i_branch)
                self.startend[branch,1]=r
                self.startend_b=np.append(self.startend_b, [self.startend[branch]],axis=0)
                self.startend=np.append(self.startend,[[r,0]],axis=0)
                self.startend=np.append(self.startend,[[r,0]],axis=0)
                surviving=np.append(surviving, [nindexing-2, nindexing-1])
            else:
                i_term=rand.randrange(n) #index of the terminating dendrite
                n=n-1 # terminating event
                self.r_terminating=np.append(self.r_terminating,r)#r_terminating.append(r)
                term=surviving[i_term]
                surviving=np.delete(surviving,i_term)
                self.startend[term,1]=r
                self.startend_t=np.append(self.startend_t, [self.startend[term]],axis=0)
                    #n_terminating+=1
#            if n!=0:
            self.n_traj=np.append(self.n_traj,n)
            self.r_traj=np.append(self.r_traj,r)
        #end while loop here
        return
        

#%% Moving average function
def movingMeanParam(probability, neuronType, labname, neuronName, dendrite):
    workingFolder = 'Working_' + str(neuronType).capitalize()
    dataPath = ['D:/Yu-Chuan/' + workingFolder+ '/figure/NewKbAndKt/matrixOutput/' + \
                str(labname)]
    neuronPath = []    
    dataName = str(probability) + '*' +str(neuronType) + '_' + str(labname) + '_' + str(neuronName) +\
                '_' + str(dendrite) + '.csv'

    for path in dataPath:
        for root, dires, files in os.walk(path):
            for f in files: 
                if fnmatch.fnmatch( f, dataName):
                    fullpath = os.path.join(root, f)
                    neuronPath.append(fullpath)

    if str(probability) == 'KB':
         KB_param = pd.read_csv(neuronPath[0])
         KB_param.columns = ['KB']
         parameterSet = KB_param
    elif str(probability) == 'KT':
        KT_param = pd.read_csv(neuronPath[0])
        KT_param.columns = ['KT']
        parameterSet = KT_param
    else:
        print('Input have to be KB or KT')
        
    return parameterSet

#%%
class Neuron_two_Tm1_movemean(Neuron):
    type='two-stage growing neoron'
    
    def nrtraj(self):
        #diffk = kb_in - kt_in
        #rmax=1/abs(diffk)*10
        self.neuronType = 'Tm'
        self.labname = 'LeeCH'
        self.neuronName = 'Tm1'
        self.dendrite = 'Dendrite'

        r=0.
        n=self.n0
        self.n_traj=np.array([n])
        self.r_traj=np.array([r])
        self.dr=np.zeros(0) #[0]
        self.r_terminating=np.zeros(0) #[]
        #n_terminating=0
        n0=self.n0
        self.startend=np.zeros((n0,2)) #keep the starting and ending r for each dendrite
        self.startend_b=np.zeros((0,2))
        self.startend_t=np.zeros((0,2))
        surviving=np.arange(n0)
        nindexing=n0
        
        KB_total = movingMeanParam(probability = 'KB',neuronType = self.neuronType, labname = self.labname, \
                                 neuronName = self.neuronName, dendrite = self.dendrite) 
        KT_total = movingMeanParam(probability = 'KT',neuronType = self.neuronType, labname = self.labname, \
                                 neuronName = self.neuronName, dendrite = self.dendrite)
        self.rmax = KB_total.shape[0]
        self.tmax = 2 * self.rmax
        
        while n>0 and r < self.tmax:
              
            lower = int(r)
            upper = int(r) + 1
            index_max = self.rmax - 1 
            
            if upper > index_max :
                upper = index_max 
            if lower > index_max:
                lower = index_max
            
            """ for kb """
            kb_up =  KB_total['KB'][upper]
            kb_low =  KB_total['KB'][lower]
            r_up =  KB_total.index[upper]
            r_low =  KB_total.index[lower]        
               
            if np.isnan(kb_up) or np.isnan(kb_low):
                # final kb is constant due to the experiment result 
                select = np.isnan(KB_total['KB']) == False 
                kb = float(KB_total[select]['KB'].tail(1))
                       
            else:
                kb = (kb_up - kb_low) / (r_up - r_low) * (r - r_low) + kb_low
            
            select_0 = KB_total['KB'] != 0
            minKB = min(KB_total[select_0]['KB'].index)
            if r <= minKB and kb == 0:
                kb = 1/3 * float(KB_total[select_0]['KB'].head(1))
                
            """ for kt """
            kt_up = KT_total['KT'][upper]
            kt_low = KT_total['KT'][lower]
            r_up = KT_total.index[upper]
            r_low = KT_total.index[lower]
            
            if np.isnan(kt_up) or np.isnan(kt_low):
                # final kb is constant due to the experiment result 
                select = np.isnan(KT_total['KT']) == False
                kt = float(KT_total[select]['KT'].tail(1))
                if kt == 0:
                    KT_rnan = KT_total[select]['KT']
                    select_0 = KT_rnan != 0.0
                    kb = float(KT_rnan[select_0].tail(1))
                
            else:
                kt = (kt_up - kt_low) / (r_up - r_low) * (r - r_low) + kt_low         
            #print('[kb, kt]  = ', [kb, kt])
                     
            prop = [1*n*kb, 1*n*kt]
            #print(prop)
            a=sum(prop)
            #print(a)
            deltar =-1./a*np.log(rand.random())
            #print(deltar)
         
            r = r + deltar
            self.dr=np.append(self.dr,deltar)#dr.append(deltar)
            #print(r)
            path_length = []
            if rand.random()<prop[0]/a:
                i_branch=rand.randrange(n)
                n = n+1 # branching event
                nindexing += 2
                branch = surviving[i_branch]
                surviving = np.delete(surviving, i_branch)
                self.startend[branch,1] = r
                self.startend_b = np.append(self.startend_b, [self.startend[branch]],axis=0)
                self.startend = np.append(self.startend,[[r,0]],axis=0)
                self.startend = np.append(self.startend,[[r,0]],axis=0)
                surviving = np.append(surviving, [nindexing-2, nindexing-1])
                
            else:
                i_term=rand.randrange(n) #index of the terminating dendrite
                n=n-1 # terminating event
                self.r_terminating=np.append(self.r_terminating,r)#r_terminating.append(r)
                term=surviving[i_term]
                surviving=np.delete(surviving,i_term)
                self.startend[term,1]=r
                self.startend_t=np.append(self.startend_t, [self.startend[term]],axis=0)
                path_length = path_length.append(r)
                #n_terminating+=1
            
            delta = 0.33 
            
            if path_length is []:
                mean_path_length = sum(path_length) / len(path_length)
                if mean_path_length** delta * sum(self.dr) >= 1:
                    print('cell stops to grow')
                    break
                
#            if n!=0:
            self.n_traj=np.append(self.n_traj,n)
            self.r_traj=np.append(self.r_traj,r)
            #print('r and E' + str(r) + '_' + str(E))
        #end while loop here
                 
        return 
class Neuron_two_Tm2_movemean(Neuron_two_Tm1_movemean):
    type='two-stage growing neoron'
    def nrtraj(self):
        super(Neuron_two_Tm2_movemean, self).nrtraj()
        self.neuronName = 'Tm2'

class Neuron_two_Tm9_movemean(Neuron_two_Tm1_movemean):
    type='two-stage growing neoron'
    def nrtraj(self):
        super(Neuron_two_Tm9_movemean, self).nrtraj()
        self.neuronName = 'Tm9'

class Neuron_two_Tm20_movemean(Neuron_two_Tm1_movemean):
    type='two-stage growing neoron'
    def nrtraj(self):
        super(Neuron_two_Tm20_movemean, self).nrtraj()
        self.neuronName = 'Tm20'   

