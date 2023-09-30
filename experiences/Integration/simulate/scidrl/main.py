#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import random
import matplotlib.pyplot as plt
import umap.umap_ as umap
import time
from sklearn.preprocessing import minmax_scale
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from pandas import Series,DataFrame
import pandas as pd
import os
import tensorflow as tf
import keras
from copy import deepcopy
from .models import *
import argparse


# In[ ]:


class SCIDRL_train:
    def __init__(self,params,data_file,meta_file):
        super(SCIDRL_train, self).__init__()
        self.params=params
        self.data_file=data_file
        self.meta_file=meta_file
        self.data1,self.batch1=self.read_data()
        self.inputs=np.hstack((self.data1,self.batch1,np.zeros((self.data1.shape[0],self.params.zdim))))
        self.model=SCIDRL(self.data1.shape[1],self.batch1.shape[1],self.params.zdim,self.params.znoise_dim,
                          self.params.lr,self.params.gamma,self.params.fg_lambda,self.params.acts)
    def train(self):
        inputs=self.inputs
        self.model.build()
        class History(Callback):
            def on_train_begin(self, logs={}):
                self.history = {}
                self.epoch = [] 
            def on_train_end(self, logs={}):
                return 
            def on_epoch_begin(self, epoch, logs={}):
                return 
            def on_epoch_end(self, epoch, logs={}):
                self.epoch.append(epoch)
                # Get the output layer names (we only have loss layers as output)
                output_layer_names = [output_layer.name for output_layer in self.model.output]
                # Get the mean error over all batches in this epoch
                output_layer_values = np.mean(self.model.predict(inputs), axis = 1)
                # Store it to the history
                for k, v in zip(output_layer_names, output_layer_values):
                    self.history.setdefault(k, []).append(v)
                return
            def on_batch_begin(self, batch, logs={}):
                return
            def on_batch_end(self, batch, logs={}):
                return
        myhistory=History()
        loss_=self.model.my_model.fit(x=inputs,epochs=self.params.epochs,batch_size=self.params.batch_size,shuffle=True,verbose=2,callbacks=[myhistory])
        total_loss=loss_.history['loss']
        myhistory.history['total_loss']=total_loss
        return myhistory.history
    
    def infer(self):
        embed=self.model.my_result.predict(self.inputs)[0]
        a=np.hstack((np.zeros((self.data1.shape[0],self.params.znoise_dim)),embed[:,self.params.znoise_dim:]))
        inputs_test_1=np.hstack((self.data1,self.batch1,a))
        correct=self.model.my_result.predict(inputs_test_1)[1]
        return embed,correct
#     def model_save(self,model_file):
#         self.model.save_weights(model_file)
#     def model_load(self,model_file):
#         sefl.model.load_weights(model_file)
    def read_data(self):
        self.data=pd.read_csv(self.data_file,header=0,index_col=0)
        self.meta=pd.read_csv(self.meta_file,header=0,index_col=0)
        cells=np.intersect1d(self.data.index,self.meta.index)
        self.meta=self.meta.loc[cells,:]
        self.data=self.data.loc[cells,:]
        if self.params.minmaxscale:
            data1=minmax_scale(self.data.values,axis=1)
        else:
            data1=self.data.values
        batch=self.meta.loc[:,'batch'].values
        ub=np.unique(batch)
        if len(ub)==2:
            batch1=np.zeros((batch.shape[0],1))
            batch1[np.where(batch==ub[1])[0],:]=1
        else:
            batch1=np.zeros((len(batch),len(ub)))
            for i in np.arange(len(ub)):
                batch1[np.where(batch==ub[i])[0],i]=1
        return data1,batch1

