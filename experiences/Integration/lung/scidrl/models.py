#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from keras.layers import Dense,Input,Dropout,Lambda,Layer,BatchNormalization,Activation
from keras import initializers,regularizers,metrics
import keras.backend as K
from keras.models import Model
from keras.optimizers import Adam
from keras.callbacks import Callback,EarlyStopping
from keras.utils import plot_model
from .flip_gradient import GradientReversal
import tensorflow as tf
import numpy as np
tf.compat.v1.disable_eager_execution()


# In[ ]:



# In[ ]:


class SCIDRL():
    def __init__(self,in_dim,in_category=1,latent_dim=10,class_dim=2,lr=0.001,gamma=1,fg_lambda=1,acts='sigmoid'):
        self.latent=latent_dim
        self.in_dim=in_dim
        self.in_category=in_category
        self.class_dim=class_dim
        self.learning_rate=lr
        self.gamma=gamma
        self.fg_lambda=fg_lambda
        self.acts=acts
    def build(self):
        INPUTS=Input(shape=(self.in_dim+self.in_category+self.latent,))
        expr_in=Lambda(lambda x:x[:,:self.in_dim])(INPUTS)
        batch_in=Lambda(lambda x:x[:,self.in_dim:(self.in_dim+self.in_category)])(INPUTS)
        latent_in=Lambda(lambda x:x[:,(self.in_dim+self.in_category):])(INPUTS)
        h0=Dropout(0.5)(expr_in)
        in_dim=self.in_dim
        gamma=self.gamma
        acts=self.acts
        #####encoder######
        h1=Dense(units=256, activation=None, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros', kernel_regularizer=None,name='encoder3')(h0)
        h1=BatchNormalization()(h1)
        h1=Activation('relu')(h1)
        #####latent#######
        z=Dense(units=self.latent, activation=None, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1),
                bias_initializer='zeros', kernel_regularizer=None,name='latent_mean')(h1)
        z_noise=Lambda(lambda x:x[:,:self.class_dim])(z)
        z_bio=Lambda(lambda x:x[:,self.class_dim:])(z)
        #####adversarial layer######
        flip=GradientReversal(self.fg_lambda)
        ad_h1=flip(z_bio)
        ad_h2=Dense(units=32, activation=None, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros', kernel_regularizer=None,name='ad_hidden2')(ad_h1)
        ad_h2=BatchNormalization()(ad_h2)
        ad_h2=Activation('relu')(ad_h2)
        ad_h3=Dense(units=16, activation=None, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros', kernel_regularizer=None,name='ad_hidden3')(ad_h2)
        ad_h3=BatchNormalization()(ad_h3)
        ad_h3=Activation('relu')(ad_h3)
        ad_output=Dense(units=self.in_category, activation=self.acts, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros', kernel_regularizer=None,name='ad_output')(ad_h3)
        ###classifier layer
        class_output=Dense(units=self.in_category, activation=self.acts, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros',name='classifier_final')(z_noise)
        #####Decoder####
        decoder_h1=Dense(units=256,activation=None, kernel_initializer=initializers.TruncatedNormal(mean=0.0, stddev=0.1), 
                 bias_initializer='zeros',name='decoder_1')
        recon_x=Dense(units=self.in_dim,activation='sigmoid')
        
        dh1=decoder_h1(z)
        dh1=BatchNormalization()(dh1)
        dh1=Activation('relu')(dh1)
        dh_output=recon_x(dh1)
        ####correct output###
        cdh1=decoder_h1(latent_in)
        cdh1=BatchNormalization()(cdh1)
        cdh1=Activation('relu')(cdh1)
        cdh_output=recon_x(cdh1)
        
        ##计算loss层
        class CoustomLossLayerReconstruction(Layer):
            def __init__(self,**kwargs):
                self.is_placeholder = True
                self.in_dim=in_dim
                super(CoustomLossLayerReconstruction, self).__init__(**kwargs)
            def call(self,inputs):
                x=inputs[0]
                x_decoded_mean=inputs[1]
                xent_loss=in_dim*metrics.binary_crossentropy(x, x_decoded_mean)
                self.add_loss(K.mean(xent_loss), inputs=inputs)
                return xent_loss
        class CoustomLossLayerAdversarial(Layer):
            def __init__(self,**kwargs):
                self.is_placeholder = True
                self.acts=acts
                super(CoustomLossLayerAdversarial, self).__init__(**kwargs)
            def call(self, inputs):
                y=inputs[0]
                y_pred=inputs[1]
                if self.acts=='sigmoid':
                    adversarial_loss=metrics.binary_crossentropy(y, y_pred)
                else:
                    adversarial_loss=metrics.categorical_crossentropy(y, y_pred)
                self.add_loss(K.mean(adversarial_loss), inputs=inputs)
                return adversarial_loss       
        class CoustomLossLayerClassifier(Layer):
            def __init__(self,**kwargs):
                self.is_placeholder = True
                self.gamma=gamma
                self.acts=acts
                super(CoustomLossLayerClassifier, self).__init__(**kwargs)
            def call(self, inputs):
                y=inputs[0]
                y_pred=inputs[1]
                if self.acts=='sigmoid':
                    classifier_loss=metrics.binary_crossentropy(y, y_pred)
                else:
                    classifier_loss=metrics.categorical_crossentropy(y, y_pred)
                self.add_loss(self.gamma*K.mean(classifier_loss), inputs=inputs)
                return classifier_loss
    
        
        y_ent=CoustomLossLayerReconstruction(name = "myCustomLossLayerReconstruction")([expr_in,dh_output])
        y_classifier=CoustomLossLayerClassifier(name="myCoustomLossLayerClassifier")([batch_in,class_output])
        y_adversarial=CoustomLossLayerAdversarial(name="myCoustomLossLayerAdversarial")([batch_in,ad_output])
        my_model=Model(inputs=INPUTS, outputs=[y_ent,y_classifier,y_adversarial])
        my_result=Model(inputs=INPUTS,outputs=[z,cdh_output])
        op=Adam(lr=self.learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
        my_model.compile(optimizer=op, loss=None)
        self.my_model=my_model
        self.my_result=my_result

