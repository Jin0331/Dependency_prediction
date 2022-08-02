import os
from os.path import exists
import pickle
from tqdm import tqdm
import pandas as pd
import datatable as dt
import datetime
import numpy as np
import time
import pandas as pd
import datatable as dt
import matplotlib.pyplot as plt
import gc
from pathlib import Path

import tensorflow as tf
import pickle
from keras.callbacks import ModelCheckpoint
from keras import models
from keras import backend as K
from keras.layers import Dense, Merge
from keras.callbacks import EarlyStopping
from sklearn.model_selection import cross_val_score, KFold


### global variable
activation_func = 'relu' # for all middle layers
activation_func2 = 'linear' # for output layer to output unbounded gene-effect scores
init = 'he_uniform'
dense_layer_dim = 250
batch_size = 200


def load_data(filename):
    df_dt = dt.fread(filename, sep="\t")
    gene_names = df_dt[0].to_list()[0]
    sample_names = list(df_dt.names[1:])
    data_np = df_dt[:, 1:].to_numpy()
    data_np = np.transpose(data_np)
    gc.collect()
    
    return data_np, sample_names, gene_names

def load_data_prediction(filename):
    data = []
    gene_names = []
    data_labels = []
    lines = open(filename).readlines()
    sample_names = lines[0].replace('\n', '').split('\t')[1:]
    dx = 1

    for line in lines[dx:]:
        values = line.replace('\n', '').split('\t')
        gene = str.upper(values[0])
        gene_names.append(gene)
        data.append(values[1:])
    data = np.array(data, dtype='float32')
    data = np.transpose(data)

    return data, data_labels, sample_names, gene_names


def trainset_load(npz_path="prediction/data/ccl_complete_data_501CCL_1298DepOI_614727samples.npz", 
                  filename="training_custom", mut=True, exp=True, cna=True, meth=True):
    
    if exists(npz_path) is False:
        # expression
        if exp:
            if exists(TEMP_PATH + filename + "_exp_training.npy") is False:
                data_exp, sample_names_exp, property_names_exp = load_data(TRAIN_PATH + filename + "_exp_training.txt")
                np.save(TEMP_PATH + filename + "_exp_training.npy", data_exp)
            else:
                data_exp = np.load(TEMP_PATH + filename + "_exp_training.npy")
            print("exp load.. completed")    

        # mutation
        if mut:
            if exists(TEMP_PATH + filename + "_mut_training.npy") is False:
                data_mut, sample_names_mut, property_names_mut = load_data(TRAIN_PATH + filename + "_mut_training.txt")
                np.save(TEMP_PATH + filename + "_mut_training.npy", data_mut)
            else:
                data_mut = np.load(TEMP_PATH + filename + "_mut_training.npy")
            print("mut load.. completed")

        # copy number alteration
        if cna:
            if exists(TEMP_PATH + filename + "_cna_training.npy") is False:
                data_cna, sample_names_cna, property_names_cna = load_data(TRAIN_PATH + filename + "_cna_training.txt")
                np.save(TEMP_PATH + filename + "_cna_training.npy", data_cna)
            else:
                data_cna = np.load(TEMP_PATH + filename + "_cna_training.npy")
            print("cna load.. completed")

        # methylation
        if meth:
            if exists(TEMP_PATH + filename + "_meth_training.npy") is False:
                data_meth, sample_names_meth, property_names_meth = load_data(TRAIN_PATH + filename + "_meth_training.txt")
                np.save(TEMP_PATH + filename + "_meth_training.npy", data_meth)
            else:
                data_meth = np.load(TEMP_PATH + filename + "_meth_training.npy")
            print("meth load.. completed")

        # dependency score
        if exists(TEMP_PATH + filename + "_DepScore_training.npy") is False:
            data_dep, sample_names_dep, property_names_dep = load_data(TRAIN_PATH + filename + "_DepScore_training.txt")
            np.save(TEMP_PATH + filename + "_DepScore_training.npy", data_dep)
        else:
            data_dep = np.load(TEMP_PATH + filename + "_DepScore_training.npy")
        print("dep load.. completed")

        # fingerprint
        if exists(TEMP_PATH + filename + "_fingerprint_training.npy") is False:
            data_fprint, sample_names_fprint, property_names_fprint = load_data(TRAIN_PATH + filename + "_fingerprint_training.txt")
            np.save(TEMP_PATH + filename + "_fingerprint_training.npy", data_fprint)
        else:
            data_fprint = np.load(TEMP_PATH + filename + "_fingerprint_training.npy")
        print("fingerprint load.. completed")

        np.savez_compressed(npz_path, 
                            data_exp=data_exp, data_mut=data_mut, data_cna=data_cna, data_meth=data_meth,
                            data_dep=data_dep, data_fprint=data_fprint)
    else :
        data = np.load(npz_path)
        data_exp = data['data_exp']
        data_mut = data['data_mut']
        data_cna = data['data_cna']
        data_meth = data['data_meth']
        data_fprint = data['data_fprint']
        data_dep = data['data_dep']
    
    # gc load
    gc.collect()
    
    return data_exp, data_mut, data_cna, data_meth, data_fprint, data_dep


def full_model(data_mut, data_exp, data_cna, data_meth,
               data_fprint, data_dep, id_train, id_test, 
               premodel_mut, premodel_exp, premodel_cna, premodel_meth,
               save_path):
    t = time.time()
    filepath=save_path + "full_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_mut = models.Sequential()
            model_mut.add(Dense(output_dim=1000, input_dim=premodel_mut[0][0].shape[0], activation=activation_func,
                                weights=premodel_mut[0], trainable=True))
            model_mut.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, weights=premodel_mut[1],
                                trainable=True))
            model_mut.add(Dense(output_dim=50, input_dim=100, activation=activation_func, weights=premodel_mut[2],
                                trainable=True))

            # subnetwork of expression
            model_exp = models.Sequential()
            model_exp.add(Dense(output_dim=500, input_dim=premodel_exp[0][0].shape[0], activation=activation_func,
                                weights=premodel_exp[0], trainable=True))
            model_exp.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_exp[1],
                                trainable=True))
            model_exp.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_exp[2],
                                trainable=True))

            # subnetwork of copy number alterations
            model_cna = models.Sequential()
            model_cna.add(Dense(output_dim=500, input_dim=premodel_cna[0][0].shape[0], activation=activation_func,
                                weights=premodel_cna[0], trainable=True))
            model_cna.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_cna[1],
                                trainable=True))
            model_cna.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_cna[2],
                                trainable=True))

            # subnetwork of DNA methylations
            model_meth = models.Sequential()
            model_meth.add(Dense(output_dim=500, input_dim=premodel_meth[0][0].shape[0], activation=activation_func,
                                 weights=premodel_meth[0], trainable=True))
            model_meth.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_meth[1],
                                 trainable=True))
            model_meth.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_meth[2],
                                 trainable=True))

            # subnetwork of gene fingerprints
            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            # prediction network
            model_final = models.Sequential()
            model_final.add(Merge([model_mut, model_exp, model_cna, model_meth, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=250, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init,
                                  trainable=True))
        else :
            model_final = models.load_model(filepath)
            
        # callback list
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=1, mode='min') # early stopping       
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')

        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_mut[id_train], data_exp[id_train], data_cna[id_train], data_meth[id_train],
                         data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_mut[id_test], data_exp[id_test], data_cna[id_test], data_meth[id_test],
                         data_fprint[id_test]], data_dep[id_test], verbose=0, batch_size=batch_size)
        print("\n\nFull-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs   

def mut_exp_cna_model(data_mut, data_exp, data_cna, 
               data_fprint, data_dep, id_train, id_test, 
               premodel_mut, premodel_exp, premodel_cna,
               save_path):
    t = time.time()
    filepath=save_path + "mut_exp_cna_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:        
            model_mut = models.Sequential()
            model_mut.add(Dense(output_dim=1000, input_dim=premodel_mut[0][0].shape[0], activation=activation_func,
                                weights=premodel_mut[0], trainable=True))
            model_mut.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, weights=premodel_mut[1],
                                trainable=True))
            model_mut.add(Dense(output_dim=50, input_dim=100, activation=activation_func, weights=premodel_mut[2],
                                trainable=True))

            # subnetwork of expression
            model_exp = models.Sequential()
            model_exp.add(Dense(output_dim=500, input_dim=premodel_exp[0][0].shape[0], activation=activation_func,
                                weights=premodel_exp[0], trainable=True))
            model_exp.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_exp[1],
                                trainable=True))
            model_exp.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_exp[2],
                                trainable=True))

            # subnetwork of copy number alterations
            model_cna = models.Sequential()
            model_cna.add(Dense(output_dim=500, input_dim=premodel_cna[0][0].shape[0], activation=activation_func,
                                weights=premodel_cna[0], trainable=True))
            model_cna.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_cna[1],
                                trainable=True))
            model_cna.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_cna[2],
                                trainable=True))

            # subnetwork of gene fingerprints
            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            # prediction network
            model_final = models.Sequential()
            model_final.add(Merge([model_mut, model_exp, model_cna, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=200, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init,
                                  trainable=True))
        else :
            model_final = models.load_model(filepath)
            
        # callback
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=1, mode='min') # early stopping
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')

        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_mut[id_train], data_exp[id_train], data_cna[id_train],
                         data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_mut[id_test], data_exp[id_test], data_cna[id_test],
                         data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nMut_Exp_CNA-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs   

def mut_exp_meth_model(data_mut, data_exp, data_meth, 
               data_fprint, data_dep, id_train, id_test, 
               premodel_mut, premodel_exp, premodel_meth,
               save_path):
    t = time.time()
    filepath=save_path + "mut_exp_meth_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_mut = models.Sequential()
            model_mut.add(Dense(output_dim=1000, input_dim=premodel_mut[0][0].shape[0], activation=activation_func,
                                weights=premodel_mut[0], trainable=True))
            model_mut.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, weights=premodel_mut[1],
                                trainable=True))
            model_mut.add(Dense(output_dim=50, input_dim=100, activation=activation_func, weights=premodel_mut[2],
                                trainable=True))

            # subnetwork of expression
            model_exp = models.Sequential()
            model_exp.add(Dense(output_dim=500, input_dim=premodel_exp[0][0].shape[0], activation=activation_func,
                                weights=premodel_exp[0], trainable=True))
            model_exp.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_exp[1],
                                trainable=True))
            model_exp.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_exp[2],
                                trainable=True))

            # subnetwork of DNA methylations
            model_meth = models.Sequential()
            model_meth.add(Dense(output_dim=500, input_dim=premodel_meth[0][0].shape[0], activation=activation_func,
                                 weights=premodel_meth[0], trainable=True))
            model_meth.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_meth[1],
                                 trainable=True))
            model_meth.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_meth[2],
                                 trainable=True))

            # subnetwork of gene fingerprints
            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            # prediction network
            model_final = models.Sequential()
            model_final.add(Merge([model_mut, model_exp, model_meth, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=200, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init,
                                  trainable=True))
        else :
            model_final = models.load_model(filepath)
            
        # callback list
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=1, mode='min') # early stopping
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')

        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_mut[id_train], data_exp[id_train], data_meth[id_train],
                         data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_mut[id_test], data_exp[id_test], data_meth[id_test],
                         data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nMut_Exp_Meth-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs   
    

    
def mut_exp_model(data_mut, data_exp, data_fprint, data_dep, id_train, id_test, 
                  premodel_mut, premodel_exp, save_path):
    t = time.time()
    filepath=save_path + "mut_exp_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_mut = models.Sequential()
            model_mut.add(Dense(output_dim=1000, input_dim=premodel_mut[0][0].shape[0], activation=activation_func,
                                weights=premodel_mut[0], trainable=True))
            model_mut.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, weights=premodel_mut[1],
                                trainable=True))
            model_mut.add(Dense(output_dim=50, input_dim=100, activation=activation_func, weights=premodel_mut[2],
                                trainable=True))

            # subnetwork of expression
            model_exp = models.Sequential()
            model_exp.add(Dense(output_dim=500, input_dim=premodel_exp[0][0].shape[0], activation=activation_func,
                                weights=premodel_exp[0], trainable=True))
            model_exp.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_exp[1],
                                trainable=True))
            model_exp.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_exp[2],
                                trainable=True))

            # subnetwork of gene fingerprints
            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            # prediction network
            model_final = models.Sequential()
            model_final.add(Merge([model_mut, model_exp, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=150, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init,
                                  trainable=True))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init,
                                  trainable=True))
        else :
            model_final = models.load_model(filepath)
            
        # callback
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=1, mode='min') # early stopping
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')

        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_mut[id_train], data_exp[id_train], data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_mut[id_test], data_exp[id_test], data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nMut_Exp-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs   


def exp_model(data_exp, data_fprint, data_dep, id_train, id_test, premodel_exp, save_path):
    t = time.time()
    filepath=save_path + "exp_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_exp = models.Sequential()
            model_exp.add(Dense(output_dim=500, input_dim=premodel_exp[0][0].shape[0], activation=activation_func,
                                weights=premodel_exp[0], trainable=True))
            model_exp.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_exp[1],
                                trainable=True))
            model_exp.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_exp[2],
                                trainable=True))

            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            model_final = models.Sequential()
            model_final.add(Merge([model_exp, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=100, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init))
        else :
            model_final = models.load_model(filepath)

        # callback-list
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=0, mode='min')
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')

        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_exp[id_train], data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_exp[id_test], data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nExp-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs     
    
def mut_model(data_mut, data_fprint, data_dep, id_train, id_test, premodel_mut, save_path):
    t = time.time()
    filepath=save_path + "mut_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_mut = models.Sequential()
            model_mut.add(Dense(output_dim=1000, input_dim=premodel_mut[0][0].shape[0], activation=activation_func,
                                weights=premodel_mut[0], trainable=True))
            model_mut.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, weights=premodel_mut[1],
                                trainable=True))
            model_mut.add(Dense(output_dim=50, input_dim=100, activation=activation_func, weights=premodel_mut[2],
                                trainable=True))

            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            model_final = models.Sequential()
            model_final.add(Merge([model_mut, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=100, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init))
        else :
            model_final = models.load_model(filepath)
            
        # callback list
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=0, mode='min')
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
        
        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_mut[id_train], data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_mut[id_test], data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nMut-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs

def meth_model(data_meth, data_fprint, data_dep, id_train, id_test, premodel_meth, save_path):
    t = time.time()
    filepath=save_path + "meth_weights.best.hdf5"
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_meth = models.Sequential()
            model_meth.add(Dense(output_dim=500, input_dim=premodel_meth[0][0].shape[0], activation=activation_func,
                                 weights=premodel_meth[0], trainable=True))
            model_meth.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_meth[1],
                                 trainable=True))
            model_meth.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_meth[2],
                                 trainable=True))

            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            model_final = models.Sequential()
            model_final.add(Merge([model_meth, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=100, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init))
        else :
            model_final = models.load_model(filepath)
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=0, mode='min')
        
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_meth[id_train], data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_meth[id_test], data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nMeth-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs

    
def cna_model(data_cna, data_fprint, data_dep, id_train, id_test, premodel_cna, save_path):
    t = time.time()
    filepath=save_path + "cna_weights.best.hdf5"
    
    with tf.device('/cpu:0'):
        if exists(filepath) is False:
            model_cna = models.Sequential()
            model_cna.add(Dense(output_dim=500, input_dim=premodel_cna[0][0].shape[0], activation=activation_func,
                                weights=premodel_cna[0], trainable=True))
            model_cna.add(Dense(output_dim=200, input_dim=500, activation=activation_func, weights=premodel_cna[1],
                                trainable=True))
            model_cna.add(Dense(output_dim=50, input_dim=200, activation=activation_func, weights=premodel_cna[2],
                                trainable=True))

            model_gene = models.Sequential()
            model_gene.add(Dense(output_dim=1000, input_dim=data_fprint.shape[1], activation=activation_func, init=init,
                                 trainable=True))
            model_gene.add(Dense(output_dim=100, input_dim=1000, activation=activation_func, init=init, trainable=True))
            model_gene.add(Dense(output_dim=50, input_dim=100, activation=activation_func, init=init, trainable=True))

            model_final = models.Sequential()
            model_final.add(Merge([model_cna, model_gene], mode='concat'))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=100, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=dense_layer_dim, input_dim=dense_layer_dim, activation=activation_func, init=init))
            model_final.add(Dense(output_dim=1, input_dim=dense_layer_dim, activation=activation_func2, init=init))
        else :
            model_final = models.load_model(filepath)
            
        # callback-list
        history = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=0, mode='min')
        checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
        
        model_final.compile(loss='mse', optimizer='adam')
        hs = model_final.fit([data_cna[id_train], data_fprint[id_train]], data_dep[id_train], nb_epoch=30,
                    validation_split=2/9, batch_size=batch_size, shuffle=True, callbacks=[history, checkpoint])
        cost_testing = model_final.evaluate([data_cna[id_test], data_fprint[id_test]], data_dep[id_test], verbose=0,
                                        batch_size=batch_size)
        print("\n\nCNA-DeepDEP model training completed in %.1f mins.\nloss:%.4f valloss:%.4f testloss:%.4f" % ((time.time() - t)/60, history.model.model.history.history['loss'][history.stopped_epoch], history.model.model.history.history['val_loss'][history.stopped_epoch], cost_testing))
        
        return model_final, hs



class generator(tf.keras.utils.Sequence):
    def __init__(self, x_mut, x_exp, x_cna, x_meth, x_fprint, y_dep, batch_size):
        self.x_mut, self.x_exp, self.x_cna, self.x_meth, self.x_fprint, self.y_dep = x_mut, x_exp, x_cna, x_meth, x_fprint, y_dep
        self.bs = batch_size
    
    def __len__(self):
        return (len(self.y_dep) - 1) // self.bs + 1
    
    def __getitem__(self,idx):
        start, end = idx * self.bs, (idx+1) * self.bs
        return (self.x_mut[start:end], self.x_exp[start:end], self.x_cna[start:end], self.x_meth[start:end], self.x_fprint[start:end]), self.y_dep[start:end]

class generator3(generator):
    def __init__(self, x_mut, x_exp, x_cna, x_fprint, y_dep, batch_size):
        self.x_mut, self.x_exp, self.x_cna, self.x_fprint, self.y_dep = x_mut, x_exp, x_cna, x_fprint, y_dep
        self.bs = batch_size
        
    def __getitem__(self, idx):
        start, end = idx * self.bs, (idx+1) * self.bs
        return (self.x_mut[start:end], self.x_exp[start:end], self.x_cna[start:end], self.x_fprint[start:end]), self.y_dep[start:end]
    
class generator2(generator):
    def __init__(self, x_mut, x_exp, x_fprint, y_dep, batch_size):
        self.x_mut, self.x_exp, self.x_fprint, self.y_dep = x_mut, x_exp, x_fprint, y_dep
        self.bs = batch_size
        
    def __getitem__(self, idx):
        start, end = idx * self.bs, (idx+1) * self.bs
        return (self.x_mut[start:end], self.x_exp[start:end], self.x_fprint[start:end]), self.y_dep[start:end]

class generator1(generator):
    def __init__(self, x_mut, x_fprint, y_dep, batch_size):
        self.x_mut, self.x_fprint, self.y_dep = x_mut, x_fprint, y_dep
        self.bs = batch_size
        
    def __getitem__(self, idx):
        start, end = idx * self.bs, (idx+1) * self.bs
        return (self.x_mut[start:end], self.x_fprint[start:end]), self.y_dep[start:end]
    
def root_mean_squared_error(y_true, y_pred):
        return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1)) 

def rmse(y_true, y_pred):
    """
    Root Mean Squared Error
    Args:
        y_true ([np.array]): test samples
        y_pred ([np.array]): predicted samples
    Returns:
        [float]: root mean squared error
    """
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

    
def adj_r2(y_true, y_pred):
    """
    Adjusted R2 regression score function with default inputs.
    Best possible score is 1.0, lower values are worse.
    Args:
        y_true ([np.array]): test samples
        y_pred ([np.array]): predicted samples
    Returns:
        [float]: adjusted R2
    """
    SS_res =  tf.reduce_sum(tf.square(y_true - y_pred), axis=-1)
    SS_tot = tf.reduce_sum(tf.square(y_true - tf.reduce_mean(y_true, axis=-1)), axis=-1)
    return (1 - SS_res/(SS_tot + tf.keras.backend.epsilon())) * (1 - (1 - r2(y_true, y_pred)) * (tf.cast(tf.size(y_true), tf.float32) - 1) / (tf.cast(tf.size(y_true), tf.float32) - tf.cast(tf.rank(y_true), tf.float32) - 1))    

def coeff_determination(y_true, y_pred):
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

def model_train_vis(history, save_path, filename):
    loss = history.history['loss']
    val_loss = history.history['val_loss']

    plt.figure()
    plt.plot(loss, 'bo', label='Training')
    plt.plot(val_loss, 'b', label='Validation')
    plt.title('Training and validation loss(MSE)')
    plt.legend()
    plt.savefig(save_path + 'model_vis-' + filename + '.png', dpi=300)
