from __future__ import division
import numpy as np
import h5py
import scipy.io as sio
import hdf5storage
import random
import matplotlib.pyplot as plt
import cv2
from skimage import color,transform
from scipy import ndimage
from PIL import Image 

def Data_Division(list_truth,list_meas,mask,sample_num):
  
    (train_num,valid_num,test_num)  = sample_num
    pair_train = (list_meas[:train_num],list_truth[:train_num])
    pair_valid = (list_meas[train_num:train_num+valid_num],list_truth[train_num:train_num+valid_num])                     
    pair_test =  (list_meas[-test_num:],list_truth[-test_num:])   
    

    return pair_train, pair_test, pair_valid, mask

def Data_Generator_1U1(dataset, mask, batch_size, is_training=True):
    """
    :param dataset: the raw data, containing groups of frames for train/test/valid 
    :param cube_size: tuple, contrain the hyperparameter of [hcube, hstride, wcube, wstride]
    :param batch_size: 
    :return: 
    """
    (measure, ground) = dataset
    num_sample = len(measure)   
    (height, width) = measure[0].shape
    nF = ground[0].shape[2]
    index = np.random.choice(num_sample, size=num_sample, replace=False).astype(np.int16)
    sample_cnt,batch_cnt,list_measure,list_ground,list_index = 0,0,[],[],[]
    H, W= 512, 512
    
    while True:
        if (sample_cnt < num_sample):
            if is_training is True:
                ind_set = index[sample_cnt]
            else:
                ind_set = sample_cnt
            meas = measure[ind_set]
            ### shot noise injection
            #QE, bit = 0.4, 512
            #meas_max = meas.max()
            #meas = meas/meas_max
            #meas = np.random.binomial((meas*bit/QE).astype(int),QE)
            #meas = meas/bit
            #meas = meas*meas_max
            
            #  initialization
            C = np.sum(mask**2,2)
            #C[C==0]=1
            meas_temp = meas/C
            meas_temp = np.tile(meas_temp[:,:,np.newaxis],(1,1,nF))
            meas = meas_temp*mask
            list_measure.append(meas)
            list_ground.append(ground[ind_set])
            batch_cnt += 1
            sample_cnt += 1
            
            if batch_cnt == batch_size:
                batch_measure,batch_ground = np.stack(list_measure,0),np.stack(list_ground,0)
                height_init,batch_cnt,list_measure,list_ground = 0,0,[],[]
                #print batch_measure.shape,batch_ground.shape
                yield batch_measure,mask,batch_ground
        else:            
            sample_cnt = 0
            index = np.random.choice(num_sample, size=num_sample, replace=False).astype(np.int16)


def Data_generation(dataset_name, H, W, nF, is_training=True):
    (data_name,mask_name) = dataset_name                        
    mask_file = sio.loadmat(mask_name+'.mat')
    mask = mask_file['mask']
     #sence cave
    mask=mask[:,:,:nF]
    truth_all,meas_all = [],[]
    train_num,valid_num = 0,0
    
    if is_training==True:
        # training data generation
        # validation data generation
        print('Please generate training data.')
    else:
        # real test data generation
        path = data_name[2]+'_cr_%s.mat'%str(nF)
        meas = sio.loadmat(path)['meas_all']
        test_num = meas.shape[2]
        img2 = np.zeros((H,W,nF))
        for i in range(test_num):
            meas_temp = meas[:,:,i]/255           # meas[:,:,i].max()
            meas_temp = meas_temp*nF/1.5
            truth_all.append(img2)
            meas_all.append(meas_temp)
    sample_num = (train_num,valid_num,test_num)   
    return truth_all, meas_all, mask, sample_num
    
    



         
