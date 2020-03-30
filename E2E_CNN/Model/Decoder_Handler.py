import tensorflow as tf
import numpy as np
import time
import sys
import math
from scipy import ndimage
import matplotlib.pyplot as plt

from E2E_CNN.Lib.Data_Processing import *
from E2E_CNN.Lib.Utility import *
from E2E_CNN.Model.E2E_CNN_model import Depth_Decoder
from E2E_CNN.Model.Base_Handler import Basement_Handler


class Decoder_Handler(Basement_Handler):
    def __init__(self, dataset_name, model_config, Cr, sess, is_training=True):
        
        # Initialization of Configuration, Parameter and Datasets
        super(Decoder_Handler, self).__init__(sess=sess, model_config=model_config, is_training=is_training)
        self.nF = Cr
        self.initial_parameter()
        self.data_assignment(dataset_name)

        # Data Generator
        self.gen_train = Data_Generator_1U1(self.set_train,self.sense_mask,self.batch_size,is_training=True)
        self.gen_valid = Data_Generator_1U1(self.set_valid,self.sense_mask,self.batch_size,is_training=False)
        self.gen_test = Data_Generator_1U1(self.set_test,self.sense_mask,self.batch_size,is_training=False)
        
        # Define the general model and the corresponding input
        shape_meas = (self.batch_size, self.sense_mask.shape[0], self.sense_mask.shape[1],self.nF)
        shape_sense = self.sense_mask.shape
        shape_truth = (self.batch_size,) + self.sense_mask.shape
        self.meas_sample = tf.placeholder(tf.float32, shape=shape_meas, name='input_meas')
        self.sense_matrix = tf.placeholder(tf.float32, shape=shape_sense, name='input_mat')
        self.truth_seg = tf.placeholder(tf.float32, shape=shape_truth, name='output_truth')
        
        # Initialization for the model training procedure.
        self.learning_rate = tf.get_variable('learning_rate', shape=(), initializer=tf.constant_initializer(self.lr_init),trainable=False)
        self.lr_new = tf.placeholder(tf.float32, shape=(), name='lr_new')
        self.lr_update = tf.assign(self.learning_rate, self.lr_new, name='lr_update')
        self.train_test_valid_assignment()
        self.trainable_parameter_info()
        self.saver = tf.train.Saver(tf.global_variables())

    def initial_parameter(self):
        # Configuration Set
        config = self.model_config
        
        # Model Input Initialization
        self.batch_size = int(config.get('batch_size',1))
        self.upbound = float(config.get('upbound',1))
        
        # Initialization for Training Controler
        self.epochs = int(config.get('epochs',100))
        self.patience = int(config.get('patience',30))
        self.lr_init = float(config.get('learning_rate',0.001))
        self.lr_decay_coe = float(config.get('lr_decay',0.1))
        self.lr_decay_epoch = int(config.get('lr_decay_epoch',20))
        self.lr_decay_interval = int(config.get('lr_decay_interval',10))

    def data_assignment(self,dataset_name):
        # Division for train, test and validation
        model_config = self.model_config
        truth_all,meas_all,mask,valid_num = Data_generation(dataset_name, 512, 512, self.nF, is_training=False)        #True
        self.set_train, self.set_test, self.set_valid, self.sense_mask = Data_Division(truth_all,meas_all,mask,valid_num)
        
        disp_train = len(self.set_train[0])
        disp_valid = len(self.set_valid[0])
        disp_test = len(self.set_test[0])
        
        self.train_size = int(np.ceil(float(disp_train)/self.batch_size))
        self.test_size  = int(np.ceil(float(disp_test)/self.batch_size))
        self.valid_size = int(np.ceil(float(disp_valid)/self.batch_size))
        
        # Display the data structure of Training/Testing/Validation Dataset
        print('Available samples (batch) train %d(%d), valid %d(%d), test %d(%d)' % (
            disp_train,self.train_size,disp_valid,self.valid_size,disp_test,self.test_size))
        
    def train_test_valid_assignment(self):#, is_training = True, reuse = False
        
        value_set = (self.meas_sample,
                     tf.expand_dims(self.sense_matrix,0),
                     self.truth_seg)
        
        with tf.name_scope('Train'):
            with tf.variable_scope('Depth_Decoder', reuse=False):
                self.Decoder_train = Depth_Decoder(value_set,self.learning_rate,self.sess,self.model_config,is_training=True)
        with tf.name_scope('Val'):
            with tf.variable_scope('Depth_Decoder', reuse=True):
                self.Decoder_valid = Depth_Decoder(value_set,self.learning_rate,self.sess,self.model_config,is_training=False)
        with tf.name_scope('Test'):
            with tf.variable_scope('Depth_Decoder', reuse=True):
                self.Decoder_test = Depth_Decoder(value_set,self.learning_rate,self.sess,self.model_config,is_training=False)
                
                
    def train(self):
        self.sess.run(tf.global_variables_initializer())
        print ('Training Started')
        if self.model_config.get('model_filename',None) is not None:
            self.restore()
            print('Pretrained Model Downloaded')
        else:
            print('New Model Training')
        epoch_cnt,wait,min_val_loss = 0,0,float('inf')
        
        while epoch_cnt <= self.epochs:
            
            # Training Preparation: Learning rate pre=setting, Model Interface summary.
            start_time = time.time()
            cur_lr = self.calculate_scheduled_lr(epoch_cnt)
            train_fetches = {'global_step': tf.train.get_or_create_global_step(), 
                             'train_op':self.Decoder_train.train_op,
                             'metrics':self.Decoder_train.metrics,
                             'pred_orig':self.Decoder_train.decoded_image,
                             'loss':self.Decoder_train.loss}
            valid_fetches = {'global_step': tf.train.get_or_create_global_step(),
                            'pred_orig':self.Decoder_valid.decoded_image,
                             'metrics':self.Decoder_valid.metrics,
                            'loss':self.Decoder_valid.loss}
            Tresults,Vresults = {"loss":[],"psnr":[],"ssim":[],"mse":[]},{"loss":[],"psnr":[],"ssim":[],"mse":[]}
            
            # Framework and Visualization SetUp for Training 
            #list_truth,list_pred,list_meas = [],[],[]
            for trained_batch in range(0,self.train_size):
                (measure_train,mask_train,ground_train) = self.gen_train.__next__()
                feed_dict_train = {self.meas_sample: measure_train, 
                                   self.sense_matrix: mask_train,
                                   self.truth_seg: ground_train}
                train_output = self.sess.run(train_fetches,feed_dict=feed_dict_train)
                Tresults["loss"].append(train_output['loss'])
                Tresults["psnr"].append(train_output['metrics'][0])
                Tresults["ssim"].append(train_output['metrics'][1])
                Tresults["mse"].append(train_output['metrics'][2])
                message = "Train Epoch [%2d/%2d] Batch [%d/%d] lr: %.4f, loss: %.8f psnr: %.4f" % (
                    epoch_cnt, self.epochs, trained_batch, self.train_size, cur_lr, Tresults["loss"][-1], Tresults["psnr"][-1])
                if trained_batch%10 == 0:
                    print(message)
                
                
            # Framework and Visualization SetUp for Validation 
            validation_time = []
            for valided_batch in range(0,self.valid_size):
                (measure_valid,mask_valid,ground_valid) = self.gen_valid.__next__()
                feed_dict_valid = {self.meas_sample: measure_valid,
                                   self.sense_matrix: mask_valid,
                                   self.truth_seg: ground_valid}
                start_time = time.time()
                valid_output = self.sess.run(valid_fetches,feed_dict=feed_dict_valid)
                end_time = time.time()
                validation_time.append(end_time-start_time)
                valid_output = self.sess.run(valid_fetches,feed_dict=feed_dict_valid)
                Vresults["loss"].append(valid_output['loss'])
                Vresults["psnr"].append(valid_output['metrics'][0])
                Vresults["ssim"].append(valid_output['metrics'][1])
                Vresults["mse"].append(valid_output['metrics'][2])
                message = "Valid Epoch [%2d/%2d] Batch [%d/%d] lr: %.4f, loss: %.8f psnr: %.4f" % (
                    epoch_cnt, self.epochs, valided_batch, self.valid_size, cur_lr, Vresults["loss"][-1], Vresults["psnr"][-1])

            
            # Information Logging for Model Training and Validation (Maybe for Curve Plotting)
            Tloss,Vloss = np.mean(Tresults["loss"]),np.mean(Vresults["loss"])
            train_psnr,valid_psnr = np.mean(Tresults["psnr"]),np.mean(Vresults["psnr"])
            train_ssim,valid_ssim = np.mean(Tresults["ssim"]),np.mean(Vresults["ssim"])
            train_mse, valid_mse  = np.mean(Tresults["mse"]), np.mean(Vresults["mse"])
            summary_format = ['loss/train_loss','loss/valid_loss','metric/train_psnr','metric/train_ssim',
                              'metric/valid_psnr','metric/valid_ssim']
            summary_data = [Tloss, Vloss, train_psnr, train_ssim, valid_psnr, valid_ssim]
            self.summary_logging(train_output['global_step'], summary_format, summary_data)
            end_time = time.time()
            message = 'Epoch [%3d/%3d] Train(Valid) loss: %.4f(%.4f), T PSNR(MSE) %s(%s), V PSNR(MSE) %s(%s), time %s' % (
                epoch_cnt, self.epochs, Tloss, Vloss, train_psnr, train_mse, valid_psnr, valid_mse, np.mean(validation_time))
            self.logger.info(message)
            
            
            if epoch_cnt%50 == 0: 
                model_filename = self.save_model(self.saver, epoch_cnt, Vloss)
                self.logger.info('Val loss decrease from %.4f to %.4f, saving to %s' % (min_val_loss, Vloss, model_filename))
                
            if Vloss <= min_val_loss:
                model_filename = self.save_model(self.saver, epoch_cnt, Vloss)
                self.logger.info('Val loss decrease from %.4f to %.4f, saving to %s' % (min_val_loss,Vloss, model_filename))
                min_val_loss,wait = Vloss,0
            else:
                wait += 1
                if wait > self.patience:
                    model_filename = self.save_model(self.saver, epoch_cnt, Vloss)
                    self.logger.info('Val loss decrease from %.4f to %.4f, saving to %s' % (min_val_loss,Vloss, model_filename))
                    self.logger.warn('Early stopping at epoch: %d' % epoch_cnt)
                    break
            
            epoch_cnt += 1
            sys.stdout.flush()
            
    def test(self,test_data_name):
        print("Testing Started")
        self.restore()
        start_time = time.time()
        test_fetches = {'global_step': tf.train.get_or_create_global_step(),
                        'pred_orig':   self.Decoder_valid.decoded_image
                        #'metrics':     self.Decoder_test.metrics,'loss':        self.Decoder_test.loss
                       }
        k = self.test_size*self.batch_size
        pred = np.zeros((k,512,512,self.nF))
        matcontent_v = {}
        for tested_batch in range(self.test_size):
            (measure_test,mask_train,ground_test) = self.gen_test.__next__()
            feed_dict_test = {self.meas_sample: measure_test,self.sense_matrix: mask_train,self.truth_seg: ground_test}
            test_output = self.sess.run(test_fetches,feed_dict=feed_dict_test)
            pred[tested_batch*self.batch_size:(tested_batch+1)*self.batch_size,:,:,:]=test_output['pred_orig']
        ## data processing and show
        
        rsize,csize =20, 4*self.nF/5
        for m in range(k):
            temp1 = np.squeeze(pred[m,:,:,:])
            temp2 = ndimage.rotate(temp1, -135, axes=(0, 1), reshape=True)
            temp3 = temp2[181:182+363,181:182+363,:];
            plt.figure(figsize=(rsize,csize))
            for i in range(self.nF):
                plt.subplot(self.nF/5,5,i+1)
                plt.imshow(temp3[:,:,i],cmap='gray',vmin=0,vmax=1)
                plt.axis('off')
                plt.title('frame%d'%(i+1))
            plt.show()

        matcontent_v[u'pred']= pred
        hdf5storage.write(matcontent_v,'.',(self.log_dir+'/Recon_'+test_data_name+'_Cr%d'+'.mat') %(self.nF),store_python_metadata=False,matlab_compatible=True)

        print("Testing Finished")
    
        
    def calculate_scheduled_lr(self, epoch, min_lr=1e-8):
        decay_factor = int(math.ceil((epoch - self.lr_decay_epoch) / float(self.lr_decay_interval)))
        new_lr = self.lr_init * (self.lr_decay_coe ** max(0, decay_factor))
        new_lr = max(min_lr, new_lr)
        
        self.logger.info('Current learning rate to: %.6f' % new_lr)
        sys.stdout.flush()
        
        self.sess.run(self.lr_update, feed_dict={self.lr_new: new_lr})
        self.Decoder_train.set_lr(self.learning_rate) 
        return new_lr
