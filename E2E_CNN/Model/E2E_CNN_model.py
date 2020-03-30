import tensorflow as tf
import tensorflow.contrib.layers as layers
import tensorflow.contrib.slim as slim
import numpy as np

from E2E_CNN.Lib.Utility import *
from E2E_CNN.Model.Base_TFModel import Basement_TFModel

class Depth_Decoder(Basement_TFModel):
    
    def __init__(self, value_sets, init_learning_rate, sess, config, is_training=True, *args, **kwargs):
        
        super(Depth_Decoder, self).__init__(sess=sess, config=config, learning_rate=init_learning_rate,is_training=is_training)

        (measurement, mat_sense, truth_seg) = value_sets
        #print(measurement.shape, mat_sense.shape, truth_seg.shape)
        self.depth = mat_sense.get_shape().as_list()[-1]
        self.mask = None
        
        # Initialization of the model hyperparameter, enc-dec structure, evaluation metric & Optimizier
        self.initial_parameter()
        self.decoded_image = self.encdec_handler(measurement, mat_sense)
        self.metric_opt(self.decoded_image, truth_seg)
        
    def encdec_handler(self, measurement, mat_sense):

        self.hyper_structure = [(2,64,3,2),(2,64,3,2),(2,64,3,2),(2,64,3,2),(2,64,3,2)]
        self.end_encoder = (2,64,3)
        
        encoder_in = measurement
        #print(encoder_in.get_shape().as_list())
        output = self.inference(encoder_in,0.8,phase_train = True)    #self.is_training
        return output
    
    def inference(self, images, keep_probability,phase_train=True, bottleneck_layer_size=128, weight_decay=0.0005, reuse=None):
        batch_norm_params = {
            # Decay for the moving averages.
            'decay': 0.995,
            # epsilon to prevent 0s in variance.
            'epsilon': 0.001,
            # force in-place updates of mean and variance estimates
            'updates_collections': None,
            'scale':True,
            'is_training':phase_train,
            # Moving averages ends up in the trainable variables collection
            'variables_collections': [tf.GraphKeys.TRAINABLE_VARIABLES],}
        
        with slim.arg_scope([slim.conv2d, slim.fully_connected,slim.conv2d_transpose],
                            weights_initializer=slim.initializers.xavier_initializer(),
                            weights_regularizer=slim.l2_regularizer(weight_decay),
                            normalizer_fn=slim.batch_norm,normalizer_params=batch_norm_params):
            return self.encoder_decoder(images, is_training=phase_train,dropout_keep_prob=keep_probability,reuse=reuse)
    

    def EncConv_module(self, net, module_ind, hyper_struc, PoolValid=True):
        (lnum,knum,ksize,pstr) = hyper_struc
        temp = net
        net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME',scope='en_%d_%d'%(module_ind,lnum-1))
        net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME', activation_fn=None, scope='en_%d_%d'%(module_ind,lnum))
        net = net+temp
        net = tf.nn.relu(net)
        self.end_points['encode_%d'%(module_ind)] = net
        #print(net.get_shape().as_list())
        if PoolValid is True:
            return slim.max_pool2d(net,pstr,stride=pstr,padding='SAME',scope='Pool%d'%(module_ind))
        else:
            return net

    def DecConv_module(self, net, module_ind, hyper_struc, PoolValid=True):
        (lnum,knum,ksize,pstr) = hyper_struc
        if PoolValid is True:
            net = slim.conv2d_transpose(net, knum, pstr, pstr, padding='SAME')
        #print(net.get_shape().as_list())
        net=tf.add(net,self.end_points['encode_%d'%(module_ind)])
        temp = net
        net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME',scope='de_%d_%d'%(module_ind,lnum-1))
        net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME', activation_fn=None, scope='de_%d_%d'%(module_ind,lnum))
        net = net+temp
        net = tf.nn.relu(net)
        return net
            
    def encoder_decoder(self, inputs, is_training=True, dropout_keep_prob=0.8, reuse=None, scope='generator'):
        self.end_points = {}
        with tf.variable_scope(scope, 'generator', [inputs], reuse=reuse):
            with slim.arg_scope([slim.batch_norm, slim.dropout],is_training=is_training):
                with slim.arg_scope([slim.conv2d, slim.max_pool2d, slim.avg_pool2d],stride=1, padding='SAME'):
                    ############################# encoder ##############################################
                    self.end_points['inputs'] = inputs
                    net = slim.conv2d(inputs, 64, 3, stride=1, padding='SAME', scope='en_0')
                    net = self.EncConv_module(net,1,self.hyper_structure[0],PoolValid=False)
                    net = self.EncConv_module(net,2,self.hyper_structure[1],PoolValid=False)
                    net = self.EncConv_module(net,3,self.hyper_structure[2],PoolValid=False)
                    net = self.EncConv_module(net,4,self.hyper_structure[3],PoolValid=False)
                    net = self.EncConv_module(net,5,self.hyper_structure[4],PoolValid=False)
                    (lnum,knum,ksize) = self.end_encoder

                    net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME', scope='en_7')
                    net = slim.conv2d(net, knum, ksize, stride=1, padding='SAME', scope='en_8')
                    net = self.DecConv_module(net, 5, self.hyper_structure[4],PoolValid=False)
                    ############################# decoder ##############################################
                    net = self.DecConv_module(net,4,self.hyper_structure[3],PoolValid=False)
                    net = self.DecConv_module(net,3,self.hyper_structure[2],PoolValid=False)
                    net = self.DecConv_module(net,2,self.hyper_structure[1],PoolValid=False)
                    net = self.DecConv_module(net,1,self.hyper_structure[0],PoolValid=False)
                    
                    net=slim.conv2d(net,self.depth,1,stride=1,padding='SAME',activation_fn=tf.nn.tanh)#tf.nn.sigmoid
                    net = tf.add(net,self.end_points['inputs'])
        return net
    
    def metric_opt(self, model_output, ground_truth):
        
        if self.loss_func == 'MSE':
            self.loss = loss_mse(model_output, ground_truth, self.mask)
        elif self.loss_func == 'RMSE':
            self.loss = loss_rmse(model_output, ground_truth, self.mask)#+ 0.7*loss_SSIM(model_output, ground_truth, self.mask)
        elif self.loss_func == 'MAE':
            self.loss = loss_mae(model_output, ground_truth, self.mask)
        elif self.loss_func == 'SSIM':
            self.loss = loss_SSIM(model_output, ground_truth, self.mask)
        else:
            self.loss = loss_rmse(model_output, ground_truth, self.mask)
            
        self.metrics = calculate_metrics(model_output, ground_truth, self.mask)
        global_step = tf.train.get_or_create_global_step()
            
        if self.is_training:
            optimizer = tf.train.AdamOptimizer(self.learning_rate)
            tvars = tf.trainable_variables()
            grads = tf.gradients(self.loss, tvars)
            grads, _ = tf.clip_by_global_norm(grads, self.max_grad_norm)
            self.train_op = optimizer.apply_gradients(zip(grads, tvars), global_step=global_step, name='train_op')
        self.info_merge = tf.summary.merge_all()
             
    def initial_parameter(self):

        config = self.config
        # Parameter Initialization of Data Assignment
        self.num_heads = int(config.get('num_heads',10))
        self.batch_size = int(config.get('batch_size',12))
        self.hcube,self.wcube = int(config.get('cube_height',50)),int(config.get('cube_width',50))
        value_units,weight_units = int(config.get('dim_value',20)),int(config.get('dim_weight',20)) 
        self.att_unit = (value_units, weight_units)
        self.num_space,self.num_spec = int(config.get('num_spatial',100)),int(config.get('num_spectral',201))
            
        
        # Parameter Initialization of Model Framework
        self.ResLayer = int(config.get('Depth_Residual',28))
        self.atte_lcoe = int(config.get('atte_length_coe',8))
        self.model_structure = self.config.get('model_structure')
        
        # label of mask
        self.flag_identity = config.get('flag_identity',False)
        
 