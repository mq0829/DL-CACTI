import os
import tensorflow as tf

class Basement_TFModel(object):
    '''Define and Initialize the basic/necessary element of a tensorflow model '''
    def __init__(self, sess, config, learning_rate, is_training):

        # Initialization of General value for model building
        # For all the other component (Hyperparameter for model structure), Please refer to the function "initial_parameter"
        self.sess = sess
        self.config = config
        self.is_training = is_training
        self.model_name = config.get('model_name','MDAnalyzer')

        # Model training SetUp
        self.train_op = None
        self.learning_rate = learning_rate
        self.max_grad_norm = float(config.get('max_grad_norm', 5.0))
        
        # Loss function & Metric
        self.loss = None
        self.loss_func = config.get('loss_func','RMSE')
        self.maximum_type = int(config.get('upbound',1))

    def set_lr(self,new_learning_rate):
        self.learning_rate = new_learning_rate
    
    def save_checkpoint(self, step=None):
        #print ''
        #ckpt_file = os.path.join(self.checkpoint_dir, self.model_name)
        #self.saver.save(self.sess, ckpt_file, global_step=step)
        pass
    
    def load_checkpoint(self,step=None):
        pass
    
    def initial_parameter(self):
        pass