from __future__ import absolute_import

import tensorflow as tf
import yaml
import os
import h5py

from Model.Decoder_Handler import Decoder_Handler


config_filename = './Model/Config.yaml'

def main():
    with open(config_filename) as handle:
        model_config = yaml.load(handle,Loader=yaml.FullLoader)
    data_name = []
    data_name.append(os.path.join(os.path.abspath('.'), model_config['category'], model_config['data_name']))
    data_name.append(os.path.join(os.path.abspath('.'), model_config['category_valid'], model_config['data_name']))
    data_name.append(os.path.join(os.path.abspath('.'), model_config['category_test'], model_config['data_name']))

    mask_name = os.path.join(os.path.abspath('.'), model_config['category'],model_config['mask_name'])
        
    dataset_name = (data_name,mask_name)
    
    tf_config = tf.ConfigProto()
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    # tf.device("GPU:1")
    tf_config = tf.ConfigProto()
    tf_config.gpu_options.allow_growth = True

    with tf.Session(config=tf_config) as sess:
        Cube_Decoder = Decoder_Handler(dataset_name=dataset_name, model_config=model_config, sess = sess, is_training=True)
        Cube_Decoder.train()

if __name__ == '__main__':
    main()
