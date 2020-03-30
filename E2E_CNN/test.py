########
# Please choose the scene you'd like to reconstruct in line 19. The names of measurements include 'duomino', 
# 'pendulumBall' and 'waterBalloon'.
# Please set the compression ratio(Cr) in line 18. You can set Cr to be 10, 20 or 30.
########
from __future__ import absolute_import

import tensorflow as tf
import yaml
import os
import h5py

from Model.Decoder_Handler import Decoder_Handler

config_filename = './Model/Config.yaml'

def main():
    Cr = 10            
    test_data_name = 'waterBalloon'
    if Cr ==10:
        ave_folder,ave_config = 'Cr10_model','config_Cr10_model.yaml'
    elif Cr ==20:
        ave_folder,ave_config = 'Cr20_model','config_Cr20_model.yaml'
    else:
        ave_folder,ave_config = 'Cr30_model','config_Cr30_model.yaml'
    folder_id,config_id = ave_folder,ave_config
    with open(config_filename) as handle:
        model_config = yaml.load(handle,Loader=yaml.FullLoader)  
    data_name = []
    data_name.append(os.path.join(os.path.abspath('.'), model_config['category'], model_config['data_name']))
    data_name.append(os.path.join(os.path.abspath('.'), model_config['category_valid'], model_config['data_name']))
    data_name.append(os.path.join(os.path.abspath('..'), 'dataset', 'meas_'+test_data_name))
    log_dir = os.path.join(os.path.abspath('.'),model_config['result_dir'],model_config['result_model'],folder_id)

    with open(os.path.join(log_dir, config_id)) as handle:
        model_config = yaml.load(handle,Loader=yaml.FullLoader)

    mask_name = os.path.join(os.path.abspath('..'), 'dataset','mask')
        
    dataset_name = (data_name,mask_name)
    
    tf_config = tf.ConfigProto()
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    tf_config = tf.ConfigProto()
    tf_config.gpu_options.allow_growth = True

    with tf.Session(config=tf_config) as sess:
        Cube_Decoder = Decoder_Handler(dataset_name=dataset_name, model_config=model_config, sess = sess, is_training=False, Cr=Cr)
        Cube_Decoder.test(test_data_name)

if __name__ == '__main__':
    main()
    
    
    