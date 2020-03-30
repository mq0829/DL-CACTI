########
# 'test_E2ECNN.py' tests the trained convolutional neural network 
# for video reconstruction in 'coded aperture compressive temporal imaging (CACTI)'

# Reference
#   [1] M. Qiao, Z. Meng, J. Ma, X. Yuan, Deep learning for video compressive
#       sensing, APL Photonics 5, 030801 (2020).
#   [2] Yuan, Xin. "Generalized alternating projection based total variation minimization for compressive sensing." 
#       2016 IEEE International Conference on Image Processing (ICIP). IEEE, 2016.


# Contact
#   Xin Yuan, Bell Labs, xyuan@bell-labs.com
#   Mu Qiao, New Jersey Institute of Technology, muqiao@njit.edu
#   Update Mar 13, 2020.

# Test steps
# [0] Specify 'test_data_name' in line 42 as 'waterBalloon', 'hand', 'duomino', or 'pendulumBall', as found in the file 'dataset'
# [1] Specify 'Cr' (compression ratio) in line 41 as 10, 20 or 30
# [2] Run the code
# [3] Results will be stored in 'E2E_CNN_algorithm/Result/Validation-Result'


# Environment requirement
# [0] Tensorflow-gpu==1.13.1 (conda install tensorflow-gpu=1.13.1)
# [1] Packages: numpy, yaml, scipy, hdf5storage, matplotlib, math

########
from __future__ import absolute_import

import tensorflow as tf
import yaml
import os
import h5py

from E2E_CNN.Model.Decoder_Handler import Decoder_Handler

config_filename = './E2E_CNN/Model/Config.yaml'

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
    data_name.append(os.path.join(os.path.abspath('.'), 'dataset', 'meas_'+test_data_name))
    log_dir = os.path.join(os.path.abspath('.'),model_config['result_dir'],model_config['result_model'],folder_id)

    with open(os.path.join(log_dir, config_id)) as handle:
        model_config = yaml.load(handle,Loader=yaml.FullLoader)

    mask_name = os.path.join(os.path.abspath('.'), 'dataset','mask')
        
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
    
    
