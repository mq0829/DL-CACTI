import tensorflow as tf
import numpy as np
from E2E_CNN.Lib.ms_ssim import *

    
def loss_mse(decoded, ground, label=None):
    loss_pixel = tf.square(tf.subtract(decoded, ground))
    if label is None:
        return tf.reduce_mean(loss_pixel)
    else:
        label /= tf.reduce_mean(label,axis=(1,2,3),keep_dims=True)
        loss_pixel = tf.multiply(loss_pixel,label)
        return tf.reduce_mean(loss_pixel)
    
def loss_rmse(decoded, ground, label=None):
    loss_pixel = tf.square(tf.subtract(decoded, ground))
    if label is None:
        loss_pixel = tf.sqrt(tf.reduce_mean(loss_pixel,axis=(1,2,3)))
        return tf.reduce_mean(loss_pixel)
    else:
        label /= tf.reduce_mean(label,axis=(1,2,3),keep_dims=True)
        loss_pixel = tf.multiply(loss_pixel,label)
        loss_pixel = tf.sqrt(tf.reduce_mean(loss_pixel,axis=(1,2,3)))
        return tf.reduce_mean(loss_pixel)
    
def tensor_log10(x):
    numerator = tf.log(x)
    denominator = tf.log(tf.constant(10, dtype=numerator.dtype))
    return  tf.divide(numerator,denominator)

def metric_psnr(decoded, ground, label=None):
    loss_pixel = tf.reduce_mean(tf.square(tf.subtract(decoded, ground)),axis=(1,2,3))
    psnr_s = tf.constant(-10.0)*tensor_log10(loss_pixel)
    return tf.reduce_mean(psnr_s)
'''
def metric_psnr(decoded, ground):
    loss_pixel = tf.reduce_mean(tf.square(tf.subtract(decoded, ground)),axis=(1,2,3))
    psnr_s = tf.constant(10.0)*tensor_log10(tf.square(tf.reduce_max(ground,axis=(1,2,3)))/loss_pixel)
    return tf.reduce_mean(psnr_s)
'''
def metric_ssim(decoded, ground, label=None):
    ssim_s = MultiScaleSSIM(decoded,ground)
    return tf.reduce_mean(ssim_s)

def calculate_metrics(decoded, ground, label):
    psnr = metric_psnr(decoded, ground, label) 
    ssim = metric_ssim(decoded, ground, label)
    mse  = loss_mse(decoded, ground, label)
    return psnr, ssim, mse

def loss_SSIM(decoded,ground,label=None):
    #return tf.constant(1.0)-tf.image.ssim_multiscale(decoded,ground,1.0)
    return tf.abs(tf.constant(1.0)-MultiScaleSSIM(decoded,ground))
    
def loss_mae(decoded,ground,label=None):
    loss_pixel = tf.abs(tf.subtract(decoded, ground))
    if label is None:
        return tf.reduce_mean(loss_pixel)
    else:
        label /= tf.reduce_mean(label,axis=(1,2,3),keep_dims=True)
        loss_pixel = tf.multiply(loss_pixel,label)
        return tf.reduce_mean(loss_pixel)

def loss_TV(decoded):
    TV_loss = tf.image.total_variation(decoded,name=None)
    return tf.reduce_mean(TV_loss)

def loss_spec(decoded, ground):
    grad_ground = tf.subtract(ground[:,:,:,1:],ground[:,:,:,:-1])
    grad_decode = tf.subtract(decoded[:,:,:,1:],decoded[:,:,:,:-1])
    loss_pixel = tf.square(tf.subtract(grad_ground, grad_decode))
    return tf.reduce_mean(loss_pixel)

def loss_grad(decoded, ground):
    grad_truth = tf.subtract(ground[:,:,1:,:],ground[:,:,:-1,:])
    grad_pred  = tf.subtract(decoded[:,:,1:,:],decoded[:,:,:-1,:])
    grad_label = tf.where(tf.equal(grad_truth,0), tf.zeros_like(grad_truth), tf.ones_like(grad_truth))
    grad_label /= tf.reduce_mean(grad_label,axis=(1,2,3),keepdims=True)
    diff_grad = tf.multiply(tf.square(grad_truth-grad_pred),grad_label)/tf.reduce_max(tf.abs(ground),axis=(1,2,3),keepdims=True)
    return tf.reduce_mean(diff_grad)
