#! /usr/bin/env python

import hera_cal
import aipy
from pyuvdata import UVData
import tensorflow as tf
import numpy as np
import pylab as plt
import sys, os

print 'imported'
files = sys.argv[1:]
uv = UVData()
uv.read_miriad(files[0])
print 'make aa'
aa = hera_cal.utils.get_aa_from_uv(uv)
info = hera_cal.omni.aa_to_info(aa, pols='x', tol=1)
reds = info.get_reds()
reds = hera_cal.omni.filter_reds(reds, ex_ants=[51])

print 'read data'
uv_info, data, flgs = aipy.miriad.read_files(files, antstr='cross', polstr='xx', verbose=True)
NCHAN = uv_info['freqs'].size
NREDS = len(reds)

#train_x = np.array([np.array([data[k+('xx',)][0].real,data[k+('xx',)][0].imag]).T for r in reds for k in r[:-1]])
train_x = np.array([np.array([data[k+('xx',)].T.real,data[k+('xx',)].T.imag]).T for r in reds for k in r[:-1]])
train_x.shape = (-1, NCHAN, 2)
train_y = np.zeros((train_x.shape[0], NREDS))
for cnt,gp in enumerate([i for i,r in enumerate(reds) for k in r[:-1]]):
    train_y[cnt,gp] = 1. 

# XXX could also try another integration here
test_x = np.array([np.array([data[k+('xx',)].T.real,data[k+('xx',)].T.imag]).T for r in reds for k in r])
test_x.shape = (-1, NCHAN, 2)
test_y = np.zeros((test_x.shape[0], NREDS))
for cnt,gp in enumerate([i for i,r in enumerate(reds) for k in r]):
    test_y[cnt,gp] = 1. 

def weight_variable(shape, name=None):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial, name=name)

def bias_variable(shape, name=None):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial, name=name)

def conv1d(x, W):
  return tf.nn.conv1d(x, W, stride=1, padding='SAME')

def max_pool_4(x):
    return tf.nn.pool(x, window_shape=[4], strides=[4], pooling_type='MAX', padding='SAME')


x = tf.placeholder(tf.float32, shape=[None, NCHAN, 2]) # input layer
W_conv1 = weight_variable([32, 2, 64], name='W_conv1') # coeffs for conv, (patch_in,layers,features_out)
b_conv1 = bias_variable([64], name='b_conv1') # bias for conv, (features_out)
h_conv1 = tf.nn.relu(conv1d(x, W_conv1) + b_conv1)
h_pool1 = max_pool_4(h_conv1) # reduce spectrum size to 256
W_conv2 = weight_variable([32, 64, 128], name='W_conv2') # coeffs for conv, (patch_in,layers,features_out)
b_conv2 = bias_variable([128], name='b_conv2') # bias for conv, (features_out)
h_conv2 = tf.nn.relu(conv1d(h_pool1, W_conv2) + b_conv2)
h_pool2 = max_pool_4(h_conv2) # reduce spectrum size to 256
W_fc1 = weight_variable([NCHAN/16 * 128, 1024], name='W_fc1') # coeffs for fully connected layer
b_fc1 = bias_variable([1024], name='b_fc1') # bias for fully connected layer
h_pool2_flat = tf.reshape(h_pool2, [-1, NCHAN/16*128])
h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1) # fully connected layer
# prevent overfitting with dropout layer
keep_prob = tf.placeholder(tf.float32)
h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)
W_fc2 = weight_variable([1024, NREDS], name='W_fc2')
b_fc2 = bias_variable([NREDS], name='b_fc2')
y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2 # final voting on red membership

y_ = tf.placeholder(tf.float32, shape=[None, NREDS]) # output layer

cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)
correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

#import IPython; IPython.embed()
init_op = tf.global_variables_initializer()
saver = tf.train.Saver()
FILENAME = os.path.abspath('./ml_baseline001_model.ckpt')

with tf.Session() as sess:
  sess.run(init_op)
  if os.path.exists(FILENAME): saver.restore(sess, FILENAME)
  for i in range(20000):
    batch = np.random.choice(train_x.shape[0], 50, replace=False)
    if i % 100 == 0:
      train_accuracy = accuracy.eval(feed_dict={x: train_x[batch], y_: train_y[batch], keep_prob: 1.0})
      print('step %d, training accuracy %g' % (i, train_accuracy))
      save_path = saver.save(sess, FILENAME)
      #import IPython; IPython.embed()
    train_step.run(feed_dict={x: train_x[batch], y_: train_y[batch], keep_prob: 0.5})

  print('test accuracy %g' % accuracy.eval(feed_dict={x: test_x, y_: test_y, keep_prob: 1.0}))

  import IPython; IPython.embed()

# The convolution will compute 32 features for each 5x5 patch. 
# Its weight tensor will have a shape of [5, 5, 1, 32]. 
# The first two dimensions are the patch size, the next is the number of 
# input channels, and the last is the number of output channels. 

# To apply the layer, we first reshape x to a 4d tensor, 
# with the second and third dimensions corresponding to image width and height, 
# and the final dimension corresponding to the number of color channels.

# We then convolve x_image with the weight tensor, add the bias, apply the ReLU function, 
# and finally max pool. The max_pool_2x2 method will reduce the image size to 14x14.

# In order to build a deep network, we stack several layers of this type. 
# The second layer will have 64 features for each 5x5 patch.

# Now that the image size has been reduced to 7x7, we add a fully-connected 
# layer with 1024 neurons to allow processing on the entire image. We reshape 
# the tensor from the pooling layer into a batch of vectors, multiply 
# by a weight matrix, add a bias, and apply a ReLU.

# To reduce overfitting, we will apply dropout before the readout layer. We
# create a placeholder for the probability that a neuron's output is kept
# during dropout. This allows us to turn dropout on during training, and turn
# it off during testing. TensorFlow's tf.nn.dropout op automatically handles
# scaling neuron outputs in addition to masking them, so dropout just works
# without any additional scaling.

