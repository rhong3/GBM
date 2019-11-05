"""
Simple model for GBM

Created on 11/05/2019

@author: RH
"""
import tensorflow as tf
from keras.layers.convolutional import Conv2D
from keras.layers.pooling import MaxPooling2D, AveragePooling2D, GlobalAveragePooling2D
from keras.layers.core import Dense, Dropout, Flatten, Activation, Lambda
from keras.layers.normalization import BatchNormalization
from keras.layers.merge import concatenate, add
from keras.regularizers import l2


def resnet_v2_stem(input):
    '''The stem of the pure Inception-v4 and Inception-ResNet-v2 networks. This is input part of those networks.'''

    # Input shape is 299 * 299 * 3 (Tensorflow dimension ordering)
    x = Conv2D(32, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", strides=(2, 2))(input)  # 149 * 149 * 32
    x = Conv2D(32, (3, 3), kernel_regularizer=l2(0.0002), activation="relu")(x)  # 147 * 147 * 32
    x = Conv2D(64, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", padding="same")(x)  # 147 * 147 * 64

    x1 = MaxPooling2D((3, 3), strides=(2, 2))(x)
    x2 = Conv2D(96, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", strides=(2, 2))(x)

    x = concatenate([x1, x2], axis=3)  # 73 * 73 * 160

    x1 = Conv2D(64, (1, 1), kernel_regularizer=l2(0.0002), activation="relu", padding="same")(x)
    x1 = Conv2D(96, (3, 3), kernel_regularizer=l2(0.0002), activation="relu")(x1)

    x2 = Conv2D(64, (1, 1), kernel_regularizer=l2(0.0002), activation="relu", padding="same")(x)
    x2 = Conv2D(64, (7, 1), kernel_regularizer=l2(0.0002), activation="relu", padding="same")(x2)
    x2 = Conv2D(64, (1, 7), kernel_regularizer=l2(0.0002), activation="relu", padding="same")(x2)
    x2 = Conv2D(96, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x2)

    x = concatenate([x1, x2], axis=3)  # 71 * 71 * 192

    x1 = Conv2D(192, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", strides=(2, 2))(x)

    x2 = MaxPooling2D((3, 3), strides=(2, 2))(x)

    x = concatenate([x1, x2], axis=3)  # 35 * 35 * 384

    x = BatchNormalization(axis=3)(x)
    x = Activation("relu")(x)

    return x


def simplefusionv1(input, demographics, dropout_keep_prob=0.8, num_classes=1000,
                  is_training=True, scope='SimpleFusionV1'):
    x = resnet_v2_stem(input)
    x = Conv2D(150, (5, 5), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x)
    x = Conv2D(500, (5, 5), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x)
    x = Conv2D(750, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x)
    x = MaxPooling2D((3, 3), strides=(2, 2))(x)  # 17 * 17 * 750
    x = Conv2D(1000, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x)
    x = Conv2D(1250, (3, 3), kernel_regularizer=l2(0.0002), activation="relu", padding="valid")(x)
    x = MaxPooling2D((3, 3), strides=(2, 2))(x)  # 8 * 8 * 1250

    net = x

    # Average Pooling
    x = GlobalAveragePooling2D(name='avg_pool')(x)  # Output: 1250

    pool5_drop_10x10_s1 = Dropout(dropout_keep_prob)(x, training=is_training)

    demographics = Dense(5, name='demographic_fc1', activation="relu", kernel_regularizer=l2(0.0002))(demographics)

    merged = concatenate([pool5_drop_10x10_s1, demographics])

    loss3_classifier_W = Dense(num_classes, name='loss3/classifier', kernel_regularizer=l2(0.0002))

    loss3_classifier = loss3_classifier_W(merged)

    w_variables = loss3_classifier_W.get_weights()

    logits = loss3_classifier

    return logits, net, tf.convert_to_tensor(w_variables[0])
