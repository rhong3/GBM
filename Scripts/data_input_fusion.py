"""
Data input preparation from decoding TFrecords, onehot encoding, augmentation, and batching 2.0

Created on 10/30/2019

@author: RH
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tensorflow as tf


class DataSet(object):
    # bs is batch size; ep is epoch; images are images; mode is test/train; filename is tfrecords
    def __init__(self, bs, count, ep=1, cls=2, images=None, mode=None, filename=None):
        self._batchsize = bs
        self._index_in_epoch = 0
        self._num_examples = count
        self._images = images
        self._mode = mode
        self._filename = filename
        self._epochs = ep
        self._classes = cls

    # decoding tfrecords; return images and labels
    def decode(self, serialized_example):
        features = tf.parse_example(
            serialized_example,
            # Defaults are not specified since both keys are required.
            features={self._mode + '/image': tf.FixedLenFeature([], tf.string),
                      self._mode + '/label': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/weight': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/nuclei': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/cellularity': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/necrosis': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/age': tf.FixedLenFeature([], tf.int64), })

        image = tf.decode_raw(features[self._mode + '/image'], tf.float32)
        image = tf.reshape(image, [-1, 299, 299, 3])

        # Convert label from a scalar uint8 tensor to an int32 scalar.
        label = tf.cast(features[self._mode + '/label'], tf.int32)
        weight = tf.cast(features[self._mode + '/weight'], tf.int32)
        nuclei = tf.cast(features[self._mode + '/nuclei'], tf.int32)
        cellularity = tf.cast(features[self._mode + '/cellularity'], tf.int32)
        necrosis = tf.cast(features[self._mode + '/necrosis'], tf.int32)
        age = tf.cast(features[self._mode + '/age'], tf.int32)
        demographic = tf.concat([weight, nuclei, cellularity, necrosis, age], -1)

        return image, label, demographic

    # decoding tfrecords for real test
    def Real_decode(self, serialized_example):
        features = tf.parse_example(
            serialized_example,
            # Defaults are not specified since both keys are required.
            features={self._mode + '/image': tf.FixedLenFeature([], tf.string),
                      self._mode + '/weight': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/nuclei': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/cellularity': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/necrosis': tf.FixedLenFeature([], tf.int64),
                      self._mode + '/age': tf.FixedLenFeature([], tf.int64), })

        image = tf.decode_raw(features[self._mode + '/image'], tf.float32)
        image = tf.reshape(image, [-1, 299, 299, 3])

        weight = tf.cast(features[self._mode + '/weight'], tf.int32)
        nuclei = tf.cast(features[self._mode + '/nuclei'], tf.int32)
        cellularity = tf.cast(features[self._mode + '/cellularity'], tf.int32)
        necrosis = tf.cast(features[self._mode + '/necrosis'], tf.int32)
        age = tf.cast(features[self._mode + '/age'], tf.int32)
        demographic = tf.concat([weight, nuclei, cellularity, necrosis, age], -1)

        return image, demographic

    # augmentation including onehot encoding
    def augment(self, images, labels, demographics):

        angles = tf.cast(tf.random_uniform([], 0, 4), tf.int32)
        images = tf.image.rot90(images, k=angles)
        images = tf.image.random_flip_left_right(images)
        images = tf.image.random_flip_up_down(images)
        images = tf.image.random_hue(images, 0.02)
        images = tf.image.random_brightness(images, 0.02)
        images = tf.image.random_contrast(images, 0.9, 1.1)
        images = tf.image.random_saturation(images, 0.9, 1.1)

        labels = tf.one_hot(indices=tf.cast(labels, tf.int32), depth=self._classes)

        return images, labels, demographics

    # onehot encoding only; for test set
    def onehot_only(self, images, labels, demographic):
        with tf.name_scope('onehot_only'):
            labels = tf.one_hot(indices=tf.cast(labels, tf.int32), depth=self._classes)
        return images, labels, demographic

    # dataset preparation; batching; Real test or not; train or test
    def data(self, Not_Realtest=True, train=True):
        batch_size = self._batchsize
        ep = self._epochs
        filenames = tf.placeholder(tf.string, shape=None)
        dataset = tf.data.TFRecordDataset(filenames)
        dataset = dataset.repeat(ep)
        if Not_Realtest:
            if train:
                batched_dataset = dataset.batch(batch_size, drop_remainder=True)
                batched_dataset = batched_dataset.map(self.decode)
                batched_dataset = batched_dataset.map(self.augment)
            else:
                batched_dataset = dataset.batch(batch_size, drop_remainder=False)
                batched_dataset = batched_dataset.map(self.decode)
                batched_dataset = batched_dataset.map(self.onehot_only)
            iterator = batched_dataset.make_initializable_iterator()
            return iterator, self._filename, filenames
        else:
            batched_dataset = dataset.batch(batch_size, drop_remainder=False)
            batched_dataset = batched_dataset.map(self.Real_decode)
        iterator = batched_dataset.make_initializable_iterator()
        return iterator, self._filename, filenames

    @property
    def images(self):
        return self._images

    @property
    def num_examples(self):
        return self._num_examples
