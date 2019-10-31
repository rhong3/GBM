"""
Main method for preparing tfrecords for Prince on tensorflow 1 code and Inception architectures

Created on 11/26/2018

@author: RH
"""
import os
import sys
import numpy as np
import tensorflow as tf
import pandas as pd
import cv2
import Sample_prep_fusion as Sample_prep
import time
import matplotlib
matplotlib.use('Agg')

dirr = sys.argv[1]  # output directory
pdmd = sys.argv[2]  # feature to predict

try:
    level = sys.argv[3]  # magnification of tiles to use
except IndexError:
    level = None

if pdmd == 'telomere':
    classes = 3
else:
    classes = 2

# paths to directories
img_dir = '../tiles/'
LOG_DIR = "../Results/{}".format(dirr)
METAGRAPH_DIR = "../Results/{}".format(dirr)
data_dir = "../Results/{}/data".format(dirr)
out_dir = "../Results/{}/out".format(dirr)


# read images
def load_image(addr):
    img = cv2.imread(addr)
    img = img.astype(np.float32)
    return img


# used for tfrecord labels generation
def _int64_feature(value):
    return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))


# used for tfrecord images generation
def _bytes_feature(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))


# loading images for dictionaries and generate tfrecords
def loader(totlist_dir, ds):
    if ds == 'train':
        slist = pd.read_csv(totlist_dir + '/tr_sample.csv', header=0)
    elif ds == 'validation':
        slist = pd.read_csv(totlist_dir + '/va_sample.csv', header=0)
    elif ds == 'test':
        slist = pd.read_csv(totlist_dir + '/te_sample.csv', header=0)
    else:
        slist = pd.read_csv(totlist_dir + '/te_sample.csv', header=0)
    imlist = slist['path'].values.tolist()
    lblist = slist['label'].values.tolist()
    wtlist = slist['weight'].values.tolist()
    tnlist = slist['percent_tumor_nuclei'].values.tolist()
    tclist = slist['percent_total_cellularity'].values.tolist()
    nclist = slist['percent_necrosis'].values.tolist()
    aglist = slist['age'].values.tolist()

    filename = data_dir + '/' + ds + '.tfrecords'
    writer = tf.python_io.TFRecordWriter(filename)
    for i in range(len(imlist)):
        if not i % 1000:
            sys.stdout.flush()
        try:
            # Load the image
            img = load_image(imlist[i])
            label = lblist[i]
            wt = wtlist[i]
            tn = tnlist[i]
            tc = tclist[i]
            nc = nclist[i]
            ag = aglist[i]
            # Create a feature
            feature = {ds + '/label': _int64_feature(label),
                       ds + '/weight': _int64_feature(wt),
                       ds + '/nuclei': _int64_feature(tn),
                       ds + '/cellularity': _int64_feature(tc),
                       ds + '/necrosis': _int64_feature(nc),
                       ds + '/age': _int64_feature(ag),
                       ds + '/image': _bytes_feature(tf.compat.as_bytes(img.tostring()))}
            # Create an example protocol buffer
            example = tf.train.Example(features=tf.train.Features(feature=feature))

            # Serialize to string and write on the file
            writer.write(example.SerializeToString())
        except AttributeError:
            print('Error image:' + imlist[i])
            pass

    writer.close()
    sys.stdout.flush()


if __name__ == "__main__":
    tf.reset_default_graph()
    # make directories if not exist
    for DIR in (LOG_DIR, METAGRAPH_DIR, data_dir, out_dir):
        try:
            os.mkdir(DIR)
        except FileExistsError:
            pass
    # get counts of testing, validation, and training datasets;
    # if not exist, prepare testing and training datasets from sampling
    try:
        tes = pd.read_csv(data_dir+'/te_sample.csv', header=0)
        vas = pd.read_csv(data_dir+'/va_sample.csv', header=0)
    except FileNotFoundError:
        alll = Sample_prep.big_image_sum(pmd=pdmd, path=img_dir)
        trs, tes, vas = Sample_prep.set_sep(alll, path=data_dir, cls=classes, level=level)
        loader(data_dir, 'train')
        loader(data_dir, 'validation')
        loader(data_dir, 'test')

