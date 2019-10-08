"""
Prepare training and testing datasets as CSV dictionaries (Further modification required for GBM)

Created on 11/26/2018

@author: RH
"""
import os
import pandas as pd
import sklearn.utils as sku
import numpy as np

tile_path = "../tiles/"


# get all full paths of images
def image_ids_in(root_dir, ignore=['.DS_Store','dict.csv', 'all.csv']):
    ids = []
    for id in os.listdir(root_dir):
        if id in ignore:
            print('Skipping ID:', id)
        else:
            ids.append(id)
    return ids


# Get intersection of 2 lists
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def tile_ids_in(slide, level, root_dir, label, sldnum):
    ids = []
    try:
        for id in os.listdir(root_dir):
            if '_{}.png'.format(str(sldnum)) in id:
                ids.append([slide, level, root_dir+'/'+id, label])
            else:
                print('Skipping ID:', id)
    except FileNotFoundError:
        print('Ignore:', root_dir)

    return ids


# Get all svs images with its label as one file; level is the tile resolution level
def big_image_sum(pmd, path='../tiles/', dict_file='../tcia_pathology_slides.tsv', ref_file='../dummy_His_MUT_joined.csv'):
    dct = pd.read_csv(dict_file, sep='\t', header=0)
    dct = dct.loc[dct['used_in_proteome'] == 'TRUE']
    ref = pd.read_csv(ref_file, header=0)
    big_images = []
    for level in range(3):
        level = str(level)
        if pmd == 'subtype':
            print(None)
        else:
            negimg = intersection(ref.loc[ref[pmd] == 0]['name'].tolist(), dct['case_id'].tolist())
            negsld = dct[dct['case_id'].isin(negimg)]['slide_id'].tolist()
            posimg = intersection(ref.loc[ref[pmd] == 1]['name'].tolist(), dct['case_id'].tolist())
            possld = dct[dct['case_id'].isin(posimg)]['slide_id'].tolist()
            for i in negsld:
                sldnum = i.split('-')[-1]
                pctnum = i[:-4]
                big_images.append([i, level, path + "{}/level{}".format(pctnum, level), sldnum, 0])
            for i in possld:
                sldnum = i.split('-')[-1]
                pctnum = i[:-4]
                big_images.append([pctnum, level, path + "{}/level{}".format(pctnum, level), sldnum, 1])

    datapd = pd.DataFrame(big_images, columns=['slide', 'level', 'path', 'sldnum', 'label'])

    return datapd


# seperate into training and testing; each type is the same separation ratio on big images
# test and train csv files contain tiles' path.
def set_sep(alll, path, cls, level=None, cut=0.2, batchsize=64):
    trlist = []
    telist = []
    valist = []
    if level:
        alll = alll[alll.level == level]

    CPTAC = alll
    for i in range(cls):
        subset = CPTAC.loc[CPTAC['label'] == i]
        unq = list(subset.slide.unique())
        np.random.shuffle(unq)
        validation = unq[:int(len(unq) * cut / 2)]
        valist.append(subset[subset['slide'].isin(validation)])
        test = unq[int(len(unq) * cut / 2):int(len(unq) * cut)]
        telist.append(subset[subset['slide'].isin(test)])
        train = unq[int(len(unq) * cut):]
        trlist.append(subset[subset['slide'].isin(train)])

    test = pd.concat(telist)
    train = pd.concat(trlist)
    validation = pd.concat(valist)
    test_tiles_list = []
    train_tiles_list = []
    validation_tiles_list = []
    for idx, row in test.iterrows():
        tile_ids = tile_ids_in(row['slide'], row['level'], row['path'], row['label'], row['sldnum'])
        test_tiles_list.extend(tile_ids)
    for idx, row in train.iterrows():
        tile_ids = tile_ids_in(row['slide'], row['level'], row['path'], row['label'], row['sldnum'])
        train_tiles_list.extend(tile_ids)
    for idx, row in validation.iterrows():
        tile_ids = tile_ids_in(row['slide'], row['level'], row['path'], row['label'], row['sldnum'])
        validation_tiles_list.extend(tile_ids)
    test_tiles = pd.DataFrame(test_tiles_list, columns=['slide', 'level', 'path', 'label'])
    train_tiles = pd.DataFrame(train_tiles_list, columns=['slide', 'level', 'path', 'label'])
    validation_tiles = pd.DataFrame(validation_tiles_list, columns=['slide', 'level', 'path', 'label'])

    # No shuffle on test set
    train_tiles = sku.shuffle(train_tiles)
    validation_tiles = sku.shuffle(validation_tiles)
    if train_tiles.shape[0] > int(batchsize*30000):
        train_tiles = train_tiles.sample(int(batchsize*30000), replace=False)
        print('Truncate training set!')
    if validation_tiles.shape[0] > int(batchsize*3000):
        validation_tiles = validation_tiles.sample(int(batchsize*3000), replace=False)
        print('Truncate validation set!')
    if test_tiles.shape[0] > int(batchsize*30000):
        test_tiles = test_tiles.sample(int(batchsize*30000), replace=False)
        print('Truncate test set!')

    test_tiles.to_csv(path+'/te_sample.csv', header=True, index=False)
    train_tiles.to_csv(path+'/tr_sample.csv', header=True, index=False)
    validation_tiles.to_csv(path+'/va_sample.csv', header=True, index=False)

    return train_tiles, test_tiles, validation_tiles

