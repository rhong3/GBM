"""
Prepare training and testing datasets as CSV dictionaries

Created on 10/30/2019

@author: RH
"""
import os
import pandas as pd
import sklearn.utils as sku
import numpy as np


def tile_ids_in(inp):
    ids = []
    try:
        for id in os.listdir(inp['path']):
            if '_{}.png'.format(str(inp['sldnum'])) in id:
                ids.append([inp['slide'], inp['level'], inp['path']+'/'+id, inp['weight'], inp['percent_tumor_nuclei'],
                            inp['percent_total_cellularity'], inp['percent_necrosis'], inp['age'], inp['label']])
    except FileNotFoundError:
        print('Ignore:', inp['path'])

    return ids


# Get all svs images with its label as one file; level is the tile resolution level
def big_image_sum(pmd, path='../tiles/', ref_file='../feature_summary.csv'):
    ref = pd.read_csv(ref_file, sep=',', header=0)
    ref = ref.loc[ref['used_in_proteome'] == True]
    ref = ref.rename(columns={pmd: 'label'})
    ref = ref.dropna(subset=['label'])
    ref['sldnum'] = ref['slide_id'].str.split("-", n=2, expand=True)[-1]
    ref = ref[['case_id', 'sldnum', 'weight', 'percent_tumor_nuclei', 'percent_total_cellularity', 'percent_necrosis',
               'age', 'label']]
    ref = ref.rename(columns={'case_id': 'slide'})
    ref1 = ref
    ref2 = ref
    ref['level'] = 0
    ref['path'] = path + "{}/level{}".format(ref['pctnum'], ref['level'])
    ref1['level'] = 1
    ref1['path'] = path + "{}/level{}".format(ref1['pctnum'], ref1['level'])
    ref2['level'] = 2
    ref2['path'] = path + "{}/level{}".format(ref2['pctnum'], ref2['level'])
    datapd = pd.concat([ref, ref1, ref2])

    return datapd


# seperate into training and testing; each type is the same separation ratio on big images
# test and train csv files contain tiles' path.
def set_sep(alll, path, cls, level=None, cut=0.3, batchsize=64):
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
        tile_ids = tile_ids_in(row)
        test_tiles_list.extend(tile_ids)
    for idx, row in train.iterrows():
        tile_ids = tile_ids_in(row)
        train_tiles_list.extend(tile_ids)
    for idx, row in validation.iterrows():
        tile_ids = tile_ids_in(row)
        validation_tiles_list.extend(tile_ids)
    test_tiles = pd.DataFrame(test_tiles_list, columns=['slide', 'level', 'path', 'weight', 'percent_tumor_nuclei',
                                                        'percent_total_cellularity', 'percent_necrosis', 'age',
                                                        'label'])
    train_tiles = pd.DataFrame(train_tiles_list, columns=['slide', 'level', 'path', 'weight', 'percent_tumor_nuclei',
                                                          'percent_total_cellularity', 'percent_necrosis', 'age',
                                                          'label'])
    validation_tiles = pd.DataFrame(validation_tiles_list, columns=['slide', 'level', 'path', 'weight',
                                                                    'percent_tumor_nuclei', 'percent_total_cellularity',
                                                                    'percent_necrosis', 'age', 'label'])

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

