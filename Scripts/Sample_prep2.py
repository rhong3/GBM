"""
Prepare training and testing datasets as CSV dictionaries 2.0 (Further modification required for GBM)

Created on 04/26/2019

@author: RH
"""
import os
import pandas as pd
import sklearn.utils as sku
import numpy as np
import re

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


# pair tiles of 20x, 10x, 5x of the same area
def paired_tile_ids_in(slide, label, root_dir, sldnum):
    dira = os.path.isdir(root_dir + 'level0')
    dirb = os.path.isdir(root_dir + 'level1')
    dirc = os.path.isdir(root_dir + 'level2')
    if dira and dirb and dirc:
        fac = 500
        ids = []
        for level in range(3):
            dirr = root_dir + 'level{}'.format(str(level))
            for id in os.listdir(dirr):
                if '_{}.png'.format(str(sldnum)) in id:
                    x = int(float(id.split('x-', 1)[1].split('-', 1)[0]) / fac)
                    y = int(float(re.split('.p| |_', id.split('y-', 1)[1])[0]) / fac)
                    try:
                        dup = int(re.split('.p', re.split('_', id.split('y-', 1)[1])[1])[0])
                    except IndexError:
                        dup = np.nan
                    ids.append([slide, label, level, dirr + '/' + id, x, y, dup])

        ids = pd.DataFrame(ids, columns=['slide', 'label', 'level', 'path', 'x', 'y', 'dup'])
        idsa = ids.loc[ids['level'] == 0]
        idsa = idsa.drop(columns=['level'])
        idsa = idsa.rename(index=str, columns={"path": "L0path"})
        idsb = ids.loc[ids['level'] == 1]
        idsb = idsb.drop(columns=['slide', 'label', 'level'])
        idsb = idsb.rename(index=str, columns={"path": "L1path"})
        idsc = ids.loc[ids['level'] == 2]
        idsc = idsc.drop(columns=['slide', 'label', 'level'])
        idsc = idsc.rename(index=str, columns={"path": "L2path"})
        idsa = pd.merge(idsa, idsb, on=['x', 'y', 'dup'], how='left', validate="many_to_many")
        idsa['x'] = idsa['x'] - (idsa['x'] % 2)
        idsa['y'] = idsa['y'] - (idsa['y'] % 2)
        idsa = pd.merge(idsa, idsc, on=['x', 'y', 'dup'], how='left', validate="many_to_many")
        idsa = idsa.drop(columns=['x', 'y', 'dup'])
        idsa = idsa.dropna()
        idsa = sku.shuffle(idsa)
    else:
        idsa = pd.DataFrame(columns=['slide', 'label', 'L0path', 'L1path', 'L2path'])

    return idsa


# Get all svs images with its label as one file; level is the tile resolution level
def big_image_sum(pmd, path='../tiles/', dict_file='../tcia_pathology_slides.tsv',
                  ref_file='../gbm_all_subtype_collections.2019-10-07.tsv'):
    refdict = {'low': 0, 'high': 1, False: 0, True: 1, 'normal': 0, 'short': 1, 'long': 2}
    dct = pd.read_csv(dict_file, sep='\t', header=0)
    dct = dct.loc[dct['used_in_proteome'] == True]
    ref = pd.read_csv(ref_file, sep='\t', header=0)
    ref = ref.dropna(subset=[pmd])
    ref[pmd] = ref[pmd].replace(refdict)
    big_images = []
    if pmd == 'telomere':
        normalimg = intersection(ref.loc[ref[pmd] == 0]['case'].tolist(), dct['case_id'].tolist())
        normalsld = dct[dct['case_id'].isin(normalimg)]['slide_id'].tolist()
        shortimg = intersection(ref.loc[ref[pmd] == 1]['case'].tolist(), dct['case_id'].tolist())
        shortsld = dct[dct['case_id'].isin(shortimg)]['slide_id'].tolist()
        longimg = intersection(ref.loc[ref[pmd] == 2]['case'].tolist(), dct['case_id'].tolist())
        longsld = dct[dct['case_id'].isin(longimg)]['slide_id'].tolist()
        for i in normalsld:
            sldnum = i.split('-')[-1]
            pctnum = i[:-3]
            big_images.append([pctnum, 0, path + "{}/".format(pctnum), sldnum])
        for i in shortsld:
            sldnum = i.split('-')[-1]
            pctnum = i[:-3]
            big_images.append([pctnum, 1, path + "{}/".format(pctnum), sldnum])
        for i in longsld:
            sldnum = i.split('-')[-1]
            pctnum = i[:-3]
            big_images.append([pctnum, 2, path + "{}/".format(pctnum), sldnum])
    else:
        negimg = intersection(ref.loc[ref[pmd] == 0]['case'].tolist(), dct['case_id'].tolist())
        negsld = dct[dct['case_id'].isin(negimg)]['slide_id'].tolist()
        posimg = intersection(ref.loc[ref[pmd] == 1]['case'].tolist(), dct['case_id'].tolist())
        possld = dct[dct['case_id'].isin(posimg)]['slide_id'].tolist()
        for i in negsld:
            sldnum = i.split('-')[-1]
            pctnum = i[:-3]
            big_images.append([pctnum, 0, path + "{}/".format(pctnum), sldnum])
        for i in possld:
            sldnum = i.split('-')[-1]
            pctnum = i[:-3]
            big_images.append([pctnum, 1, path + "{}/".format(pctnum), sldnum])

    datapd = pd.DataFrame(big_images, columns=['slide', 'label', 'path', 'sldnum'])

    return datapd


# seperate into training and testing; each type is the same separation ratio on big images
# test and train csv files contain tiles' path.
def set_sep(alll, path, cls, cut=0.3, batchsize=24):
    trlist = []
    telist = []
    valist = []
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

    test_tiles = pd.DataFrame(columns=['slide', 'label', 'L0path', 'L1path', 'L2path'])
    train_tiles = pd.DataFrame(columns=['slide', 'label', 'L0path', 'L1path', 'L2path'])
    validation_tiles = pd.DataFrame(columns=['slide', 'label', 'L0path', 'L1path', 'L2path'])
    for idx, row in test.iterrows():
        tile_ids = paired_tile_ids_in(row['slide'], row['label'], row['path'], row['sldnum'])
        test_tiles = pd.concat([test_tiles, tile_ids])
    for idx, row in train.iterrows():
        tile_ids = paired_tile_ids_in(row['slide'], row['label'], row['path'], row['sldnum'])
        train_tiles = pd.concat([train_tiles, tile_ids])
    for idx, row in validation.iterrows():
        tile_ids = paired_tile_ids_in(row['slide'], row['label'], row['path'], row['sldnum'])
        validation_tiles = pd.concat([validation_tiles, tile_ids])

    # No shuffle on test set
    train_tiles = sku.shuffle(train_tiles)
    validation_tiles = sku.shuffle(validation_tiles)
    if train_tiles.shape[0] > int(batchsize*80000/3):
        train_tiles = train_tiles.sample(int(batchsize*80000/3), replace=False)
        print('Truncate training set!')
    if validation_tiles.shape[0] > int(batchsize*80000/30):
        validation_tiles = validation_tiles.sample(int(batchsize*80000/30), replace=False)
        print('Truncate validation set!')
    if test_tiles.shape[0] > int(batchsize*80000/3):
        test_tiles = test_tiles.sample(int(batchsize*80000/3), replace=False)
        print('Truncate test set!')

    test_tiles.to_csv(path+'/te_sample.csv', header=True, index=False)
    train_tiles.to_csv(path+'/tr_sample.csv', header=True, index=False)
    validation_tiles.to_csv(path+'/va_sample.csv', header=True, index=False)

    return train_tiles, test_tiles, validation_tiles
