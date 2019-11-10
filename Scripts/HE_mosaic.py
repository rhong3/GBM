#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make tSNE mosaic

Created on Wed Jul 17 10:16:21 2019

@author: lwk, RH
"""

import pandas as pd
from PIL import Image
import sys
import os

filename = sys.argv[1]   # DEFAULT='tSNE_P_N'
bin = int(sys.argv[2])   # DEFAULT=50
size = int(sys.argv[3])  # DEFAULT=299
pdmd = sys.argv[4]
dirr = sys.argv[5]
outim = sys.argv[6]


# random select representative images and output the file paths
def sample(dat, md, bins):
    if md == 'telomere':
        classes = 3
        redict = {0: 'normal_score', 1: 'short_score', 2: 'long_score'}
    elif md == 'immune':
        redict = {0: 'low_score', 1: 'high_score'}
        classes = 2
    else:
        redict = {0: 'NEG_score', 1: 'POS_score'}
        classes = 2
    sampledls = []
    for m in range(classes):
        for i in range(bins):
            for j in range(bins):
                try:
                    sub = dat.loc[(dat['x_int'] == i) & (dat['y_int'] == j)
                                    & (dat[redict[m]] > 0.51) & (dat['True_label'] == m)]
                    picked = sub.sample(1, replace=False)
                    for idx, row in picked.iterrows():
                        sampledls.append([row['path'], row['x_int'], row['y_int']])
                except ValueError:
                    pass
    samples = pd.DataFrame(sampledls, columns=['impath', 'x_int', 'y_int'])
    return samples


if __name__ == "__main__":
    dirls = dirr.split(',')

    for i in dirls:
        try:
            ipdat = pd.read_csv('../Results/{}/out/{}.csv'.format(i, filename))
            imdat = sample(ipdat, pdmd, bin)
            imdat.to_csv('../Results/{}/out/tsne_selected.csv'.format(i), index=False)

            new_im = Image.new(mode='RGB', size=(size*(bin+1), size*(bin+1)), color='white')

            for rows in imdat.itertuples():
                impath = rows.impath
                x = rows.x_int
                y = rows.y_int
                try:
                    im = Image.open(impath)
                    im.thumbnail((size, size))
                    new_im.paste(im, ((x-1)*size, (bin-y)*size))
                except FileNotFoundError:
                    print(impath)
                    pass
            new_im.save(os.path.abspath('../Results/{}/out/{}.jpeg'.format(i, outim)), "JPEG")
            print('{} done'.format(i))
        except FileNotFoundError:
            print('{} passed'.format(i))
            pass


