"""
Tile svs/scn files; modified for GBM

Created on 10/02/2019

@author: RH
"""

import time
import matplotlib
import os
import shutil
import pandas as pd
import numpy as np
matplotlib.use('Agg')
import Slicer
import staintools


# Get all images in the root directory
def image_ids_in(root_dir, ignore=['.DS_Store', 'dict.csv']):
    ids = []
    for id in os.listdir(root_dir):
        if id in ignore:
            print('Skipping ID:', id)
        else:
            imgname = id.split('.')[0]
            dirname = imgname.split('-2')[0]
            sldnum = imgname.split('-')[-1]
            ids.append((id, dirname, sldnum))
    return ids


# cut; each level is 2 times difference (20x, 10x, 5x)
def cut():
    # load standard image for normalization
    std = staintools.read_image("../colorstandard.png")
    std = staintools.LuminosityStandardizer.standardize(std)
    CPTACpath = '../images/'
    ref = pd.read_csv('../tcia_pathology_slides.tsv', sep='\t', header=0)

    # cut tiles with coordinates in the name (exclude white)
    start_time = time.time()
    CPTAClist = image_ids_in(CPTACpath)

    # CPTAC
    for i in CPTAClist:
        matchrow = ref.loc[ref['case_id'] == i[1]]
        if matchrow.empty:
            continue
        try:
            os.mkdir("../tiles/{}".format(i[1]))
        except(FileExistsError):
            pass
        for m in range(3):
            if m == 0:
                tff = 1
                level = 0
            elif m == 1:
                tff = 2
                level = 0
            else:
                tff = 1
                level = 1
            otdir = "../tiles/{}/level{}".format(i[1], str(m))
            try:
                os.mkdir(otdir)
            except(FileExistsError):
                continue
            try:
                n_x, n_y, raw_img, resx, resy, ct = Slicer.tile(image_file=i[0], outdir=otdir,
                                                                level=level, std_img=std, dp=i[2], ft=tff)
            except(IndexError):
                pass
            if len(os.listdir(otdir)) < 2:
                shutil.rmtree(otdir, ignore_errors=True)

    print("--- %s seconds ---" % (time.time() - start_time))

    # # Time measure tool
    # start_time = time.time()
    # print("--- %s seconds ---" % (time.time() - start_time))


# Run as main
if __name__ == "__main__":
    if not os.path.isdir('../tiles'):
        os.mkdir('../tiles')
    cut()

