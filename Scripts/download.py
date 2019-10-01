"""
Download SVS slides from TCIA

Created on 10/01/2019

@author: RH
"""
import urllib.request
import pandas as pd

dictkey = pd.read_csv('../tcia_pathology_slides.tsv', sep='\t', header=0)

for idx, row in dictkey.iterrows():
    urllib.request.urlretrieve(row['slide_url'], '../images/{}.svs'.format(row['slide_id']))

print("done!")
