#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Created By  : Giulio Spinozzi, PhD
# Created Date: 31/01/2024
# version ='1.0'
# ---------------------------------------------------------------------------

import argparse
import os
import cv2
from scipy import signal
import tifffile as tf
import numpy as np
import shutil

def process_array(img, filtertype, filtervalue):
    if filtertype == 'tophat':
        filterSize =(filtervalue, filtervalue)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT,filterSize)
        tophat_data = cv2.morphologyEx(img,cv2.MORPH_TOPHAT,kernel)
        processed_img = img - tophat_data
    elif filtertype=='median':
        processed_img = signal.medfilt2d(img, filtervalue)

    return processed_img

parser = argparse.ArgumentParser()

parser.add_argument('steinbock_path', help='Path to the steinbock sample folder')
parser.add_argument('--type', help='Type of filtering to perform',required=True)
parser.add_argument('--value', help='Filtering kernel value',required=True,type=int)
parser.add_argument('--channels', help='List of channels to process separated by commas, if omitted processes all channels')

args = parser.parse_args()

rootpath = args.steinbock_path
filtertype = args.type
filtervalue = args.value
channels = args.channels

inpath = os.path.join(rootpath,'img')
renamed_path = os.path.join(rootpath,'img_old')

os.rename(inpath,renamed_path)
os.makedirs(inpath,exist_ok=True)

print('[OPBG] PYTHON HOT PIXEL FILTER -> START')
for f in os.listdir(renamed_path):
    if f.endswith('.tiff'):
        filename = os.path.join(renamed_path,f)
        print(f'processing {f}')
        img = tf.imread(filename)
        #create new image
        img_new = np.zeros_like(img,dtype='uint16')
        
        for c in range(img.shape[0]):
            if channels:
                if c in channels:
                    data = img[c,:,:]
                    data = data / (2**32)
                    data = data * (2**16)
                    data.dtype = 'uint16'
                    data_processed = process_array(data, filtertype, filtervalue)
                    img_new[c,:,:] = data_processed
                else:
                    data = img[c,:,:]
                    data = data / (2**32)
                    data = data * (2**16)
                    data.dtype = 'uint16'
                    img_new[c,:,:] = data
            else:
                data = img[c,:,:]
                data_processed = process_array(data, filtertype, filtervalue)
                img_new[c,:,:] = data_processed
        
        newfilename = os.path.join(inpath,f)
        tf.imwrite(newfilename,img_new)
        print('Done')
print('Done!')

print('[OPBG] PYTHON HOT PIXEL FILTER -> Renaming Folder in External')
external_path = os.path.join(rootpath,'external')
os.rename(inpath,external_path)
print('Done!')

print('[OPBG] PYTHON HOT PIXEL FILTER -> Compressing original img folder')
shutil.make_archive(renamed_path, 'zip', root_dir=renamed_path)
shutil.rmtree(renamed_path)
print('[OPBG] PYTHON HOT PIXEL FILTER -> END')
