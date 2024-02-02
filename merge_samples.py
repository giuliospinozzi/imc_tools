import os
import pandas as pd
import glob
import argparse
import shutil
from pathlib import Path

def copy_all(path,newpath):
    files = os.listdir(path)
    for f in files:
        print(f)
        shutil.copy(os.path.join(path,f),newpath)

parser = argparse.ArgumentParser()

parser.add_argument("in_path", help="Path of the steinbock data")
parser.add_argument("out_path", help="Path to save merged data")
parser.add_argument("--has_neighbors", help="If selected, copies the neighbourhood analysis too", action='store_true')
parser.add_argument("--no_copy", help="skip the actual copy", action='store_true')

args = parser.parse_args()

in_path = args.in_path
out_path = args.out_path
nocopy = args.no_copy
has_neighbors = args.has_neighbors

if nocopy:
    do_copy = False
else:
    do_copy = True

Path(os.path.join(out_path,'img')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(out_path,'masks')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(out_path,'raw')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(out_path,'intensities')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(out_path,'regionprops')).mkdir(parents=True, exist_ok=True)

if has_neighbors:
    Path(os.path.join(out_path,'neighbors')).mkdir(parents=True, exist_ok=True)

dirs = os.listdir(in_path)
good_dirs = []

for d in dirs:
    if os.path.isdir(os.path.join(in_path,d)):
        good_dirs.append(os.path.join(in_path,d))

if do_copy:
    for d in good_dirs:
        print(f'Copying files from sample {d}')
        print('img')
        copy_all(os.path.join(d,'img'),os.path.join(out_path,'img/'))
        print('masks')
        copy_all(os.path.join(d,'masks'),os.path.join(out_path,'masks/'))
        print('intensities')
        copy_all(os.path.join(d,'intensities'),os.path.join(out_path,'intensities/'))
        print('raw')
        copy_all(os.path.join(d,'raw'),os.path.join(out_path,'raw/'))
        print('regionprops')
        copy_all(os.path.join(d,'regionprops'),os.path.join(out_path,'regionprops/'))
        print('panel.csv')
        shutil.copy(os.path.join(d,'panel.csv'),os.path.join(out_path,'panel.csv'))
        if has_neighbors:
            print('neighbors')
            copy_all(os.path.join(d,'neighbors'),os.path.join(out_path,'neighbors/'))
        

print('Merging images.csv')
    
files = glob.glob(in_path + '/*/images.csv')

images = []

for f in files:
    images.append(pd.read_csv(f))

all_images = pd.concat(images)
all_images.to_csv(out_path + '/images.csv',index=None)
