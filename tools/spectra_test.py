import yt 
import trident 
import numpy as np 
import sys

if len(sys.argv)<5:
    print('usage:python spectra_test.py datafile start_position{left_edge} end_position{right_edge} savename')

datafile = sys.argv[1]

start_position = sys.argv[2]

end_position = sys.argv[3]

savename = sys.argv[4]

line_list = ['C','N','O']

ds = yt.load(datafile) # Loading the dataset

if start_position == 'left_edge':
    start_position = ds.domain_left_edge
if end_position == 'right_edge':
    end_position = ds.domain_right_edge


trident.add_ion_fields(ds,ions=line_list) 

actual_ray = trident.make_simple_ray(ds,start_position = np.array([0.,0.,0.]),end_position = np.array([128.,128.,128.]),data_filename="mray.h5",lines=line_list,ftype='gas') 

#lr = trident.LightRay(ds)
#lr.make_light_ray(ds,start_position = start_position,end_position = end_position,data_filename="mray.h5",fields=['temperature','density'])

sg = trident.SpectrumGenerator('COS')
#sg = trident.SpectrumGenerator(lambda_min='auto', lambda_max='auto',dlambda=0.01)
sg.make_spectrum(actual_ray, lines=line_list) 
sg.save_spectrum(savename + '.txt') 
sg.plot_spectrum(savename + '.png') 
