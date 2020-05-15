import yt as yt
import numpy as np
import matplotlib.pyplot as plt
from derived_field_CII import *

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rC2e','rC2a','C2_e_cooling','C2_a_cooling', 'C2_HeI_cooling']

def make_slice_plot(dataset,list_or_array):
    all_data_at_z_0 = dataset.r[:,:,0]
    for p in list_or_array:
        minimum = np.amin(all_data_at_z_0[p])
        maximum = np.amax(all_data_at_z_0[p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio))
        if ratio < 100:
            slc = yt.SlicePlot(ds, 'z',p)
            slc.set_log(p, False)
            slc.save(p + '_sliceplot.png')
        else: 
            slc = yt.SlicePlot(ds, 'z',p)
            slc.save(p + '_sliceplot.png')

make_slice_plot(ds,plot_list)
