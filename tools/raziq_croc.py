import yt as yt
import numpy as np
import matplotlib.pyplot as plt
from derived_field_CII import *


def make_slice_plot(dataset,list_or_array):
    for p in list_or_array:
        minimum = np.amin(dataset[p])
        maximum = np.amax(dataset[p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio))
        if ratio < 100:
            slc = yt.SlicePlot(ds, 'z',p)
            slc.set_log(p, False)
            slc.save("../SlicePlots/"+ p + '_sliceplot.png')
        else: 
            slc = yt.SlicePlot(ds, 'z',p)
            slc.save("../SlicePlots/"+ p + '_sliceplot.png')

def make_histogram_slice(dataset,list_or_array):

    for p in list_or_array:
        minimum = np.amin(dataset[p])
        maximum = np.amax(dataset[p])
        ratio = abs(maximum/minimum)
        if ratio < 100:
            x,bin_edges = np.histogram(dataset[p],bins = 100) 
            plt.bar(bin_edges[:-1], x,width = bin_edges[1]-bin_edges[0])
            plt.title(p)
            plt.xlabel(p)
            plt.ylabel('Frequency')
            plt.savefig("../Histograms/"+ p + '_histogram.png')
            plt.clf()

        else:

            x,bin_edges = np.histogram(np.log10(dataset[p]+0j),bins = np.arange(np.around(np.log10(minimum+0j))-1,np.around(np.log10(maximum+0j))+1)) 
            plt.bar(bin_edges[:-1], x, width = 1) 
            plt.xlim(min(bin_edges), max(bin_edges))
            plt.title(p )
            plt.xlabel(p + " order of magnitude")
            plt.ylabel('Frequency')
            plt.savefig("../Histograms/"+ p + '_histogram.png')
            plt.clf()



ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")
all_data_at_z_0 = ds.r[:,:,0]

plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rC2e','rC2a','C2_e_cooling','C2_a_cooling', 'C2_HeI_cooling']


make_slice_plot(all_data_at_z_0,plot_list)

make_histogram_slice(all_data_at_z_0,plot_list)
