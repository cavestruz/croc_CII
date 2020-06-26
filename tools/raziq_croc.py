import yt as yt
import numpy as np
import matplotlib.pyplot as plt
from derived_field_CII import *


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
            slc.save("../SlicePlots/"+ p + '_sliceplot.png')
        else: 
            slc = yt.SlicePlot(ds, 'z',p)
            slc.save("../SlicePlots/"+ p + '_sliceplot.png')

def make_histogram_slice(dataset,list_or_array):
    all_data_at_z_0 = dataset.r[:,:,0]
    for p in list_or_array:
        minimum = np.amin(all_data_at_z_0[p])
        maximum = np.amax(all_data_at_z_0[p])
        ratio = abs(maximum/minimum)
        if ratio < 100:
            x,bin_edges = np.histogram(all_data_at_z_0[p],bins = 100) 
            plt.bar(bin_edges[:-1], x,width = bin_edges[1]-bin_edges[0])
            plt.title(p)
            plt.xlabel(p)
            plt.ylabel('Frequency')
            plt.savefig("../Histograms/"+ p + '_histogram.png')
            plt.clf()

        else:
            neg_data = all_data_at_z_0[p][all_data_at_z_0[p]<0]

            non_neg_data = all_data_at_z_0[p][all_data_at_z_0[p]>=0]


            if len(neg_data)>0:
                neg_min = np.amin(neg_data)
                neg_max = np.amax(neg_data)
                x,bin_edges = np.histogram(np.log10(np.absolute(neg_data)),bins = np.arange(np.around(np.log10(neg_min))-1,np.around(np.log10(neg_max))+1)) 
                plt.bar(bin_edges[:-1], x, width = 1) 
                plt.xlim(min(bin_edges), max(bin_edges))
                plt.title(p + 'Negative Values')
                plt.xlabel(p + " order of magnitude")
                plt.ylabel('Frequency')
                plt.savefig("../Histograms/"+ p + '_neg_histogram.png')
                plt.clf()

            if len(non_neg_data)>0:
                non_neg_min = np.amin(non_neg_data)
                non_neg_max = np.amax(non_neg_data)
                x,bin_edges = np.histogram(np.log10(np.absolute(non_neg_data)),bins = np.arange(np.around(np.log10(non_neg_min))-1,np.around(np.log10(non_neg_max))+1)) 
                plt.bar(bin_edges[:-1], x, width = 1) 
                plt.xlim(min(bin_edges), max(bin_edges))
                plt.title(p+ 'Positive Values')
                plt.xlabel(p + " order of magnitude")
                plt.ylabel('Frequency')
                plt.savefig("../Histograms/"+ p + '_pos_histogram.png')
                plt.clf()


ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rC2e','rC2a','C2_e_cooling','C2_a_cooling', 'C2_HeI_cooling']


#make_slice_plot(ds,plot_list)

make_histogram_slice(ds,plot_list)
