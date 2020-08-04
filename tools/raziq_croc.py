from derived_field_CII import * # Contains yt and numpy
import matplotlib.pyplot as plt
from astropy.table import Table

def make_slice_plot(dataset,list_or_array,full_box=False):
    for p in list_or_array:
        minimum = np.amin(dataset[p])
        maximum = np.amax(dataset[p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio))
        if ratio < 100:
            if full_box == True:
                slc = yt.SlicePlot(ds, 'z',p)
            elif full_box == False:
                slc = yt.SlicePlot(ds, 'z',p,data_source = dataset)
            else:
                raise ValueError("full_box must be either True or False")
           
            slc.set_log(p, False)
        else: 
            if full_box == True:
                slc = yt.SlicePlot(ds, 'z',p)
            elif full_box == False:
                slc = yt.SlicePlot(ds, 'z',p,data_source = dataset)
            else:
                raise ValueError("full_box must be either True or False")
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

#all_data_at_z_0 = ds.r[:,:,0]

#plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rCIIe','rCIIa','CII_e_cooling','CII_a_cooling', 'CII_HeI_cooling', 'CII_CMB_emission','CII_H2_ortho', 'CII_H2_para']
plot_list = ['CII_H2_ortho', 'CII_H2_para']

#make_slice_plot(all_data_at_z_0,plot_list)

#make_histogram_slice(all_data_at_z_0,plot_list)

halo_table = Table.read('/home/rnoorali/Data/halo_catalogs/out_14.list',format = "ascii.commented_header")

halo_table.add_index('Mvir')

largest_mass = halo_table[halo_table.loc_indices[np.amax(halo_table['Mvir'])]]

a = 1/(1+ds.current_redshift)
h = ds.hubble_constant
x = largest_mass['X']*a/h
y = largest_mass['Y']*a/h
z = largest_mass['Z']*a/h
r = largest_mass['Rvir']*a/h
sphere = ds.sphere([x, y, z],(2*r, "kpc"))

make_slice_plot(sphere,plot_list)
make_histogram_slice(sphere,plot_list)
