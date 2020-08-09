from derived_field_CII import * # Contains yt and numpy
import matplotlib.pyplot as plt
from astropy.table import Table

def make_slice_plot(dataset,list_or_array,full_box=False):
    for p in list_or_array:
        minimum = np.amin(dataset['data_object'][p])
        maximum = np.amax(dataset['data_object'][p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio))
        if ratio < 100:
            slc = yt.SlicePlot(ds, 'z',p,data_source = dataset['data_object'],center = dataset['data_object'].center,width = dataset['width'])
            slc.set_log(p, False)
        else: 
            slc = yt.SlicePlot(ds, 'z',p,data_source = dataset['data_object'],center = dataset['data_object'].center,width = dataset['width'])

        slc.save("../SlicePlots/"+ p + '_sliceplot.png')

def make_histogram_slice(dataset,list_or_array,log_frequency=False):

    for p in list_or_array:
        minimum = np.amin(dataset['data_object'][p])
        maximum = np.amax(dataset['data_object'][p])
        ratio = abs(maximum/minimum)
        if ratio < 100:
            x,bin_edges = np.histogram(dataset['data_object'][p],bins = 100) 
            plt.bar(bin_edges[:-1], x,width = bin_edges[1]-bin_edges[0])
            plt.title(p)
            plt.xlabel(p+' ('+str(dataset['data_object'][p].units)+')')
            plt.ylabel('Frequency')
            if log_frequency == True:
                plt.yscale('log')

        else:

            x,bin_edges = np.histogram(np.log10(dataset['data_object'][p]+0j),bins = np.arange(np.around(np.log10(minimum+0j))-1,np.around(np.log10(maximum+0j))+1)) 
            plt.bar(bin_edges[:-1], x, width = 1) 
            plt.xlim(min(bin_edges), max(bin_edges))
            plt.title(p )
            plt.xlabel('log10('+p+') ('+str(dataset['data_object'][p].units)+')')
            plt.ylabel('Frequency')
            if log_frequency == True:
                plt.yscale('log')
        plt.savefig("../Histograms/"+ p + '_histogram.png')
        plt.clf()
def make_sphere_region(row_of_table):
    a = 1/(1+ds.current_redshift)
    h = ds.hubble_constant
    x = row_of_table['X']*a/h
    y = row_of_table['Y']*a/h
    z = row_of_table['Z']*a/h
    r = row_of_table['Rvir']*a/h
    sphere = ds.sphere(yt.YTArray([x, y, z],"Mpc"),(2*r, "kpc"))
    return {'data_object':sphere,'width':2*sphere.radius}

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

#entire_box = {'data_object':ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge),'width':128}#Equivilently all_data()

#all_data_at_z_0 = ds.r[:,:,0]

#plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rCIIe','rCIIa','CII_e_cooling','CII_a_cooling', 'CII_HeI_cooling', 'CII_CMB_emission','CII_H2_ortho', 'CII_H2_para']
plot_list = ['temperature','metallicity']

halo_table = Table.read('/home/rnoorali/Data/halo_catalogs/out_14.list',format = "ascii.commented_header")

halo_table.add_index('Mvir')

largest_mass = halo_table[halo_table.loc_indices[np.amax(halo_table['Mvir'])]]

sphere = make_sphere_region(largest_mass)

make_slice_plot(sphere,plot_list)
make_histogram_slice(sphere,plot_list)
