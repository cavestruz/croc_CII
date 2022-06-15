from derived_field_CII import * # Contains yt, unyt_array and numpy
import matplotlib.pyplot as plt
from astropy.table import Table


def make_slice_plot(dataset,list_or_array):
    """
    Takes a given dataset and width and makes a slice plot along the z axis. 
    Includes infrastructure to plot log based or linear based on overall ratio 
    of highest to lowest values. 
    Future development includes including a way to specify which axis to slice
    on and how to appropriately determine when to take the log

    **Parameters**

    :dataset: dictionary
    
         dictionary in the form of {'data_object':<yt.data_object>,'width':<float in kpc>}
    
    :list_or_array: iterable python data structure, e.g. list, array, etc.. 
         
         This represents the fields in yt in which we want to make sliceplots en
         masse. Choosing a field name not in yt will return an error. Future
         development includes adding a config file to choose fields manually
         

    """
    for p in list_or_array:
        minimum = np.amin(dataset['data_object'][p])
        maximum = np.amax(dataset['data_object'][p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio)) # For testing of log vs linear
        if ratio < 100:
            slc = yt.SlicePlot(ds, 'z',p,data_source = dataset['data_object'],center = dataset['data_object'].center,width = dataset['width'])
            slc.set_log(p, False)
        else: 
            slc = yt.SlicePlot(ds, 'z',p,data_source = dataset['data_object'],center = dataset['data_object'].center,width = dataset['width'])

        slc.save("../SlicePlots/"+ p + '_sliceplot.png')
    return 


def make_masked_slice_plot(dataset,list_or_array,mask_list_or_array):
    """
    Takes a given dataset and width and makes a slice plot along the z axis. 
    Includes infrastructure to plot log based or linear based on overall ratio 
    of highest to lowest values. 
    Future development includes including a way to specify which axis to slice
    on and how to appropriately determine when to take the log

    **Parameters**

    :dataset: dictionary
    
         dictionary in the form of {'data_object':<yt.data_object>,'width':<float in kpc>}
    
    :list_or_array: iterable python data structure, e.g. list, array, etc.. 
         
         This represents the fields in yt in which we want to make sliceplots en
         masse. Choosing a field name not in yt will return an error. Future
         development includes adding a config file to choose fields manually
         

    """
    if len(list_or_array) !=  len(mask_list_or_array):
        raise NameError('Check your array sizes')
    
    for p,mask in zip(list_or_array,mask_list_or_array):
#         print(type(dataset['data_object'][p]))
        masked_dataset = ds.cut_region(dataset['data_object'],["obj[('gas', '"+str(p)+"')] > "+str(mask)])
        minimum = np.amin(masked_dataset[p])
        maximum = np.amax(masked_dataset[p])
        ratio = maximum/minimum
#        print(p+"_"+str(ratio)) # For testing of log vs linear
        if ratio < 100:
            slc = yt.SlicePlot(ds, 'z',p,data_source = masked_dataset,center = dataset['data_object'].center,width = dataset['width'])
            slc.set_log(p, False)
        else: 
            slc = yt.SlicePlot(ds, 'z',p,data_source = masked_dataset,center = dataset['data_object'].center,width = dataset['width'])

        slc.save("../SlicePlots/"+ p + '_sliceplot.png')
        
def make_histogram_slice(dataset,list_or_array,log_frequency=False):

    """
    Takes a given dataset and makes a histogram of all values in the set. 
    Includes infrastructure to plot either log based or linear based values
    categorical values (ratio), as well as log based or linear based frequencies.  
    Future development includes including how to autmatically determine when to take the log of categories and log of frequencies. 

    **Parameters**

    :dataset: dictionary
    
         dictionary in the form of {'data_object':<yt.data_object>,'width':<float in kpc>}
    
    :list_or_array: iterable python data structure, e.g. list, array, etc.. 
         
         This represents the fields in yt in which we want to make sliceplots en
         masse. Choosing a field name not in yt will return an error. Future
         development includes adding a config file to choose fields manually

    :log_frequency: boolean, optional
    
         Sets the frequency axis to log10
         

    """

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
def make_sphere_region(row_of_table,zoom_factor = 1.0):
    """
    Takes a specified row in DM halo catalog and makes a spherical region file 
    around it. Auto converts from comoving coordinates. 

    **Parameters**

    :row_of_table: astropy.Table Table

         Astropy table containing at least the coordinates and the virial radius
         of the DM halo. 

    :zoom_factor: float, optional
    
         Used to create a region file zoomed in on the origin. The larger the
         factor, the more zoomed in the region. Defaults to creating a region
         twice the radius of the virial radius


    """
    a = 1/(1+ds.current_redshift)
    h = ds.hubble_constant
    x = row_of_table['X']*a/h
    y = row_of_table['Y']*a/h
    z = row_of_table['Z']*a/h
    r = row_of_table['Rvir']*a/h
#     print(a,h,x,y,z,r,2*r/zoom_factor)
    center_vector = yt.YTArray([x, y, z],units="Mpc")
    sphere = ds.sphere(center_vector,(2*r/zoom_factor, "kpc"))
    return {'data_object':sphere,'width':2*sphere.radius}

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

# entire_box = {'data_object':ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge),'width':128}#Equivilently all_data()

#all_data_at_z_0 = ds.r[:,:,0]

#plot_list = ['temperature','metallicity','HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rCIIe','rCIIa','CII_e_cooling','CII_a_cooling', 'CII_HeI_cooling', 'CII_CMB_emission','CII_H2_ortho', 'CII_H2_para']
plot_list = ['LCII_H2_ortho','LCII_H2_para','LCII_total']

halo_table = Table.read('/home/rnoorali/Data/halo_catalogs/out_14.list',format = "ascii.commented_header")

halo_table.add_index('Mvir')

largest_mass = halo_table[halo_table.loc_indices[np.amax(halo_table['Mvir'])]]

sphere = make_sphere_region(largest_mass,zoom_factor=14.0)

make_slice_plot(sphere,plot_list)
# make_masked_slice_plot(sphere,plot_list,[10**29,10**24,10**20,10**32])
# make_histogram_slice(sphere,plot_list)

# plot_list = ['log_dust_attenuation','rCIIe','rCIIa']
# make_histogram_slice(sphere,plot_list,log_frequency=True)

