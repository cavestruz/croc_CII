from derived_field_CII import * # Contains yt, unyt_array and numpy
import matplotlib.pyplot as plt
from astropy.table import Table

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
ds.add_particle_filter("young_stars")

col_list = ['CII_e_cooling','CII_a_cooling','CII_HeI_cooling','CII_CMB_emission','CII_H2_ortho','CII_H2_para','LCII_total','sphere_vol','SFR']
multiplier = 30991
mult_range = (2,3)
which_file = 'LCII_'+str(multiplier*mult_range[0])+'_'+str(multiplier*mult_range[1])+'.txt'
#try:
halo_table = Table.read('/home/rnoorali/Data/rnoorali_turbo/' + which_file,format = "ascii.commented_header")
#print(halo_table.colnames)
#except FileNotFoundError:
#    print('Error!')
# halo_table = Table.read('/home/rnoorali/Data/halo_catalogs/out_14.list',format = "ascii.commented_header")
# for table_col in col_list:
#     halo_table[table_col] = -0.0123

for ID in halo_table['ID']:
    row = ID-multiplier*mult_range[0]
#     if halo_table['LCII_total'][row]>0:
#         print('Row',row,'is already filled')
#         continue

#    print("Now analyzing row " +str(row))
    sphere = make_sphere_region(halo_table[row])
    data = sphere['data_object']
    
    for quant in col_list[:-2]:
        halo_table[quant][row] = float(data.quantities.weighted_average_quantity(('gas',quant),('gas','cell_volume')))
    
    halo_table['sphere_vol'][row] = data.quantities.total_quantity(('gas','cell_volume'))
    
#     radius = sphere['width'].to('cm')/2.0
    halo_table['sphere_vol'][row] = data.quantities.total_quantity(('gas','cell_volume'))
    total_young_star_mass = data.quantities.total_quantity(('young_stars','particle_mass'))
    halo_table['SFR'][row] = total_young_star_mass.to('Msun')/20/10**6
#     print(str(row))
    
    if row%100 ==0:
        halo_table.write('/home/rnoorali/Data/rnoorali_turbo/' + which_file,format = "ascii.commented_header",overwrite=True)
        del sphere,data#,radius
        print('Row '+ str(row)+' Writing and clearing Cache')
#     if row>=75:
#         print('Test Completed')
#         break
halo_table.write('/home/rnoorali/Data/rnoorali_turbo/' + which_file,format = "ascii.commented_header",overwrite=True)



