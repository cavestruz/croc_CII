#from derived_field_CII import * # Contains yt and numpy
import matplotlib.pyplot as plt
import yt as yt
import numpy as np
from astropy.table import Table
from mpl_toolkits import mplot3d

ds = yt.load("/home/rnoorali/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")


all_data = ds.all_data()

ax = plt.axes(projection='3d')
ax.scatter3D(all_data['POSITION_X'],all_data['POSITION_Y'],all_data['POSITION_Z'])

plt.savefig('test_fig.png')
#ax.plot_trisurf(all_data['POSITION_X'],all_data['POSITION_Y'],all_data['POSITION_Z'],cmap='Pastels')
#plt.savefig('test_fig_tri.png')
