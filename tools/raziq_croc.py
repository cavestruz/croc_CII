import yt as yt
import numpy as np
import matplotlib.pyplot as plt
from derived_field_CII import *

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

plot_list = ['HI number density','HII number density','HeI number density','HeII number density','HeIII number density','log_dust_attenuation','rC2e','rC2a','C2_e_cooling','C2_a_cooling', 'C2_HeI_cooling']

for p in plot_list:
    slc = yt.SlicePlot(ds, 'z',p)
    slc.save(p + '.png')
