import yt as yt
from yt.units import gram, centimeter
import numpy as np
import matplotlib.pyplot as plt


ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

amu_to_gram = 1.660539040*10**(-24)*gram # atomic mass unit to grams

L_soblen = 100. #Parsec

parsec_to_cm = 1./(3.2407792700054*10**(-19)) * centimeter

dust_cross_section = 2*10**(-21) * centimeter * centimeter



def _HI_number_density(field, data):
    return data["HI density"]/(1.007825032*amu_to_gram)
ds.add_field(("gas", "HI_number_density"), function=_HI_number_density, units="1/cm**3")

def _HII_number_density(field, data):
    return data["HII density"]/(2.014101778*amu_to_gram)
ds.add_field(("gas", "HII_number_density"), function=_HII_number_density, units="1/cm**3")

def _log_dust_attenuation(field, data):
    return -(data["HI_number_density"]+2*data["HII_number_density"])*data['metallicity']*L_soblen*parsec_to_cm*dust_cross_section
ds.add_field(("gas", "log_dust_attenuation"), function=_log_dust_attenuation, units="")


#slc = yt.SlicePlot(ds, 'z','metallicity')
#slc.save('metallicity.png')
#slc = yt.SlicePlot(ds, 'z','HI_number_density')
#slc.save('HI.png')
#slc = yt.SlicePlot(ds, 'z','HII_number_density')
#slc.save('HII.png')
#slc = yt.SlicePlot(ds, 'z','log_dust_attenuation')
#slc.set_log('log_dust_attenuation',False)
#slc.set_zlim('log_dust_attenuation', -0.02, 0)
#slc.save('dust.png')


sl = ds.r[:,:,0.]

data = sl['log_dust_attenuation']

def count_elements(seq) -> dict:
     """Tally elements from `seq`."""
     hist = {}
     for i in seq:
         hist[i] = hist.get(i, 0) + 1
     return hist
dictionary = count_elements(data)

plt.bar(list(dictionary.keys()), dictionary.values(), color='g')

plt.savefig('dust_hist.png')
