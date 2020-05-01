import yt as yt
from yt.units import gram, centimeter
import numpy as np

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

amu_to_gram = 1.660539040*10**(-24)*gram # atomic mass unit to grams

L_soblen = 100. #Parsec

parsec_to_cm = 1./(3.2407792700054*10**(-19)) * centimeter

dust_cross_section = 2*10**(-21)



def _HI_number_density(field, data):
    return data["HI density"]/(1.007825032*amu_to_gram)
ds.add_field(("gas", "HI_number_density"), function=_HI_number_density, units="1/cm**3")

def _HII_number_density(field, data):
    return data["HII density"]/(2.014101778*amu_to_gram)
ds.add_field(("gas", "HII_number_density"), function=_HII_number_density, units="1/cm**3")

def _dust_attenuation(field, data):
    return np.exp(-(data["HI_number_density"]+2*data["HII_number_density"])*data['metallicity']*L_soblen*parsec_to_cm*dust_cross_section)
ds.add_field(("gas", "dust_attenuation"), function=_dust_attenuation, units="")


slc = yt.SlicePlot(ds, 'z','metallicity')
slc.save('metallicity.png')
slc = yt.SlicePlot(ds, 'z','HI_number_density')
slc.save('HI.png')
slc = yt.SlicePlot(ds, 'z','HII_number_density')
slc.save('HII.png')
slc = yt.SlicePlot(ds, 'z','dust_attenuation')
slc.save('dust.png')
