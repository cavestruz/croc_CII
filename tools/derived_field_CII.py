import yt as yt
from yt.units import Kelvin, gram, kboltz, erg, centimeter
import numpy as np

#CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL CONSTANTS: 2014
amu = 1.660539040*10**(-24)*gram # atomic mass unit to kilograms

parsec = 3.085678 * 10**(18) * centimeter # parsec to meters

electron_mass = 5.48579909070*10**(-4) * amu

L_soblen = 100. * parsec

dust_cross_section = 2*10**(-21) * centimeter * centimeter

CII_abun = 3.31*10**(-4)


### All mass values come from the CRC Handbook of Chemistry and Physics May 2020


def _HI_number_density(field, data):
    return data["HI density"]/(1.007825032*amu)
yt.add_field(("gas", "HI number density"), function=_HI_number_density, units="1/cm**3")

def _HII_number_density(field, data):
    return data["HII density"]/(1.007825032*amu-electron_mass)
yt.add_field(("gas", "HII number density"), function=_HII_number_density, units="1/cm**3")

def _H2_number_density(field, data):
    return data["H2 density"]/(2.016*amu)
yt.add_field(("gas", "H2 number density"), function=_H2_number_density, units="1/cm**3")

def _HeI_number_density(field, data):
    return data["HeI density"]/(4.002602*amu)
yt.add_field(("gas", "HeI number density"), function=_HeI_number_density, units="1/cm**3")

def _HeII_number_density(field, data):
    return data["HeII density"]/(4.00260*amu-electron_mass)
yt.add_field(("gas", "HeII number density"), function=_HeII_number_density, units="1/cm**3")

def _HeIII_number_density(field, data):
    return data["HeIII density"]/(4.00260*amu-2*electron_mass)
yt.add_field(("gas", "HeIII number density"), function=_HeIII_number_density, units="1/cm**3")

def _log_dust_attenuation(field, data):
    return -(data["HI number density"]+2*data["H2 number density"])*data['metallicity']*L_soblen*dust_cross_section # With H2
#    return -(data["HI number density"])*data['metallicity']*L_soblen*dust_cross_section # Without H2
yt.add_field(("gas", "log_dust_attenuation"), function=_log_dust_attenuation, units="")

def _rCIIe(field, data):
    return 6.67*10**(-20)*Kelvin**(0.5)*np.exp(-91.2*Kelvin/data['temperature'])/data["temperature"]**0.5
yt.add_field(("gas", "rCIIe"), function=_rCIIe, units="")

def _rCIIa(field, data):
    x = 16 + .344*(data['temperature']/Kelvin)**(0.5) - 47.7*Kelvin/data['temperature']
    x[x<0.0]=0
    return 10**(-24)*np.exp(-91.2*Kelvin/data['temperature'])*x
yt.add_field(("gas", "rCIIa"), function=_rCIIe, units="")

def _CII_e_cooling(field, data):
    return data['rCIIe'] * (data['HII number density']+data['HeII number density']+2*data['HeIII number density']) * CII_abun * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_e_cooling"), function=_CII_e_cooling, units="1/cm**6")

def _CII_a_cooling(field, data):
    return data['rCIIa'] * CII_abun * data['HI number density'] * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_a_cooling"), function=_CII_a_cooling, units="1/cm**6")

def _CII_HeI_cooling(field, data):
    return 0.38 * data['rCIIa'] * data['HeI number density'] * CII_abun * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_HeI_cooling"), function=_CII_HeI_cooling, units="1/cm**6")

def _CII_CMB_emission(field, data):
    z = data.ds.current_redshift
    return 2.0 *  CII_abun * data['HI number density'] * data['metallicity']*2.298*10**(-6)*kboltz*91.2*Kelvin*np.exp(-91.2*Kelvin/(2.725*Kelvin)*1/(1+z))
yt.add_field(("gas", "CII_CMB_emission"), function=_CII_CMB_emission, units="erg/cm**3")

def _CII_H2_para(field, data):
    rate_coefficient = 4.25*10**(-10)*erg*centimeter**3/second
    return rate_coefficient*pow(data['temperature']/(100*Kelvin),0.124-0.018*np.log(data['temperature']/(100*Kelvin)))*kboltz*(91.2*Kelvin)*CII_abun*data['HI number density']*0.25*data['H2 number density']
yt.add_field(("gas", "CII_H2_para"), function=_CII_H2_para, units="erg**2/(cm**3*second)") # Draine Table F6 pg 501

def _CII_H2_ortho(field, data):
    return 5.14*10**(-10)*pow(data['temperature']/(100*Kelvin),0.124-0.018*np.log(data['temperature']/(100*Kelvin)))*kboltz*91.2*Kelvin*CII_abun*data['HI number density']*0.75*data['H2 number density'] /erg
yt.add_field(("gas", "CII_H2_ortho"), function=_CII_H2_ortho, units="1/cm**6")

#slc = yt.SlicePlot(ds, 'z','rC2a')
#slc.save('rC2a.png')
#slc = yt.SlicePlot(ds, 'z','HI_number_density')
#slc.save('HI.png')
#slc = yt.SlicePlot(ds, 'z','HII_number_density')
#slc.save('HII.png')
#slc = yt.SlicePlot(ds, 'z','log_dust_attenuation')
#slc.set_log('log_dust_attenuation',False)
#slc.set_zlim('log_dust_attenuation', -0.02, 0)
#slc.save('dust.png')


#sl = ds.r[:,:,0.]

#data = sl['log_dust_attenuation']

#def count_elements(seq) -> dict:
#     """Tally elements from `seq`."""
#     hist = {}
#     for i in seq:
#         hist[i] = hist.get(i, 0) + 1
#     return hist
#dictionary = count_elements(data)

#plt.bar(list(dictionary.keys()), dictionary.values(), color='g')

#plt.savefig('dust_hist.png')
