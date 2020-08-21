import yt as yt
from yt.units import Kelvin, gram, kboltz, erg, centimeter, second
import numpy as np

### CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL CONSTANTS: 2014
amu = 1.660539040*10**(-24)*gram # atomic mass unit to kilograms

parsec = 3.085678 * 10**(18) * centimeter # parsec to meters

electron_mass = 5.48579909070*10**(-4) * amu

### Specific units from c_emission.c

L_soblen = 100. * parsec # Not sure how this unit was derived

dust_cross_section = 2*10**(-21) * centimeter * centimeter # Not sure how this unit was derived

CII_abun = 3.31*10**(-4) # Not sure how this unit was derived


### All mass values come from the CRC Handbook of Chemistry and Physics May 2020

<<<<<<< HEAD
H_mass = 1.007825032*amu

H2_mass = 2.016*amu

He_mass = 4.002602*amu

### Conversions to YTEP-0003 compatible fields for use in Trident

def _H_p0_density_YTEP_0003(field,data):
    return data["HI density"].copy()
yt.add_field(("gas", "H_p0_density"), function=_H_p0_density_YTEP_0003, units="g/cm**3")

def _H_p1_density_YTEP_0003(field,data):
    return data["HII density"].copy()
yt.add_field(("gas", "H_p1_density"), function=_H_p1_density_YTEP_0003, units="g/cm**3")

def _H2_p0_density_YTEP_0003(field,data):
    return data["H2 density"].copy()
yt.add_field(("gas", "H2_p0_density"), function=_H2_p0_density_YTEP_0003, units="g/cm**3")

def _He_p0_density_YTEP_0003(field,data):
    return data["HeI density"].copy()
yt.add_field(("gas", "He_p0_density"), function=_He_p0_density_YTEP_0003, units="g/cm**3")

def _He_p1_density_YTEP_0003(field,data):
    return data["HeII density"].copy()
yt.add_field(("gas", "He_p1_density"), function=_He_p1_density_YTEP_0003, units="g/cm**3")

def _He_p2_density_YTEP_0003(field,data):
    return data["HeIII density"].copy()
yt.add_field(("gas", "He_p2_density"), function=_He_p2_density_YTEP_0003, units="g/cm**3")

### Functions from cII_emission.c
=======
def _test_mass(field, data):
    return data["HI density"]/data["HI density"]*(1.007825032*amu)
yt.add_field(("gas", "test_mass"), function=_test_mass, units="g")
>>>>>>> 2b8736913fff87feb8afd4fbaafc9482faad42da

def _HI_number_density(field, data):
                 return data["HI density"]/(H_mass)
yt.add_field(("gas", "HI number density"), function=_HI_number_density, units="1/cm**3")
yt.add_field(("gas", "H_p0_number_density"), function=_HI_number_density, units="1/cm**3") # Adds additional field for YTEP 0003 and Trident compatability

def _HII_number_density(field, data):
    return data["HII density"]/(H_mass-electron_mass)
yt.add_field(("gas", "HII number density"), function=_HII_number_density, units="1/cm**3")
yt.add_field(("gas", "H_p1_number_density"), function=_HII_number_density, units="1/cm**3")# Adds additional field for YTEP 0003 and Trident compatability

def _H2_number_density(field, data):
    return data["H2 density"]/(H2_mass)
yt.add_field(("gas", "H2 number density"), function=_H2_number_density, units="1/cm**3")
yt.add_field(("gas", "H2_p0_number_density"), function=_H2_number_density, units="1/cm**3")# Adds additional field for YTEP 0003 and Trident compatability

def _HeI_number_density(field, data):
    return data["HeI density"]/(He_mass)
yt.add_field(("gas", "HeI number density"), function=_HeI_number_density, units="1/cm**3")
yt.add_field(("gas", "He_p0_number_density"), function=_HeI_number_density, units="1/cm**3")# Adds additional field for YTEP 0003 and Trident compatability

def _HeII_number_density(field, data):
    return data["HeII density"]/(He_mass-electron_mass)
yt.add_field(("gas", "HeII number density"), function=_HeII_number_density, units="1/cm**3")
yt.add_field(("gas", "He_p1_number_density"), function=_HeII_number_density, units="1/cm**3")# Adds additional field for YTEP 0003 and Trident compatability

def _HeIII_number_density(field, data):
    return data["HeIII density"]/(He_mass-2*electron_mass)
yt.add_field(("gas", "HeIII number density"), function=_HeIII_number_density, units="1/cm**3")
yt.add_field(("gas", "He_p2_number_density"), function=_HeIII_number_density, units="1/cm**3")# Adds additional field for YTEP 0003 and Trident compatability

def _log_dust_attenuation(field, data):
    return -(data["HI number density"]+2*data["H2 number density"])*data['metallicity']*L_soblen*dust_cross_section # With H2
#    return -(data["HI number density"])*data['metallicity']*L_soblen*dust_cross_section # Without H2
yt.add_field(("gas", "log_dust_attenuation"), function=_log_dust_attenuation, units="")

def _rCIIe(field, data):
    return 6.67*10**(-20)*Kelvin**(0.5)*np.exp(-91.2*Kelvin/data['temperature'])/data["temperature"]**0.5*erg*centimeter**3/second #Not sure which one contains the per second and why there is a square root of temperature
yt.add_field(("gas", "rCIIe"), function=_rCIIe, units="erg*cm**3/s")

def _rCIIa(field, data):
    x = 16 + .344*(data['temperature']/Kelvin)**(0.5) - 47.7*Kelvin/data['temperature']
    x[x<0.0]=0
    return 10**(-24)*np.exp(-91.2*Kelvin/data['temperature'])*x*erg*centimeter**3/second #Not sure which one contains the per second and why there is a square root of temperature                
yt.add_field(("gas", "rCIIa"), function=_rCIIe, units="erg*cm**3/s")

def _CII_e_cooling(field, data):
    return data['rCIIe'] * (data['HII number density']+data['HeII number density']+2*data['HeIII number density']) * CII_abun * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_e_cooling"), function=_CII_e_cooling, units="erg/cm**3/s")

def _CII_a_cooling(field, data):
    return data['rCIIa'] * CII_abun * data['HI number density'] * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_a_cooling"), function=_CII_a_cooling, units="erg/cm**3/s")

def _CII_HeI_cooling(field, data):
    return 0.38 * data['rCIIa'] * data['HeI number density'] * CII_abun * data['HI number density'] * data['metallicity']
yt.add_field(("gas", "CII_HeI_cooling"), function=_CII_HeI_cooling, units="erg/cm**3/s")

def _CII_CMB_emission(field, data):
    z = data.ds.current_redshift
    return 2.0 *  CII_abun * data['HI number density'] * data['metallicity']*2.298*10**(-6)*kboltz*91.2*Kelvin*np.exp(-91.2*Kelvin/(2.725*Kelvin)*1/(1+z))/second # not sure which one contains the per second unit
yt.add_field(("gas", "CII_CMB_emission"), function=_CII_CMB_emission, units="erg/cm**3/s")

def _CII_H2_para(field, data):
    """                                                                                                               
    These are cooling rates, i.e. Lambda*n^2 where Lambda is the cooling function. See https://www.astro.umd.edu/~richard/ASTRO620/A620_2015_Gas_lec2.pdf                  
    """
    rate_coefficient = 4.25*10**(-10)*centimeter**3/second*pow(data['temperature']/(100*Kelvin),0.124-0.018*np.log(data['temperature']/(100*Kelvin))) # Draine Table F6 pg 501
    return rate_coefficient*kboltz*(91.2*Kelvin)*CII_abun*data['HI number density']*0.25*data['H2 number density']
yt.add_field(("gas", "CII_H2_para"), function=_CII_H2_para, units="erg/cm**3/s")

def _CII_H2_ortho(field, data):
    """                                                                                                               
    These are cooling rates, i.e. Lambda*n^2 where Lambda is the cooling function.. See https://www.astro.umd.edu/~richard/ASTRO620/A620_2015_Gas_lec2.pdf                 
    """
    rate_coefficient = (5.14*10**(-10)*centimeter**3/second)*pow(data['temperature']/(100*Kelvin),0.095+0.023*np.log(data['temperature']/(100*Kelvin)))# Draine Table F6 pg 501
    return rate_coefficient*kboltz*91.2*Kelvin*CII_abun*data['HI number density']*0.75*data['H2 number density']
yt.add_field(("gas", "CII_H2_ortho"), function=_CII_H2_ortho, units="erg/cm**3/s")

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
