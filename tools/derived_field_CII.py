import yt as yt
from unyt import unyt_array
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


H_mass = 1.007825032*amu

H2_mass = 2.016*amu

He_mass = 4.002602*amu


### Added particle filed to filter out for young stars only

def young_stars(pfilter, data):
    age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
    filter = np.logical_and(age.in_units("Myr") <= 20, age >= 0)
    return filter


yt.add_particle_filter(
    "young_stars",
    function=young_stars,
    filtered_type="STAR",
    requires=["creation_time"],
)

### Conversions to YTEP-0003 compatible fields for use in Trident

# def _H_p0_density_YTEP_0003(field,data):
#     return data["HI density"].copy()
# yt.add_field(("gas", "H_p0_density"), function=_H_p0_density_YTEP_0003, units="g/cm**3",sampling_type="local")

# def _H_p1_density_YTEP_0003(field,data):
#     return data["HII density"].copy()
# yt.add_field(("gas", "H_p1_density"), function=_H_p1_density_YTEP_0003, units="g/cm**3",sampling_type="local")

# def _H2_p0_density_YTEP_0003(field,data):
#     return data["H2 density"].copy()
# yt.add_field(("gas", "H2_p0_density"), function=_H2_p0_density_YTEP_0003, units="g/cm**3",sampling_type="local")

# def _He_p0_density_YTEP_0003(field,data):
#     return data["HeI density"].copy()
# yt.add_field(("gas", "He_p0_density"), function=_He_p0_density_YTEP_0003, units="g/cm**3",sampling_type="local")

# def _He_p1_density_YTEP_0003(field,data):
#     return data["HeII density"].copy()
# yt.add_field(("gas", "He_p1_density"), function=_He_p1_density_YTEP_0003, units="g/cm**3",sampling_type="local")

# def _He_p2_density_YTEP_0003(field,data):
#     return data["HeIII density"].copy()
# yt.add_field(("gas", "He_p2_density"), function=_He_p2_density_YTEP_0003, units="g/cm**3",sampling_type="local")

### Functions from cII_emission.c

def _test_mass(field, data):
    return data["H_density"]/data["H_density"]*(1.007825032*amu)
yt.add_field(("gas", "test_mass"), function=_test_mass, units="g",sampling_type="local")

def _HI_number_density(field, data):
                 return data["H_density"]/(H_mass)
# yt.add_field(("gas", "HI number density"), function=_HI_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "H_number_density"), function=_HI_number_density, units="1/cm**3",sampling_type="local") # Adds additional field for YTEP 0003 and Trident compatability

def _HII_number_density(field, data):
    return data["H_p1_density"]/(H_mass-electron_mass)
# yt.add_field(("gas", "HII number density"), function=_HII_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "H_p1_number_density"), function=_HII_number_density, units="1/cm**3",sampling_type="local")# Adds additional field for YTE,sampling_type=",sampling_type="local"local"P 0003 and Trident compatability

def _H2_number_density(field, data):
    return data["H2_density"]/(H2_mass)
# yt.add_field(("gas", "H2 number density"), function=_H2_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "H2_number_density"), function=_H2_number_density, units="1/cm**3",sampling_type="local")# Adds additional field for YTEP 0003 and Trident compatability

def _HeI_number_density(field, data):
    return data["He_density"]/(He_mass)
# yt.add_field(("gas", "HeI number density"), function=_HeI_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "He_number_density"), function=_HeI_number_density, units="1/cm**3",sampling_type="local")# Adds additional field for YTEP 0003 and Trident compatability

def _HeII_number_density(field, data):
    return data["He_p1_density"]/(He_mass-electron_mass)
# yt.add_field(("gas", "HeII number density"), function=_HeII_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "He_p1_number_density"), function=_HeII_number_density, units="1/cm**3",sampling_type="local")# Adds additional field for YTEP 0003 and Trident compatability

def _HeIII_number_density(field, data):
    return data["He_p2_density"]/(He_mass-2*electron_mass)
# yt.add_field(("gas", "HeIII number density"), function=_HeIII_number_density, units="1/cm**3",sampling_type="local")
yt.add_field(("gas", "He_p2_number_density"), function=_HeIII_number_density, units="1/cm**3",sampling_type="local")# Adds additional field for YTEP 0003 and Trident compatability

def _log_dust_attenuation(field, data):
    return -(data["H_number_density"]+2*data["H2_number_density"])*data['metallicity']*L_soblen*dust_cross_section # With H2
#    return -(data["HI number density"])*data['metallicity']*L_soblen*dust_cross_section # Without H2
yt.add_field(("gas", "log_dust_attenuation"), function=_log_dust_attenuation, units="",sampling_type="local")

def _rCIIe(field, data):
    return 6.67*10**(-20)*Kelvin**(0.5)*np.exp(-91.2*Kelvin/data['temperature'])/data["temperature"]**0.5*erg*centimeter**3/second #Not sure which one contains the per second and why there is a square root of temperature
yt.add_field(("gas", "rCIIe"), function=_rCIIe, units="erg*cm**3/s",sampling_type="local")

def _rCIIa(field, data):
    x = 16 + .344*(data['temperature']/Kelvin)**(0.5) - 47.7*Kelvin/data['temperature']
    x[x<0.0]=0
    return 10**(-24)*np.exp(-91.2*Kelvin/data['temperature'])*x*erg*centimeter**3/second #Not sure which one contains the per second and why there is a square root of temperature                
yt.add_field(("gas", "rCIIa"), function=_rCIIe, units="erg*cm**3/s",sampling_type="local")

def _CII_e_cooling(field, data):
    return data['rCIIe'] * (data['H_p1_number_density']+data['He_p1_number_density']+2*data['He_p2_number_density']) * CII_abun * data['H_number_density'] * data['metallicity']
yt.add_field(("gas", "CII_e_cooling"), function=_CII_e_cooling, units="erg/cm**3/s",sampling_type="local")

def _CII_a_cooling(field, data):
    return data['rCIIa'] * CII_abun * data['H_number_density'] * data['H_number_density'] * data['metallicity']
yt.add_field(("gas", "CII_a_cooling"), function=_CII_a_cooling, units="erg/cm**3/s",sampling_type="local")

def _CII_HeI_cooling(field, data):
    return 0.38 * data['rCIIa'] * data['He_number_density'] * CII_abun * data['H_number_density'] * data['metallicity']
yt.add_field(("gas", "CII_HeI_cooling"), function=_CII_HeI_cooling, units="erg/cm**3/s",sampling_type="local")

def _CII_CMB_emission(field, data):
    z = data.ds.current_redshift
    return 2.0 *  CII_abun * data['H_number_density'] * data['metallicity']*2.298*10**(-6)*kboltz*91.2*Kelvin*np.exp(-91.2*Kelvin/(2.725*Kelvin)*1/(1+z))/second # not sure which one contains the per second unit
yt.add_field(("gas", "CII_CMB_emission"), function=_CII_CMB_emission, units="erg/cm**3/s",sampling_type="local")

def _CII_H2_para(field, data):
    """                                                                                                               
    These are cooling rates, i.e. Lambda*n^2 where Lambda is the cooling function. See https://www.astro.umd.edu/~richard/ASTRO620/A620_2015_Gas_lec2.pdf                  
    """
    rate_coefficient = 4.25*10**(-10)*centimeter**3/second*unyt_array(pow(data['temperature'].to_ndarray()/(100),0.124-0.018*np.log(data['temperature'].to_ndarray()/(100))))# Draine Table F6 pg 501
    return rate_coefficient*kboltz*(91.2*Kelvin)*CII_abun*data['H_number_density']*0.25*data['H2_number_density']
yt.add_field(("gas", "CII_H2_para"), function=_CII_H2_para, units="erg/cm**3/s",sampling_type="local")

def _CII_H2_ortho(field, data):
    """                                                                                                               
    These are cooling rates, i.e. Lambda*n^2 where Lambda is the cooling function.. See https://www.astro.umd.edu/~richard/ASTRO620/A620_2015_Gas_lec2.pdf                 
    """
    rate_coefficient = (5.14*10**(-10)*centimeter**3/second)*unyt_array(pow(data['temperature'].to_ndarray()/(100),0.095+0.023*np.log(data['temperature'].to_ndarray()/(100))))# Draine Table F6 pg 501
    return rate_coefficient*kboltz*91.2*Kelvin*CII_abun*data['H_number_density']*0.75*data['H2_number_density']
yt.add_field(("gas", "CII_H2_ortho"), function=_CII_H2_ortho, units="erg/cm**3/s",sampling_type="local")


def _LCII_e(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['CII_e_cooling']*(data['H_p1_number_density']+data['He_p1_number_density']+2*data['He_p2_number_density']) *data['cell_volume']
    
yt.add_field(("gas", "LCII_e"), function=_LCII_e, units="erg/cm**3/s",sampling_type="local")

def _LCII_a(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['CII_a_cooling']*(data['H_number_density'])*data['cell_volume']
    
yt.add_field(("gas", "LCII_a"), function=_LCII_a, units="erg/cm**3/s",sampling_type="local")

def _LCII_HeI(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['CII_HeI_cooling']*(data['He_number_density'])*data['cell_volume']
    
yt.add_field(("gas", "LCII_HeI"), function=_LCII_HeI, units="erg/cm**3/s",sampling_type="local")

def _LCII_CMB(field, data):
    """                                                                                                               
    CII Point Luminosity of CMB
    """
    return data['CII_CMB_emission']*(413 /(centimeter**3))*data['cell_volume']
    
yt.add_field(("gas", "LCII_CMB"), function=_LCII_CMB, units="erg/cm**3/s",sampling_type="local")



def _LCII_H2_ortho(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['CII_H2_ortho']*0.75*(data['H2_number_density'])*data['cell_volume']
    
yt.add_field(("gas", "LCII_H2_ortho"), function=_LCII_H2_ortho, units="erg/cm**3/s",sampling_type="local")

def _LCII_H2_para(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['CII_H2_para']*0.25*(data['H2_number_density'])*data['cell_volume']
    
yt.add_field(("gas", "LCII_H2_para"), function=_LCII_H2_para, units="erg/cm**3/s",sampling_type="local")

def _LCII_total(field, data):
    """                                                                                                               
    CII Point Luminosity of electrons
    """
    return data['LCII_e']+data['LCII_a']+data['LCII_HeI']+data['LCII_CMB']+data['LCII_H2_ortho']+data['LCII_H2_para']
    
yt.add_field(("gas", "LCII_total"), function=_LCII_total, units="erg/cm**3/s",sampling_type="local")




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
