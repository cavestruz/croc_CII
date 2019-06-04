//  Get dust attenuation: need HI density, HII density 
dust_atten =
  exp(-(cell_HI_density(losDataStorage[i].cell)+2*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length*soblenfloat*met_cell*2.0e-21);

//  Get temperature
tem = rtTem(losDataStorage[i].cell)*units->temperature;

//  CIIe and CIIa depends on temperature
rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
rC2a =
  1.0e-24*exp(-92./tem)*max2(0.0,16.+0.344*sqrt(tem)-47.7/tem);

//  CIIe and CIIa cooling depend on HII density, HeII density, HEIII density, and HI density
  
C2_e_cooling=rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;

C2_a_cooling=rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
CII_HeI =
  0.38*rC2a*cell_HeI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;

//  Add to the cold or warm medium
if(tem <= 1.0e3)
  
  CII_CNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
 else
   
   CII_WNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

//  Calculated CII emission
CII_CMB_emission =
  2.0*C2_abun*cell_HI_density(losDataStorage[i].cell)*met_cell*units->number_density*2.291e-6*constants->k*91.2*exp(-91.2/2.725*abox[min_level]);



pressure+=cell_gas_pressure(losDataStorage[i].cell)*units->energy_density/constants->k*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

                     A_V +=
		       (cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*5.3e-22*met_cell*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

CII_H2_para =
  4.25e-10*pow(tem/100.0,0.124-0.018*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.25*cell_H2_density(cell)*units->number_density;

CII_H2_ortho =
  5.14e-10*pow(tem/100.0,0.124+0.023*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.75*cell_H2_density(cell)*units->number_density;

tem = 100.0;
rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
rC2a =
  1.0e-24*exp(-92/tem)*max2(0.0,16+0.344*sqrt(tem)-47.7/tem);
CII_H2_basic +=
  cell_H2_density(losDataStorage[i].cell)/(cell_gas_density(losDataStorage[i].cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell
												    +
												    rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell
												    +
												    1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sfr_in_cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
