// z=5 box is pretty small.   Raziq needs account on Nanna (can use Nanna as a computer).
// Step 1: make maps/slices.  Density and temperature, then calculate CII

// This project is semi-exploratory.  If cooling radiation dominates -
// diffuse ISM - CROC is sufficient.  IF dominated by PDR, then we
// need 10pc scale.  The results of this project tell you where to go
// from here.  Is higher resolution needed?  Followup projects depend
// on results.

//  See textbook by Bruce Draine (book ~2005) ISM physics.  Use as a
//  comprehensive reference, not a learning tool.

//  Get dust attenuation: need HI density, HII density Note: soblenfloat (see Nick's paper to turn into column density - rho gradient), 2.0e-21 is the cross section for dust
dust_atten =
  exp(-(cell_HI_density(losDataStorage[i].cell)+2*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length*soblenfloat*met_cell*2.0e-21);

//  Get temperature
tem = rtTem(losDataStorage[i].cell)*units->temperature;

//  CIIe and CIIa depends on temperature - these are rates, likely units of erg*cm^3/s (emission), of electrons (e) and atoms (a).  These rates are in Simon Glover's papers - Glover+Abel(2010-2015)  Ultimately care about total emission, spectrum secondary.  Note - when doing CII, should be able to do OIII (comes from regions like PDR)
rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
rC2a =
  1.0e-24*exp(-92./tem)*max2(0.0,16.+0.344*sqrt(tem)-47.7/tem);

//  CIIe and CIIa cooling depend on HII density, HeII density, HEIII
//  density, and HI density. Multiplying rates by density of
//  electrons.  C2 abundance*HIdensity is just CII abundance per unit
//  metal, then multiplied by metallicity.
  
C2_e_cooling=rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;

// Collisions with hydrogen atoms
C2_a_cooling=rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;

// Collisions with HeI, and the rate is .38 times the rate of
// collisions with Hydrogen.  (Good to review these rates - email
// Raziq.)
CII_HeI =
  0.38*rC2a*cell_HeI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;

//  Add to the cold or warm neutral medium
/* if(tem <= 1.0e3) */
  
/*   CII_CNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1); */
/*  else */
   
/*    CII_WNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1); */



// Must be collisions with CMB photons --> excitations by CMB (see
// term with CMB temperature) Need to check this process, why is there
// no temperature dependence here(?)  May be assumption that the
// medium is always hotter than CMB, which is ok assumption at
// intermediate redshifts.
CII_CMB_emission =
  2.0*C2_abun*cell_HI_density(losDataStorage[i].cell)*met_cell*units->number_density*2.291e-6*constants->k*91.2*exp(-91.2/2.725*abox[min_level]);



/* pressure+=cell_gas_pressure(losDataStorage[i].cell)*units->energy_density/constants->k*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1); */

/* A_V += */
/*   (cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*5.3e-22*met_cell*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1); */

// Collision rates with molecular hydrogen in para and ortho state
// (aligned and misaligned spin with hydrogen).  Unless in high
// densities >1e4, these are like two different molecules.  Transition
// between these only happens at high density.  .25 - assume 1/4 in
// para state, and .75 3/4 in ortho - typical assumption for thermal
// equilibrium.
CII_H2_para =
  4.25e-10*pow(tem/100.0,0.124-0.018*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.25*cell_H2_density(cell)*units->number_density;

CII_H2_ortho =
  5.14e-10*pow(tem/100.0,0.124+0.023*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.75*cell_H2_density(cell)*units->number_density;

/* tem = 100.0; */
/* rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem); */
/* rC2a = */
/*   1.0e-24*exp(-92/tem)*max2(0.0,16+0.344*sqrt(tem)-47.7/tem); */
/* CII_H2_basic += */
/*   cell_H2_density(losDataStorage[i].cell)/(cell_gas_density(losDataStorage[i].cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell */
/* 												    + */
/* 												    rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell */
/* 												    + */
/* 												    1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sfr_in_cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1); */


sfr_in_cell = 0.0;
ipart = cell_particle_list[losDataStorage[i].cell];

while(ipart != NULL_PARTICLE)
  {
    cart_assert( ipart >= 0 && ipart < num_particles );
    if(particle_is_star(ipart))
      {
	stellar_surface_density+=particle_mass[ipart]*units->mass/constants->Msun*pow(units->length/po \
										      w(2.0,cell_level(losDataStorage[i].cell)),-3.0)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(consta \
																									      nts->pc,2.0);
	
	if((particle_t[ipart]-(double)star_tbirth[ipart])*units->time/constants->yr < 2.0e7)
	  {
	    // in Msun/yr which is what we need in the PDR calculation                                 
	    sfr_in_cell += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7);
	    sfr_per_area += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7)/( \
												      pow(units->length/constants->pc*pow(2.0,-(double)cell_level(losDataStorage[i].cell)),3.0))*(min2(losDataStorage[i].r2,len)-los \
																								  DataStorage[i].r1)*units->length/constants->pc;
	  }
      }
    ipart = particle_list_next[ipart];
  }

PDR_C2cooling = 0.;
        if((cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length/pow(2.0,cell\
_level(losDataStorage[i].cell))*met_cell > 1.0e21)
                      PDR_C2cooling=1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*sfr_in_cell/pow(units->length/pow(2.0,cell_level(losDataStorage[i].cel\
l)),3.0);


                    losHIdensity+=cell_HI_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losH2density+=cell_H2_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losgasdensity+=cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losmetdensity+=met_cell*constants->Zsun*cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);


                    //this velocity is in km/s                                                                                                                
                    projected_momentum = cell_momentum(losDataStorage[i].cell,0)*v3[0]+cell_momentum(losDataStorage[i].cell,1)*v3[1]+cell_momentum(losDataSto\
rage[i].cell,2)*v3[2];
                    projected_halo_velocity = (double)ngHalos->list[counter].vel[0]*v3[0]+(double)ngHalos->list[counter].vel[1]*v3[1]+(double)ngHalos->list[c\
ounter].vel[2]*v3[2];
   velocity = Hubble(abox[min_level])*(losDataStorage[i].r2+losDataStorage[i].r1)/2.0*units->length/constants->Mpc+projected_momentum/cell_g\
as_density(losDataStorage[i].cell)*units->velocity*1e-5-projected_halo_velocity;

                    //by convention we define the first direction away from us, and the other direction towards us.                                           

                    if(velocity > min_velocity && velocity < max_velocity)
                      {
                        optical_depth[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling)*(min2(losDataS\
torage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);

                        optical_depth_all[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling+CII_HeI+PDR\
_C2cooling+CII_CMB_emission+CII_H2_para+CII_H2_ortho)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);

                      }

                  }
              }
            for(i=0;i<(3*(int)ceil(pow(2.0,NMAXLOS)*len));i++)
              {
                losDataStorage[i].r2=0;
              }
            ///////////////// NOTE WE ARE ASSUMING THAT THE FIRST HALF OF PIX DOES NOT INCLUDE ANY ANTIPODAL PAIRS                                            
            //Now again for the antipodal points.                                                                                                             
            //We want to be careful that we do not double count the originating point for these LOS                                                           

            segmentvar.Buffer.Data = losDataStorage;
            segmentvar.Buffer.Size = sizeof(losdatatype);
            losTraverseSegment_dke(0,position,M_PI-theta,M_PI+phi,len,NMAXLOS,WorkerBee,&segmentvar);

            for(i=0;i<(3*(int)ceil(pow(2.0,NMAXLOS)*len));i++)
              {
                //This if used to be if(...>1) so that we missed the unresolved regions defaulted to level 0,                                                 
                if(losDataStorage[i].r2!=0 && cell_level(losDataStorage[i].cell) >= NMIN)
                  {
                    if(cell_gas_metal_density(losDataStorage[i].cell)/cell_gas_density(losDataStorage[i].cell)*1/constants->Zsun > met_floor)
                      met_cell = cell_gas_metal_density(losDataStorage[i].cell)/cell_gas_density(losDataStorage[i].cell)*1/constants->Zsun;
                    else
                      met_cell = met_floor;

                    rtGetCoolingRate(losDataStorage[i].cell,&coolingrate,&heatingrate);
                    soblenfloat = cell_sobolev_length(losDataStorage[i].cell);
                    rtGetPhotoRatesFS(losDataStorage[i].cell,rates);
                    dust_atten = exp(-(cell_HI_density(losDataStorage[i].cell)+2*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length\
*(double)soblenfloat*met_cell*2.0e-21);

                    tem = rtTem(losDataStorage[i].cell)*units->temperature;

                    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
                    rC2a =  1.0e-24*exp(-92./tem)*max2(0.0,16.+0.344*sqrt(tem)-47.7/tem);

                    C2_e_cooling=rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStora\
ge[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
                    C2_a_cooling=rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->nu\
mber_density*met_cell;

                    if(tem <= 1.0e3)
                      CII_CNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    else
                      CII_WNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

                    CII_CMB_emission= 2.0*C2_abun*cell_HI_density(losDataStorage[i].cell)*met_cell*units->number_density*2.291e-6*constants->k*91.2*exp(-91.2\
/2.725*abox[min_level]);

                    pressure+=cell_gas_pressure(losDataStorage[i].cell)*units->energy_density/constants->k*(min2(losDataStorage[i].r2,len)-losDataStorage[i].\
r1);

                    A_V += (cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*5.3e-22*met_cell*(min2\
(losDataStorage[i].r2,len)-losDataStorage[i].r1);

                    CII_H2_para = 4.25e-10*pow(tem/100.0,0.124-0.018*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.\
25*cell_H2_density(cell)*units->number_density;

                    CII_H2_ortho = 5.14e-10*pow(tem/100.0,0.124+0.023*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0\
.75*cell_H2_density(cell)*units->number_density;

                    tem = 100.0;
                    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
                    rC2a =  1.0e-24*exp(-92/tem)*max2(0.0,16+0.344*sqrt(tem)-47.7/tem);
     CII_H2_basic += cell_H2_density(losDataStorage[i].cell)/(cell_gas_density(losDataStorage[i].cell)*constants->XH)*(rC2e*(cell_HII_density(\
losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_densi\
ty(losDataStorage[i].cell)*units->number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDat\
aStorage[i].cell)*units->number_density*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sfr_in_cell)*(min2(losDataStorage[i].r2,len)-losDataSt\
orage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    //CII_H2_basic += cell_H2_density(cell)/(cell_gas_density(cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_\
density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->\
number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_de\
nsity*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sf_recipe1_rate(losDataStorage[i].cell)*units->density/units->time*1/constants->Msun*con\
stants->yr)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);                                     
                    sfr_in_cell = 0.0;
                    ipart = cell_particle_list[losDataStorage[i].cell];

                    while(ipart != NULL_PARTICLE)
                      {
                        cart_assert( ipart >= 0 && ipart < num_particles );
                        if(particle_is_star(ipart))
                              {
                                stellar_surface_density+=particle_mass[ipart]*units->mass/constants->Msun*pow(units->length/pow(2.0,cell_level(losDataStorage\
[i].cell)),-3.0)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
                                if((particle_t[ipart]-(double)star_tbirth[ipart])*units->time/constants->yr < 2.0e7)
                                  {
                                    // in Msun/yr                                                                                                             
                                    sfr_in_cell += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7);
                                    sfr_per_area += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7)/(pow(units->length/constants->pc\
*pow(2.0,-(double)cell_level(losDataStorage[i].cell)),3.0))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length/constants->pc;

                                  }


                              }
                        ipart = particle_list_next[ipart];
                      }

                    PDR_C2cooling = 0.;
                    // For things which are not resolved, make                                   
                    // assumptions about PDR (photodissociation                                                                                               
                    // region, outer HII region - not yet ionized, but                                                                                        
                    // molecules are already dissociating) in star                                                                                            
                    // forming regions (Carbon II).  This is a subgrid                                                                                        
                    // model treatment for less resolved physical                                                                                             
                    // pieces.  Comparison of different pieces will be                                                                                        
                    // a part of the storyline (which dominates, etc.)                                                                                        

                    //  Atomic number density criterion - likely based                                                                                        
                    //  on some sort of fitting formula, but need                                                                                             
                    //  literature source for this.                                                                                                           
                    if((cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length/pow(2.0,cell\
_level(losDataStorage[i].cell))*met_cell > 1.0e21)
                      PDR_C2cooling=1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*sfr_in_cell/pow(units->length/pow(2.0,cell_level(losDataStorage[i].cel\
l)),3.0);
                    losHIdensity+=cell_HI_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losH2density+=cell_H2_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losgasdensity+=cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
                    losmetdensity+=met_cell*constants->Zsun*cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

                    // To get spectrum                                                                                                                        
                    projected_momentum = cell_momentum(losDataStorage[i].cell,0)*v3[0]+cell_momentum(losDataStorage[i].cell,1)*v3[1]+cell_momentum(losDataSto\
rage[i].cell,2)*v3[2];


                    projected_halo_velocity = (double)ngHalos->list[counter].vel[0]*v3[0]+(double)ngHalos->list[counter].vel[1]*v3[1]+(double)ngHalos->list[c\
ounter].vel[2]*v3[2];
                    // Add hubble component to get the proper velocity                                                                                        
                    velocity = -Hubble(abox[min_level])*(losDataStorage[i].r2+losDataStorage[i].r1)/2.0*units->length/constants->Mpc + projected_momentum/cel\
l_gas_density(losDataStorage[i].cell)*units->velocity*1e-5-projected_halo_velocity;

                  // If within range, add to the optical depth - two                                                                                        
                    // different optical depths.  First is only                                                                                               
                    // cooling (atomic gas contributions), second is                                                                                          
                    // all other contributions.  CII - when gas cools,                                                                                        
                    // CII is one of the lines from cooling radiation,                                                                                        
                    // then there is molecular radiation, then also                                                                                           
                    // PDR.  Need to calculate separate cooling                                                                                               
                    // contributions to do a budget of contributions                                                                                          
                    // to CII.  --> Create CII emission in each cell.                                                                                         
                    // YT should be able to do an emission spectrum.                                                                                          

                    if(velocity > min_velocity && velocity < max_velocity)
                      {
                        optical_depth[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling)*(min2(losDataS\
torage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
                        optical_depth_all[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling+CII_HeI+PDR\
_C2cooling+CII_CMB_emission+CII_H2_para+CII_H2_ortho)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
                      }

                  }

