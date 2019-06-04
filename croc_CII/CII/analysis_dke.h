//Headers

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h> 
#include <stdlib.h>
#include "/usr/include/cfitsio/fitsio.h"

#include "cosmology.h"
#include "iterators.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "units.h"
#include "tree.h"
#include "auxiliary.h"
#include "los_dke.h"

#include "times.h"
#include "rt.h"
#include "parallel.h"
#include "particle.h"

#include "control_parameter.h"
#include "hydro.h"
#include "config.h"

//From extra
//#include "plugin.h"
#include "healpix.h"
#include "halo_finder.h"
#include "CO.h"

#ifdef RADIATIVE_TRANSFER 
#include "frt/frt_c.h"
#include "frt/frt_index.h"
#endif

//Function Prototypes

void DEboxAnalysis();
void DEBoxDistributionAnalysis();
void DELOS_alongprojgrid();
void DELOS_alongprojgrid_generic();
void DELOS_alongprojgrid_multi_d_comparison();
void DELOS_alongprojgrid_multi_d_LOS_characteristic_length();
void DELOS_cross_section_vs_Mvir_file_approach();

  /***************************************************/
  /*                  CII Block                      */
  /***************************************************/

void CII_net_heating_cooling_etc_per_primary_galaxy();
void CII_emission_individual_galaxy_map();

void CII_3d_radial_profiles_Art_Format();
void CII_projected_3d_radial_profiles_Art_Format();

void CII_phase_diagrams_each_galaxy();
void CII_phase_diagrams_each_galaxy_Photoionization_vs_GPE();
void CII_individual_primary_halo_properties();

void CII_projected_radiation_map();
void CII_sliced_radiation_map();

void CII_background_spectrum();
void CII_radial_radiation_spectrum();
void CII_radial_radiation_spectrum_HI_cut();
void CII_radial_radiation_spectrum_HI_weighted();
void CII_radial_radiation_spectrum_T_cut();
void CII_radial_radiation_spectrum_height_cut();
void CII_radial_radiation_spectrum_scatter_plot();
void CII_radiation_field_check();

  /***************************************************/
  /*                  BQC Block                      */
  /***************************************************/

void BQC_covering_fraction_versus_radius_fast_version();
void BQC_covering_fraction_versus_radius_fast_version_random_subset();
void BQC_box_covering_fraction_by_halo();

void BQC_star_formation_catalog();

  /***************************************************/
  /*                   CO Block                      */
  /* This block requires the CO flag to be turned on */
  /***************************************************/

void CO_emission_trial();
void CII_net_heating_cooling_etc_per_primary_galaxy_wCO();
void CII_phase_diagrams_each_galaxy_Photoionization_vs_GPE_wCO();
void CII_emission_individual_galaxy_map_wCO();
void CII_map_with_galaxy_face_on();
void CII_map_with_galaxy_face_on_velocity_profiles();

void LOS_test();

  /***************************************************/
  /*                  LLS Block                      */
  /***************************************************/

void LLS_3d_distance_to_closest_primary_halo();
void LLS_3d_distance_to_closest_halo();
void LLS_3d_distance_to_closest_halo_r_trunk();
void LLS_Characteristic_Length();
void LLS_cross_section_vs_Mvir_maxima_approach();
void LLS_cross_section_vs_Mvir_maxima_approach_rtrunk();
void LLS_phase_diagrams_central_cell();
void LLS_NHI_histogram_for_different_lengths();
void LLS_Photoionization_versus_Recominbation();

  /***************************************************/
  /*                Metal Absorption                 */
  /***************************************************/

void MgII_face_on_profile_within_km_s();
void DEVEL_velocity_width();
void DEVEL_DLA_velocity_width();
void DEVEL_DLA_velocity_width_vs_halo_mass_vs_metallicity();
void HI_map_with_galaxy_face_on_velocity_cube();
void HI_map_with_galaxy_perp_velocity_cube();
void Circular_Velocity_By_Component();

  /**************************************************/
  /*       Halo and Stream Fragmentation            */
  /**************************************************/

void HSF_individual_galaxy_map();
void HSF_map_in_stream_plane();
void HSF_stream_cross_section();

  /**************************************************/
  /*            Analysis for others                 */
  /**************************************************/

void Sasha_recombination_rate_1();




//Constants used in running some of Sam's code

#define  LAM0_HI   1215.6701e-8
#define  LAM0_O6   1031.9261e-8
//#define  LAM0_O6b  1037.6167e-8
#define  LAM0_C4   1548.2041e-8
#define  LAM0_C3    977.0201e-8
#define  LAM0_N5   1238.8210e-8
#define  LAM0_Si4  1393.7602e-8
#define  LAM0_Si3  1206.5000e-8
#define  LAM0_Mg2  2796.3543e-8

#define EL_HI 0
#define EL_O6 1
#define EL_C4 2
#define EL_C3 3
#define EL_N5 4
#define EL_Si4 5
#define EL_Si3 6
#define EL_Mg2 7
#define MAX_ELEMENT 8




/* DO NOT CHANGE: */
#define BAD_CGRID_Z 2
#define ngz (13)
#define cgZlo  (-5.0)
#define cloudy_step_Z (0.5)
#define cgZhi ((ngz-1 - BAD_CGRID_Z)*cloudy_step_Z+cgZlo)
//#define cgZhi (0) //cloudy grid is BAD above 0

#define ngd (46)
#define cgrholo  (-6.0)
#define cloudy_step_rho (0.2)
#define cgrhohi ((ngd-1)*cloudy_step_rho+cgrholo)

#define ngt (26)
#define cgTlo  (2)
#define cloudy_step_T (0.2)
#define cgThi ((ngt-1)*cloudy_step_T+cgTlo)

#define npcld (ngz*ngd*ngt)
#define d3cld(K,J,I)  ((K)*ngt*ngd   + (J)*ngt  + (I)) 


//extern const double lam0_els[MAX_ELEMENT] ;
//extern const double lam0_els[MAX_ELEMENT] = { LAM0_HI, LAM0_O6 ,LAM0_C4 ,LAM0_C3 ,LAM0_Si4 ,LAM0_Si3 ,LAM0_Mg2 };

extern double cloudy_grid_dr;
double cloudy_grid_dr;

extern double nden_cloudy[MAX_ELEMENT][npcld];
double nden_cloudy[MAX_ELEMENT][npcld];

//extern  double nO6[npcld],nC4[npcld],nC3[npcld],nN5[npcld],nSi4[npcld],nSi3[npcld],nMg2[npcld];
//double nO6[npcld],nC4[npcld],nC3[npcld],nN5[npcld],nSi4[npcld],nSi3[npcld],nMg2[npcld];


#define eps (.0001)

//const double lam0_els[MAX_ELEMENT] = { LAM0_HI, LAM0_O6 ,LAM0_C4 ,LAM0_C3 , LAM0_N5 , LAM0_Si4 ,LAM0_Si3 ,LAM0_Mg2 };
