//CII for carbon emission from galaxies and ISM physics

#include "analysis.h"
#include "analysis_dke.h"

// levels used in the analysis

#define NMIN 2
#define NMAX 9

// maximum level stepped through along a line of sight - should be set to NMAX

#define NMAXLOS NMAX

// root grid size

#define RGSIZE 256.0 //64 for 6Mpc, 256 for 25Mpc box

#define CII_NUM_HALOS 10 // number of halos to make maps from in the halo catalog

#define LOS_Projection_Length 50. // length of the LOS in each direction from the galaxy we will integrate along (in physical kpc), i.e. 25 would give a 50kpc line of sight
  
#define LOS_Side_Length 25. // side length of the square we will project along (in physical kpc)

int WorkerBee(int id, int cell, double r1, double r2, losBuffer *data);

typedef struct{
  double r1;  
  double r2;   
  int cell;
}losdatatype;

double variablequantity(int cell,int value);
char* variablequantityname(int value); 
char* variableaxesname(int value); 
int variableaxisone(int value);
int variableaxistwo(int value);
int projected_axis(int value);
double distance(double position1[3],double position2[3]);
double distance_cartesian(double position1[3],double position2[3]);
double distance_cartesian_2d(double position1[2],double position2[2]);
double twoDdistance(double position1[2],double position2[2]);
double projected_distance(double vector1[3], double vector2[3]);
double min3(double one, double two, double three);
double min2(double one, double two);
double max2(double one, double two);
double max3(double one, double two, double three);
double lognormaldist(double mu);
double lognormaldist_plaw(double mu);
double machnumber(int cell);
double distance_to_nearest_halo_at_or_above_level(double cell_position[3],int level,int axescounter);
int level_of_halo_at_position(double position[3]);
int dust_opacity_check(double tau);
void find_primary_halos(int *primary_halo_id, int *primary_halo_number);
double velocity_width_from_optical_depth(double *optical_depth, int RES);
void cross_product(double *x, double *y, double *result);
double kronecker_delta(int i,int j);
double mod_DE(double x, double y);
double mod_map(double x);

int double_cmp(const void *a, const void *b);
double cic3d(double x1, double x2, double x3, double *data);
void read_cloudy_grid();
void cloudy_grid_check(double xcgZ,double xcgrho,double xcgT);
double cloudy_column(double zmet, double rho_nH, double gtemp, double element[npcld]);

int number_of_cells_in_annulus(double r_in, double r_out);
double jeans_length(int cell);
double jeans_length_no_turb(int cell);

int comp_multid_array(const void* a, const void* b) {

  double* da = (double*)a;

  double* db = (double*)b;

  int diff1 = (da[4] > db[4]) - (da[4] < db[4]);

  return -diff1;

}

const double lam0_els[MAX_ELEMENT] = { LAM0_HI, LAM0_O6 ,LAM0_C4 ,LAM0_C3 , LAM0_N5 , LAM0_Si4 ,LAM0_Si3 ,LAM0_Mg2 };

/* Note: rtGetPhotoRates(x,x,0) is now rtGetPhotoRates and returns the shielded rate while rtGetPhotoRates(x,x,1) is now rtGetPhotoRatesFS and returns the free space field */

/* Careful with rtUmw2, it used to mean a shielded field */



void CII_map_with_galaxy_face_on_velocity_profiles()
{
  printf("Entering CII_map_with_galaxy_face_on_velocity_profiles(), a= %.4f\n",abox[min_level]);
  printf("Hubble constant = %f km/s/Mpc\n",Hubble(abox[min_level]));

  halo_list *ngHalos = NULL;

  char str[999];
  int closest_halo_id;
  int ipart = 0;
  
  sprintf(str,"DAT_HC/HC/hlist_%6.4f.dat",auni[min_level]);
  ngHalos = load_halo_finder_catalog(str,50,0.0,0.0,0.0,100000);

  double position[3],position_halo[3],delta_r[3];

  double momentum[3],DeltaL[3],L_tot[3];
  double theta,phi;

  double v1[3],v2[3],v3[3];

  int counter,i,j;

  int v1_counter,v2_counter;

  losSegment segmentvar; 
  losdatatype *losDataStorage;
  double len, len_side;

  double losHIdensity;
  double losH2density;
  double losgasdensity;
  double losmetdensity;

  FILE *dataout,*dataout_all;
  char fname[999];

  double L_inclusion_length=3.0; //in kpc, for angular momentum calculation, used to be 10kpc
  double Projection_length=LOS_Projection_Length; //in kpc, side of cube we will be projecting 
  double Side_length = LOS_Side_Length;

  double met_cell,met_floor = 1.0e-3;
  float heatingrate,coolingrate;
  float soblenfloat,sobvelfloat;
  float rates[FRT_RATE_DIM];
  double dust_atten;
  double tem;
  double rC2e,rC2a;
  double C2_a_cooling,C2_e_cooling;
  double CII_CNM;
  double CII_WNM;
  double PDR_C2cooling;
  double CII_CMB_emission;
  double C2_abun = 3.31e-4;
  double pressure;
  double A_V; //From Draine's book, page 240.
  double CII_H2_basic;
  double CII_H2_collision;
  double CII_H2_para;
  double CII_H2_ortho;
  double CII_HeI;
  double stellar_surface_density;
  double sfr_per_area;
  double sfr_in_cell;
  double total_sfr;
  double total_flux;
  double CII_v,CII_v2,CII_sum,CII_vwidth;



  double velocity, projected_momentum,projected_halo_velocity;
  double min_velocity = -1000.0;
  double max_velocity = 1000.0;
  int VEL_RES = 100;
  double *optical_depth = malloc(VEL_RES*sizeof(double));
  if(optical_depth == NULL)
    {
      printf("Not enough memory for optical depth\n");
      exit(1);
    }
  double *optical_depth_all = malloc(VEL_RES*sizeof(double));
  if(optical_depth_all == NULL)
    {
      printf("Not enough memory for optical_depth_all\n");
      exit(1);
    }


  double bufferout;
  
  double *flux_profile = malloc(VEL_RES*sizeof(double));
  if(flux_profile == NULL)
    {
      printf("Not enough memory for flux profile\n");
      exit(1);
    }

  double *flux_profile_all = malloc(VEL_RES*sizeof(double));
  if(flux_profile_all == NULL)
    {
      printf("Not enough memory for flux profile all\n");
      exit(1);
    }


  double grain_eff=0.001; //grain ionization efficiency
  double frac_abs = 0.3;  //fraction of FUV absorbed
  
  int ret;

  rtUpdateTables(min_level,mpi.comm.run);
  
  for(counter = 0;counter < CII_NUM_HALOS;counter++)
    {
      for(i=0;i<VEL_RES;i++)
	{
	  flux_profile[i]=0.0;
	  flux_profile_all[i] = 0.0;
	}

      position_halo[0] = ngHalos->list[counter].pos[0];
      position_halo[1] = ngHalos->list[counter].pos[1];
      position_halo[2] = ngHalos->list[counter].pos[2];
      
      L_tot[0] = L_tot[1] = L_tot[2] = 0.0;
      
      MESH_RUN_DECLARE(level, cell);
      MESH_RUN_OVER_LEVELS_BEGIN(level,NMIN,NMAX);
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
      if(cell_is_leaf(cell)||level==NMAX)
	{
	  cell_center_position(cell,position);
	  if(distance(position,position_halo) < L_inclusion_length*constants->kpc/units->length)
	    {
	      delta_r[0] = position[0]-position_halo[0];
	      delta_r[1] = position[1]-position_halo[1];
	      delta_r[2] = position[2]-position_halo[2];

	      momentum[0] = (cell_HI_fraction(cell)+2.0*cell_H2_fraction(cell))*cell_momentum(cell,0);
	      momentum[1] = (cell_HI_fraction(cell)+2.0*cell_H2_fraction(cell))*cell_momentum(cell,1);
	      momentum[2] = (cell_HI_fraction(cell)+2.0*cell_H2_fraction(cell))*cell_momentum(cell,2);
	      
	      cross_product(delta_r,momentum,DeltaL);

	      L_tot[0]+=DeltaL[0];
	      L_tot[1]+=DeltaL[1];
	      L_tot[2]+=DeltaL[2];
	    }
	}
      MESH_RUN_OVER_CELLS_OF_LEVEL_END;
      printf("#Level %d completed\n",level);
      MESH_RUN_OVER_LEVELS_END;

      phi = atan(L_tot[1]/L_tot[0]);
      theta = atan(sqrt(L_tot[0]*L_tot[0]+L_tot[1]*L_tot[1])/L_tot[2]);

      printf("Halo %d points in %f, %f\n",counter,theta,phi);

      v1[0]=sin(phi); //v1 and v2 are the two directions perpendicular to the line of sight
      v1[1]=-cos(phi);
      v1[2]=0.0;

      v2[0]=cos(theta)*cos(phi);
      v2[1]=cos(theta)*sin(phi);
      v2[2]=-sin(theta);

      v3[0]=sin(theta)*cos(phi); //v3 is the direction along the line of sight
      v3[1]=sin(theta)*sin(phi);
      v3[2]=cos(theta);

      sprintf(fname,"a%6.4f/CII",auni[min_level]);
      ret = system_mkdir(fname);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",fname,ret);
	  exit(1);
	}

      sprintf(fname,"a%6.4f/CII/by_halo_map",auni[min_level]);
      ret = system_mkdir(fname);
      if(ret != 0)
        {
          cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",fname,ret);
          exit(1);
        }


      sprintf(fname,"CII/by_halo_map/N%d_%d_face_on_velocity_map_halo_%d.bin",NMIN,NMAX,counter);
      dataout = fopen(ngOutputFile(fname),"wb");
      
      sprintf(fname,"CII/by_halo_map/N%d_%d_face_on_velocity_map_halo_%d_all.bin",NMIN,NMAX,counter);
      dataout_all = fopen(ngOutputFile(fname),"wb");

      len = Projection_length*constants->kpc/units->length;//ngHalos->list[counter].rvir;

      len_side = Side_length*constants->kpc/units->length;
      
      losDataStorage = malloc((3*(int)(ceil(pow(2.0,NMAX)*2*len)*sizeof(losdatatype))));
      //Factor of 2 just to be case since paths can pass through many extra cells if off diagonal.
      if(losDataStorage==NULL)
	{
	  printf("Failed to assign memory for losdata: Critical Failure!\n");
	  exit(1);
	}

      total_sfr = 0.0;
      
      for(v1_counter=-(int)(len_side/2.0*pow(2.0,NMAX));v1_counter< (int)(len_side/2.0*pow(2.0,NMAX));v1_counter++)
	for(v2_counter=-(int)(len_side/2.0*pow(2.0,NMAX));v2_counter< (int)(len_side/2.0*pow(2.0,NMAX));v2_counter++)
	  {
	    losHIdensity=0.0;
	    losH2density=0.0;
	    losgasdensity=0.0;
	    losmetdensity=0.0;
	    CII_CNM=0.0;
	    CII_WNM=0.0;
	    PDR_C2cooling=0.0;
	    CII_CMB_emission=0.0;
	    pressure = 0.0;
	    A_V = 0.0;
	    CII_H2_basic = 0.0;
	    CII_H2_collision = 0.0;
	    CII_H2_ortho = CII_H2_para = 0.0;
	    CII_HeI = 0.;
	    stellar_surface_density = 0.0;
	    sfr_per_area = 0.0;

	    for(i=0;i<VEL_RES;i++)
	      optical_depth[i]=optical_depth_all[i]=0.0;
	    

	    for(i=0;i<(3*(int)ceil(pow(2.0,NMAXLOS)*len));i++)
	      {
		losDataStorage[i].r2=0;
	      }
	    
	    segmentvar.Buffer.Data = losDataStorage; 
	    segmentvar.Buffer.Size = sizeof(losdatatype);
	    
	    //printf("#In direction theta = %f, phi = %f,memory starting %d\n",theta,phi,losDataStorage);
	    
	    position[0]=mod_DE(position_halo[0]+(v1_counter*v1[0]+v2_counter*v2[0])/pow(2.0,NMAX),RGSIZE);
	    position[1]=mod_DE(position_halo[1]+(v1_counter*v1[1]+v2_counter*v2[1])/pow(2.0,NMAX),RGSIZE);
	    position[2]=mod_DE(position_halo[2]+(v1_counter*v1[2]+v2_counter*v2[2])/pow(2.0,NMAX),RGSIZE);

	    losTraverseSegment_dke(0,position,theta,phi,len,NMAXLOS,WorkerBee,&segmentvar);
	    
	    //NOTE: The original cell is in losDataStorage for both forwards and backwards so don't double count it.
	    for(i=0;i<(3*(int)ceil(pow(2.0,NMAXLOS)*len));i++)
	      {
		//This if used to be if(...>1) so that we missed the unresolved regions defaulted to level 0,1
		if(losDataStorage[i].r2!=0 && cell_level(losDataStorage[i].cell) >= NMIN)
		  {
		    if(cell_gas_metal_density(losDataStorage[i].cell)/cell_gas_density(losDataStorage[i].cell)*1/constants->Zsun > met_floor)
		      met_cell = cell_gas_metal_density(losDataStorage[i].cell)/cell_gas_density(losDataStorage[i].cell)*1/constants->Zsun;
		    else
		      met_cell = met_floor;
			  
		    rtGetCoolingRate(losDataStorage[i].cell,&coolingrate,&heatingrate);
		    soblenfloat = cell_sobolev_length(losDataStorage[i].cell);
		    rtGetPhotoRatesFS(losDataStorage[i].cell,rates);
			  
		    dust_atten = exp(-(cell_HI_density(losDataStorage[i].cell)+2*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length*soblenfloat*met_cell*2.0e-21);
			  
		    tem = rtTem(losDataStorage[i].cell)*units->temperature;
			  
		    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
		    rC2a =  1.0e-24*exp(-92./tem)*max2(0.0,16.+0.344*sqrt(tem)-47.7/tem);
			  
		    C2_e_cooling=rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
		    C2_a_cooling=rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
		    CII_HeI = 0.38*rC2a*cell_HeI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
		    
		    if(tem <= 1.0e3)
		      CII_CNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    else
		      CII_WNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
			  
			  		    		    
		    CII_CMB_emission = 2.0*C2_abun*cell_HI_density(losDataStorage[i].cell)*met_cell*units->number_density*2.291e-6*constants->k*91.2*exp(-91.2/2.725*abox[min_level]);   

		    pressure+=cell_gas_pressure(losDataStorage[i].cell)*units->energy_density/constants->k*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

		    A_V += (cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*5.3e-22*met_cell*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

		    CII_H2_para = 4.25e-10*pow(tem/100.0,0.124-0.018*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.25*cell_H2_density(cell)*units->number_density;

		    CII_H2_ortho = 5.14e-10*pow(tem/100.0,0.124+0.023*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.75*cell_H2_density(cell)*units->number_density;
		    
		    tem = 100.0;
		    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
		    rC2a =  1.0e-24*exp(-92/tem)*max2(0.0,16+0.344*sqrt(tem)-47.7/tem);
		    CII_H2_basic += cell_H2_density(losDataStorage[i].cell)/(cell_gas_density(losDataStorage[i].cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sfr_in_cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    //CII_H2_basic += cell_H2_density(cell)/(cell_gas_density(cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sf_recipe1_rate(losDataStorage[i].cell)*units->density/units->time*1/constants->Msun*constants->yr)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);


		    sfr_in_cell = 0.0;
		    ipart = cell_particle_list[losDataStorage[i].cell];
		    
		    while(ipart != NULL_PARTICLE)
		      {
			cart_assert( ipart >= 0 && ipart < num_particles );
			if(particle_is_star(ipart))
			      {
				stellar_surface_density+=particle_mass[ipart]*units->mass/constants->Msun*pow(units->length/pow(2.0,cell_level(losDataStorage[i].cell)),-3.0)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);

				if((particle_t[ipart]-(double)star_tbirth[ipart])*units->time/constants->yr < 2.0e7)
				  {
				    // in Msun/yr which is what we need in the PDR calculation
				    sfr_in_cell += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7);
				    sfr_per_area += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7)/(pow(units->length/constants->pc*pow(2.0,-(double)cell_level(losDataStorage[i].cell)),3.0))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length/constants->pc;
				  }
				
			      }
			ipart = particle_list_next[ipart];
		      }

		    PDR_C2cooling = 0.;

		    if((cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length/pow(2.0,cell_level(losDataStorage[i].cell))*met_cell > 1.0e21)
		      PDR_C2cooling=1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*sfr_in_cell/pow(units->length/pow(2.0,cell_level(losDataStorage[i].cell)),3.0);


		    losHIdensity+=cell_HI_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losH2density+=cell_H2_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losgasdensity+=cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losmetdensity+=met_cell*constants->Zsun*cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);


		    //this velocity is in km/s
		    projected_momentum = cell_momentum(losDataStorage[i].cell,0)*v3[0]+cell_momentum(losDataStorage[i].cell,1)*v3[1]+cell_momentum(losDataStorage[i].cell,2)*v3[2];
		    projected_halo_velocity = (double)ngHalos->list[counter].vel[0]*v3[0]+(double)ngHalos->list[counter].vel[1]*v3[1]+(double)ngHalos->list[counter].vel[2]*v3[2];
		    velocity = Hubble(abox[min_level])*(losDataStorage[i].r2+losDataStorage[i].r1)/2.0*units->length/constants->Mpc+projected_momentum/cell_gas_density(losDataStorage[i].cell)*units->velocity*1e-5-projected_halo_velocity;

		    //by convention we define the first direction away from us, and the other direction towards us.
		    
		    if(velocity > min_velocity && velocity < max_velocity)
		      {
			optical_depth[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);

			optical_depth_all[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling+CII_HeI+PDR_C2cooling+CII_CMB_emission+CII_H2_para+CII_H2_ortho)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);

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
	    
	    //printf("#In direction theta = %f, phi = %f,Memory starting %d\n",M_PI-theta,-M_PI+phi,losDataStorage);
	    
	    
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
			  
		    dust_atten = exp(-(cell_HI_density(losDataStorage[i].cell)+2*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length*(double)soblenfloat*met_cell*2.0e-21);
			  
		    tem = rtTem(losDataStorage[i].cell)*units->temperature;
			  
		    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
		    rC2a =  1.0e-24*exp(-92./tem)*max2(0.0,16.+0.344*sqrt(tem)-47.7/tem);
			  
		    C2_e_cooling=rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
		    C2_a_cooling=rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell;
		    
		    if(tem <= 1.0e3)
		      CII_CNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    else
		      CII_WNM+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
			  
			  
		    CII_CMB_emission= 2.0*C2_abun*cell_HI_density(losDataStorage[i].cell)*met_cell*units->number_density*2.291e-6*constants->k*91.2*exp(-91.2/2.725*abox[min_level]);

		    pressure+=cell_gas_pressure(losDataStorage[i].cell)*units->energy_density/constants->k*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

		    A_V += (cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*5.3e-22*met_cell*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);

		    CII_H2_para = 4.25e-10*pow(tem/100.0,0.124-0.018*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.25*cell_H2_density(cell)*units->number_density;

		    CII_H2_ortho = 5.14e-10*pow(tem/100.0,0.124+0.023*log(tem/100.0))*constants->k*91.2*C2_abun*cell_HI_density(cell)*units->number_density*0.75*cell_H2_density(cell)*units->number_density;

		    tem = 100.0;
		    rC2e = 6.67e-20/pow(tem,0.5)*exp(-92.0/tem);
		    rC2a =  1.0e-24*exp(-92/tem)*max2(0.0,16+0.344*sqrt(tem)-47.7/tem);
		    CII_H2_basic += cell_H2_density(losDataStorage[i].cell)/(cell_gas_density(losDataStorage[i].cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sfr_in_cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    //CII_H2_basic += cell_H2_density(cell)/(cell_gas_density(cell)*constants->XH)*(rC2e*(cell_HII_density(losDataStorage[i].cell)+cell_HeII_density(losDataStorage[i].cell)+2.0*cell_HeIII_density(losDataStorage[i].cell))*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + rC2a*cell_HI_density(losDataStorage[i].cell)*units->number_density*C2_abun*cell_HI_density(losDataStorage[i].cell)*units->number_density*met_cell + 1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*(sf_recipe1_rate(losDataStorage[i].cell)*units->density/units->time*1/constants->Msun*constants->yr)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    sfr_in_cell = 0.0;
		    ipart = cell_particle_list[losDataStorage[i].cell];
		    
		    while(ipart != NULL_PARTICLE)
		      {
			cart_assert( ipart >= 0 && ipart < num_particles );
			if(particle_is_star(ipart))
			      {
				stellar_surface_density+=particle_mass[ipart]*units->mass/constants->Msun*pow(units->length/pow(2.0,cell_level(losDataStorage[i].cell)),-3.0)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
				if((particle_t[ipart]-(double)star_tbirth[ipart])*units->time/constants->yr < 2.0e7)
				  {
				    // in Msun/yr
				    sfr_in_cell += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7);
				    sfr_per_area += (double)star_initial_mass[ipart]*units->mass/constants->Msun*1.0/(2.0e7)/(pow(units->length/constants->pc*pow(2.0,-(double)cell_level(losDataStorage[i].cell)),3.0))*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length/constants->pc;

				  }


			      }
			ipart = particle_list_next[ipart];
		      }

		    PDR_C2cooling = 0.;
		    
		    if((cell_HI_density(losDataStorage[i].cell)+2.0*cell_H2_density(losDataStorage[i].cell))*units->number_density*units->length/pow(2.0,cell_level(losDataStorage[i].cell))*met_cell > 1.0e21)
		      PDR_C2cooling=1.0e-5/(4.0*M_PI)*1.0e21*grain_eff*frac_abs*8.0e27*sfr_in_cell/pow(units->length/pow(2.0,cell_level(losDataStorage[i].cell)),3.0);
		    
		    losHIdensity+=cell_HI_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losH2density+=cell_H2_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losgasdensity+=cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);
		    losmetdensity+=met_cell*constants->Zsun*cell_gas_density(losDataStorage[i].cell)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1);


		    projected_momentum = cell_momentum(losDataStorage[i].cell,0)*v3[0]+cell_momentum(losDataStorage[i].cell,1)*v3[1]+cell_momentum(losDataStorage[i].cell,2)*v3[2];


		    projected_halo_velocity = (double)ngHalos->list[counter].vel[0]*v3[0]+(double)ngHalos->list[counter].vel[1]*v3[1]+(double)ngHalos->list[counter].vel[2]*v3[2];
		    velocity = -Hubble(abox[min_level])*(losDataStorage[i].r2+losDataStorage[i].r1)/2.0*units->length/constants->Mpc + projected_momentum/cell_gas_density(losDataStorage[i].cell)*units->velocity*1e-5-projected_halo_velocity;
		    
		    

		    if(velocity > min_velocity && velocity < max_velocity)
		      {
			optical_depth[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
			optical_depth_all[(int)((velocity-min_velocity)/(max_velocity-min_velocity)*(double)VEL_RES)]+=(C2_e_cooling+C2_a_cooling+CII_HeI+PDR_C2cooling+CII_CMB_emission+CII_H2_para+CII_H2_ortho)*(min2(losDataStorage[i].r2,len)-losDataStorage[i].r1)*units->length*pow(constants->pc,2.0);
		      }

		  }
	      }
	    losHIdensity = log10(losHIdensity*units->number_density*units->length);
	    losH2density = log10(2.0*losH2density*units->number_density*units->length);
	    losmetdensity = log10(losmetdensity/losgasdensity*1.0/constants->Zsun);
	    losgasdensity = log10(losgasdensity*units->number_density*units->length);
	    CII_CNM = CII_CNM*units->length*pow(constants->pc,2.0);
	    CII_WNM = CII_WNM*units->length*pow(constants->pc,2.0);
	    PDR_C2cooling = PDR_C2cooling;
	    CII_CMB_emission = CII_CMB_emission*units->length*pow(constants->pc,2.0);
	    pressure = pressure*units->length;
	    A_V = A_V*units->length;
	    CII_H2_basic = CII_H2_basic*units->length*pow(constants->pc,2.0);
	    total_sfr += sfr_per_area*pow(units->length/pow(2.0,(double)NMAX)*1.0/constants->pc,2.0);

	    for(i=0;i<VEL_RES;i++)
	      {
		flux_profile[i]+=optical_depth[i]*pow(units->length/pow(2.0,(double)NMAX)*1.0/constants->pc,2.0)/(1.900539e12*(max_velocity-min_velocity)/((double)VEL_RES*2.9979e5))*pow(1.322e18,-2.0)*1.0e3*1.0/(4.0*M_PI);
		flux_profile_all[i]+=optical_depth_all[i]*pow(units->length/pow(2.0,(double)NMAX)*1.0/constants->pc,2.0)/(1.900539e12*(max_velocity-min_velocity)/((double)VEL_RES*2.9979e5))*pow(1.322e18,-2.0)*1.0e3*1.0/(4.0*M_PI);
		//This has the units of Jy but we must multiply by 1/((1+z)*chi^2) where chi is the dimensionless comoving distance
	      }


	    //fprintf(dataout,"%f %f %f %f %f %f %e %e %e %e %e %e %e %e %e %e\n",v1_counter/pow(2.0,NMAX)*units->length/constants->kpc,v2_counter/pow(2.0,NMAX)*units->length/constants->kpc,losHIdensity,losH2density,losgasdensity,losmetdensity,CII_CNM,CII_WNM,PDR_C2cooling,CII_CMB_emission, pressure,A_V,CII_H2_basic,CII_H2_para+CII_H2_ortho,stellar_surface_density,sfr_per_area);


	    bufferout=v1_counter/pow(2.0,NMAX)*units->length/constants->kpc;
	    fwrite(&bufferout,sizeof(double),1,dataout);
	    fwrite(&bufferout,sizeof(double),1,dataout_all);
	    	    
	    bufferout=v2_counter/pow(2.0,NMAX)*units->length/constants->kpc;
	    fwrite(&bufferout,sizeof(double),1,dataout);
	    fwrite(&bufferout,sizeof(double),1,dataout_all);
	    
	    fwrite(optical_depth,sizeof(double),VEL_RES,dataout);
	    fwrite(optical_depth_all,sizeof(double),VEL_RES,dataout_all);

	    //fprintf(dataout,"%f %f\n",v1_counter/pow(2.0,NMAX)*units->length/constants->kpc,v2_counter/pow(2.0,NMAX)*units->length/constants->kpc);
	    //for(i=0;i<VEL_RES;i++)
	    //  fprintf(dataout,"%f %e\n",min_velocity + ((double)i+0.5)/(double)VEL_RES*(max_velocity - min_velocity),optical_depth[i]);
	  }
    
      fclose(dataout);
      fclose(dataout_all);
      free(losDataStorage);
      printf("SFR for halo %d = %e Mstar/yr\n",counter,total_sfr);

      sprintf(fname,"CII/by_halo_map/N%d_%d_flux_profile_halo_%d.txt",NMIN,NMAX,counter);
      dataout = fopen(ngOutputFile(fname),"w");
      fprintf(dataout,"#This has the units of Jy but we must multiply \n");
      fprintf(dataout,"#by 1/((1+z)*chi^2) where chi is the dimensionless comoving distance \n");
      for(i=0;i<VEL_RES;i++)
	fprintf(dataout,"%f %e\n",min_velocity + ((double)i+0.5)/(double)VEL_RES*(max_velocity-min_velocity),flux_profile[i]);
      fclose(dataout);

      sprintf(fname,"CII/by_halo_map/N%d_%d_flux_profile_halo_%d_all.txt",NMIN,NMAX,counter);
      dataout = fopen(ngOutputFile(fname),"w");
      fprintf(dataout,"#This has the units of Jy but we must multiply \n");
      fprintf(dataout,"#by 1/((1+z)*chi^2) where chi is the dimensionless comoving distance \n");
      for(i=0;i<VEL_RES;i++)
	fprintf(dataout,"%f %e\n",min_velocity + ((double)i+0.5)/(double)VEL_RES*(max_velocity-min_velocity),flux_profile_all[i]);
      fclose(dataout);


    }

  free(optical_depth);
  free(optical_depth_all);
  free(flux_profile);
  free(flux_profile_all);

  cart_free(ngHalos->list);
  cart_free(ngHalos);


  printf("Exiting CII_map_with_galaxy_face_on_velocity_profiles(), a= %.4f\n",abox[min_level]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////



//Finally we have all our defined functions.

double variablequantity(int cell, int value)
{
    switch(value)
    {
    case 0: return cell_gas_density(cell)*units->number_density;
      //case 0: return cell_gas_density(cell)*units->density*1/(constants->XH*constants->mH + constants->xHe*constants->mHe);
    case 1: return cell_HI_density(cell)*units->number_density;
      //case 1: return cell_HI_density(cell)*units->density*1/constants->mH;
    case 2: return 1.0;
    case 3: return cell_gas_metal_density(cell)/cell_gas_density(cell);
    case 4: return cell_HI_density(cell)*units->number_density + 2.0*cell_H2_density(cell)*units->number_density;
    case 5: return 2*cell_H2_density(cell)*units->number_density;
    case 6: return cell_gas_density(cell)*constants->XH*units->number_density;
      //case 5: return 2*cell_H2_density(cell)*units->density*1/constants->mH;
    default: printf("#Invalid selection in variablequantity. You must restart the analysis!\n");
      return 0.0;
      break;
    }
}
  
char* variablequantityname(int value)
{
  switch(value)
    {
    case 0: return "cell_gas_density";
    case 1: return "cell_HI_density";
    case 2: return "DLA_length";
    case 4: return "cell_HI+2H2_density";
    case 5: return "cell_H2_density";
    case 6: return "cell_H_density";
    default: printf("#Invalid selection in variablequantityname. You must restart the analysis!\n");
      return "blah";
      break;
    }
}

char* variableaxesname(int value)
{
  switch(value) 
    {
    case 0: return "xy";
    case 1: return "xz";
    case 2: return "yz";
    default: printf("#Invalid selection in variableaxesname. You must restart the analysis!\n");
      return "blah";
      break;
    }
}

int variableaxisone(int value)
{
  switch(value)
    {
    case 0: return 0;
    case 1: return 0;
    case 2: return 1;
    default: printf("#Invalid selection in variableaxisone. You must restart the analysis!\n");
      return 42;
      break;
    }
}

int variableaxistwo(int value)
{
  switch(value)
    { 
    case 0: return 1;
    case 1: return 2;
    case 2: return 2;
    default: printf("#Invalid selection in variableaxistwo. You must restart the analysis!\n");
      return 42;
      break;
    } 
}

int projected_axis(int value)
{
  switch(value)
    {
    case 0: return 2;
    case 1: return 1;
    case 2: return 0;
    default: printf("#Invalid selection in projected_axis. You must restart the analysis!\n");
      return 42;
      break;
    }
}

double min3(double one, double two, double three)
{
  if(one >= two)
    {
      if(two >= three)
	return three;
      else
	return two;
    }
  else
    {
      if(three >= one)
	return one;
      else
	return three;
    }
}

double min2(double one, double two)
{
  if(one >=two)
    return two;
  else return one;
}

double max2(double one, double two)
{
  if(one >= two)
    return one;
  else return two;
}

double max3(double one, double two, double three)
{
  if(one >= two)
    {
      if(one >= three)
	return one;
      else
	return three;
    }
  else
    {
      if(two >= three)
	return two;
      else
	return three;
    }
}

double mod_map(double x)
{
  if(x>0)
    return min3(fabs(x),fabs(x+RGSIZE),fabs(x-RGSIZE));
  else
    return -min3(fabs(x),fabs(x+RGSIZE),fabs(x-RGSIZE));
}

double mod_DE(double x, double y)
{
  if(x<0)
    {
      while(x<0)
	x+=y;
      return x;
    }
  if(x >= y)
    {
      while(x>=y)
	x-=y;
      return x;
    }
  return x;

}

double distance_cartesian(double position1[3],double position2[3])
{
  return max3(sqrt(min3(pow(position1[0]-position2[0],2.0),pow(position1[0]-position2[0]-RGSIZE,2.0),pow(position1[0]-position2[0]+RGSIZE,2.0))),sqrt(min3(pow(position1[1]-position2[1],2.0),pow(position1[1]-position2[1]-RGSIZE,2.0),pow(position1[1]-position2[1]+RGSIZE,2.0))),sqrt(min3(pow(position1[2]-position2[2],2.0),pow(position1[2]-position2[2]-RGSIZE,2.0),pow(position1[2]-position2[2]+RGSIZE,2.0))));
}

double distance(double position1[3],double position2[3])
{
  return sqrt(min3(pow(position1[0]-position2[0],2.0),pow(position1[0]-position2[0]-RGSIZE,2.0),pow(position1[0]-position2[0]+RGSIZE,2.0))+min3(pow(position1[1]-position2[1],2.0),pow(position1[1]-position2[1]-RGSIZE,2.0),pow(position1[1]-position2[1]+RGSIZE,2.0))+min3(pow(position1[2]-position2[2],2.0),pow(position1[2]-position2[2]-RGSIZE,2.0),pow(position1[2]-position2[2]+RGSIZE,2.0)));
}

double twoDdistance(double position1[2],double position2[2])
{
  return sqrt(min3(pow(position1[0]-position2[0],2.0),pow(position1[0]-position2[0]-RGSIZE,2.0),pow(position1[0]-position2[0]+RGSIZE,2.0))+min3(pow(position1[1]-position2[1],2.0),pow(position1[1]-position2[1]-RGSIZE,2.0),pow(position1[1]-position2[1]+RGSIZE,2.0)));
}


void cross_product(double *x, double *y, double *result)
{
  result[0]=x[1]*y[2]-x[2]*y[1];
  result[1]=-x[0]*y[2]+x[2]*y[0];
  result[2]=x[0]*y[1]-x[1]*y[0];
}

double projected_distance(double vector1[3], double vector2[3])
{
  return fabs(vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2]);  
}


void find_primary_halos(int *primary_halo_id,int *primary_halo_number)
{
  halo_list *ngHalos = NULL;
  char str[999];
  
  sprintf(str,"DAT_HC/HC/hlist_%6.4f.dat",auni[min_level]);
  ngHalos = load_halo_finder_catalog(str,50,0.0,0.0,0.0,100000);
  
  int i,j;

  //This gives us our primary halos

  int halo_is_primary = 1;

  *primary_halo_number=0;

  for(i =0;i<ngHalos->num_halos;i++)
    {
      for(j=0;j<ngHalos->num_halos && halo_is_primary;j++)
	{
	  if(ngHalos->list[i].mvir < ngHalos->list[j].mvir && distance(ngHalos->list[i].pos,ngHalos->list[j].pos) < ngHalos->list[j].rvir && j!=i)
	    {
	      halo_is_primary = 0;
	    }
	}
      if(halo_is_primary)
	{
	  primary_halo_id[*primary_halo_number]=i;
	  *primary_halo_number=*primary_halo_number+1;
	}

      halo_is_primary = 1;
    }

  cart_free(ngHalos->list);
  cart_free(ngHalos);
}

int number_of_cells_in_annulus(double r_in, double r_out)
{
  int x,y;
  int num_elts=0;

  if(r_out == 0.0)
    {
      printf("r_out = 0.0 doesn't make sense\n");
      exit(1);
    }
  
    for(x=0;x<=(int)r_out+1;x++)
    for(y=0;y<=(int)r_out+1;y++)
      if(sqrt(pow((double)x,2.0)+pow((double)y,2.0)) >= r_in && sqrt(pow((double)x,2.0)+pow((double)y,2.0)) < r_out)
	{
	  if(x==0 || y==0)
	    num_elts=num_elts+2;
	  else
	    num_elts=num_elts+4;
	}
  if(r_in == 0.0)
    return num_elts - 1;
  else
    return num_elts;
}

int WorkerBee(int id, int cell, double r1, double r2, losBuffer *data)
{

  losdatatype *datatostore;
  
  datatostore = (*data).Data;
  
  (*datatostore).r1=r1;
  (*datatostore).r2=r2;
  
  (*datatostore).cell=cell;

  //printf("WorkerBee says: %d %d %f %f\n",cell,cell_level(cell),r1,r2);

  (*data).Data=(*data).Data+(*data).Size;

  //printf("Exiting workerbee!\n");

   return 0;
}

int double_cmp(const void *a, const void *b)
{
  double A, B;

  A = *(double *)a;
  B = *(double *)b;
  
  if(A < B)
    return -1;
  else if(A > B)
    return +1;
  return 0;
}

void read_cloudy_grid(){
  double xcgZ,xcgrho,xcgT;
  float x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
  FILE *fp;
  int ind, npcldf;
  int ix1,ix2,ix3;
  int j=0, i;

  double rfact = 1000.0 *units->length * abox[min_level] / cosmology->h;

  double convnah, hproxy,zproxy;
  char a[1000];
  
  int iel;
    for(iel = 0 ; iel < MAX_ELEMENT ; iel++){
      cart_debug("lambda0[%d] = %e",iel,lam0_els[iel]);
      cart_assert(lam0_els[iel] > 100e-8 );
    }

  if(!(fp=fopen("/home/derkal/trunk/utils/analysis_dke/cloudy/all.z-5to1.h-6to3.Te2toe7","r")))
    cart_error("couldn't open cloudy grid file") ;


  //==============================READ IN CLOUDY GRID =====================...
  npcldf=0;
  while(!feof(fp))
    {
      fgets(a,1000,fp);
      npcldf++;
    }
  rewind(fp);
  npcldf--;

  cart_debug("extra els %d npcld %d %d %d %d",
	  npcldf-npcld,npcld, ngz,ngd,ngt );

  j=0;
  for (i=0;i<npcldf;i++)  {
    fscanf(fp,"%f %f %e %f %e %e %e %e %e %e %e",
	   &x1,&x2,&x3,&x4,&x5,&x6,&x7,&x8,&x9,&x10,&x11); 
    //to invert these x1=cloudy_step_Z*ix1+cgZlo
    ix1   =(x1-cgZlo+eps)/cloudy_step_Z;
    ix2   =(x2-cgrholo+eps)/cloudy_step_rho;
    ix3   =(log10((double)x3)-cgTlo+eps)/cloudy_step_T;//x3 should change fastest so d3cld() by 1, so efficient mem
    ind = d3cld(ix1,ix2,ix3);
    

    //-------------------
    //CHECK CLOUDY GRID IS WHAT YOU EXPECT AND EPS IS OK
    //include only values from chosen grid == <1% accurate steps
    if( fabs((x1-cgZlo)/cloudy_step_Z-ix1)>.01 || 
	fabs((x2-cgrholo)/cloudy_step_rho-ix2)>.01 ||
	fabs((log10((double)x3)-cgTlo)/cloudy_step_T-ix3)>.01 ||
	ix1>=ngz||ix2>=ngd||ix3>=ngt||
	ind>=npcld){


      
      cart_debug("Doesn't fit my beautiful grid:"
		 " %f %f %e %d %d %d %d %d %d ; %f %f %f  ; %d %d", 
		 x1,x2, x3, ix1,ngz,ix2,ngd,ix3,ngt, 
		 ((x1-cgZlo)/cloudy_step_Z-ix1),((x2-cgrholo)/cloudy_step_rho-ix2),
		 ((log10((double)x3)-cgTlo)/cloudy_step_T-ix3),
		 ind,npcld );
      
      
      
      j=j+1;
    } else{
      nden_cloudy[EL_O6][ind]  = log10( (double)x5 )  - x4; //column #/cm**2 -> #/cm**3
      nden_cloudy[EL_C4][ind]  = log10( (double)x6 )  - x4; //column #/cm**2 -> code
      nden_cloudy[EL_C3][ind]  = log10( (double)x7 )  - x4; //column #/cm**2 -> code
      nden_cloudy[EL_N5][ind]  = log10( (double)x8 )  - x4; //column #/cm**2 -> code
      nden_cloudy[EL_Si4][ind] = log10( (double)x9  ) - x4; //column #/cm**2 -> code
      nden_cloudy[EL_Si3][ind] = log10( (double)x10 ) - x4; //column #/cm**2 -> code
      nden_cloudy[EL_Mg2][ind] = log10( (double)x11 ) - x4; //column #/cm**2 -> code
    }
    
    cloudy_grid_dr = pow(10,x4)/constants->kpc/rfact; //cm->code units
    //nden = 1/cm2,  1/(cm2* 1/(kpc*rfact) ) = kpc*rfac / cm2

    //-------------------
    if(i==10){cart_debug("cloudy test: %f %f %e %f %e %e %e %e %e %e %e",
		     x1,x2,x3,x4,nden_cloudy[EL_O6][ind],x6,x7,x8,x9,x10,x11 );}

  }
  if(j != npcldf-npcld){cart_error("%d should equal npcldf-npcld %d", j, npcldf-npcld);}
  fclose(fp);

  //=========================DONE READ IN CLOUDY GRID=====================...
}


double cloudy_column(double zmet, double rho_nH, double gtemp, double element[npcld]){ 
  double col;
  double xcgZ, xcgrho, xcgT;
  double model_extrap=0;


  if(      zmet >= cgZhi){
    model_extrap += zmet-cgZhi ;
    zmet = cgZhi-eps; 
  }else if(zmet <= cgZlo){
    model_extrap += zmet-cgZlo ;
    zmet = cgZlo+eps; 
  }

  if(      rho_nH >= cgrhohi){
   // model_extrap += rho_nH - cgrhohi; //not obviously, but potentially bad
    rho_nH = cgrhohi-eps; 
  }else if(rho_nH <= cgrholo){
    model_extrap += rho_nH - cgrholo;
    rho_nH = cgrholo+eps; 
  }

  if(      gtemp >= cgThi){
    gtemp = cgThi-eps; 
  }else if(gtemp <= cgTlo){
    gtemp = cgTlo+eps; 
  }

  //---------discretize to cloudy grid.---------
  xcgZ=(zmet-cgZlo)/cloudy_step_Z;
  xcgrho=(rho_nH-cgrholo)/cloudy_step_rho; 
  xcgT=(gtemp-cgTlo)/cloudy_step_T;
  
  cloudy_grid_check(xcgZ,xcgrho,xcgT);
  
  col= cic3d(xcgZ,xcgrho,xcgT, element)+ model_extrap; 
  
  return col;
}

void cloudy_grid_check(double xcgZ,double xcgrho,double xcgT){
  if(xcgZ >= ngz - BAD_CGRID_Z || xcgZ < 0){cart_error("bad metal: %f %d ",xcgZ,ngz - BAD_CGRID_Z );}
  if(xcgrho >= ngd || xcgrho < 0){cart_error("bad density: %f %d",xcgrho,ngd);}
  if(xcgT >= ngt || xcgT < 0){cart_error("bad temperature: %f %d",xcgT,ngt);}
  
}


double cic3d(double x1, double x2, double x3, double *data)
{
  int xc,yc,zc, xcn,ycn,zcn;
  double dx,dy,dz, tx,ty,tz;
  double val;

    xc=x1;
    yc=x2;
    zc=x3;

    xcn=(xc+1);
    ycn=(yc+1);
    zcn=(zc+1);
    
    dx=x1-xc;
    dy=x2-yc;
    dz=x3-zc;
    
    tx=1-dx;
    ty=1-dy;
    tz=1-dz;
    
    val=
      data[d3cld(xc ,yc ,zc )]*tx*ty*tz+
      data[d3cld(xcn,yc ,zc )]*dx*ty*tz+
      data[d3cld(xc ,ycn,zc )]*tx*dy*tz+
      data[d3cld(xc ,yc ,zcn)]*tx*ty*dz+
      data[d3cld(xcn,ycn,zc )]*dx*dy*tz+
      data[d3cld(xcn,yc ,zcn)]*dx*ty*dz+
      data[d3cld(xc ,ycn,zcn)]*tx*dy*dz+
      data[d3cld(xcn,ycn,zcn)]*dx*dy*dz;
    
    return val;
}

double velocity_width_from_optical_depth(double *optical_depth, int RES)
{
  int i;
  double interpolated_index_min;
  double interpolated_index_max;
  
  for(i=0;i<RES-1;i++)
    optical_depth[i+1]+=optical_depth[i];
  
  i=0;
  
  while(optical_depth[i]/optical_depth[RES-1] < 0.05)
    i++;

  interpolated_index_min = (double)i-1.0+(0.05*optical_depth[RES-1]-optical_depth[i-1])/(optical_depth[i]-optical_depth[i-1]);
    
  i=RES-1;
  
  while(optical_depth[i]/optical_depth[RES-1] > 0.95)
    i--;

  interpolated_index_max = (double)i+(0.95*optical_depth[RES-1]-optical_depth[i])/(optical_depth[i+1]-optical_depth[i]);

  return interpolated_index_max-interpolated_index_min;
}


double jeans_length(int cell)
{
  int i, j, nb[num_neighbors];
  float s, d, q, dv2, cs2, ljeans;
  float vc[nDim];

  if(cell_gas_density(cell) < 1.0e-30) return 0.0;

  /*
  //  Sound speed squared
  */
  cs2 = cell_gas_gamma(cell)*MAX((cell_gas_gamma(cell)-1.0)*cell_gas_internal_energy(cell),0.0)/cell_gas_density(cell);

  /*
  //  Local velocity dispersion squared
  */
  cell_all_neighbors(cell,nb);
  for(i=0; i<nDim; i++)
    {
      vc[i] = cell_var(cell,HVAR_MOMENTUM+i)/cell_gas_density(cell);
    }

  s = 0.0;
  for(j=0; j<num_neighbors; j++) if(cell_gas_density(nb[j]) > 0.0)
    {
      for(i=0; i<nDim; i++)
	{
	  d = cell_var(nb[j],HVAR_MOMENTUM+i)/cell_gas_density(nb[j]) - vc[i];
	  s += d*d;
	}
    }
  dv2 = s/num_neighbors;

  ljeans = sqrt(M_PI*(cs2+dv2)/((constants->G*units->mass*pow(units->time,2.0)/pow(units->length,3.0))*cell_gas_density(cell)));

  return ljeans;
}

double jeans_length_no_turb(int cell)
{
  float cs2, ljeans;
  
  if(cell_gas_density(cell) < 1.0e-30) return 0.0;

  /*
  //  Sound speed squared
  */
  cs2 = cell_gas_gamma(cell)*MAX((cell_gas_gamma(cell)-1.0)*cell_gas_internal_energy(cell),0.0)/cell_gas_density(cell);

  ljeans = sqrt(M_PI*(cs2)/((constants->G*units->mass*pow(units->time,2.0)/pow(units->length,3.0))*cell_gas_density(cell)));

  return ljeans;
}
