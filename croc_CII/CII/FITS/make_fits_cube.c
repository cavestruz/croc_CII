#include <string.h>
#include <stdlib.h>
#include <fitsio.h>
#include <math.h>

// these must match the values used in the simulation

#define VEL_RES 100
#define min_velocity -1000.0
#define max_velocity 1000.0

#define AXES_LENGTH 300
//382 for 6mpc at a=0.2502
//378 for 25mpc at a=0.2501
//300 for 25mpc at a=0.3200
//each cell has 102*8 bytes of memory so you can determine the length from the file size


#define halo_counter 0

#define z 2.0
#define a_act 0.3200
#define x 1.21 // d_L = (1+z)*c/H_0*x
#define l_cell 272.48 //pc 261.58 for 6mpc box
#define nu_158 1897.0e9



//gcc make_fits.c -I /usr/include/cfitsio/ -lcfitsio

//gcc -c make_fits_cube.c -I /home/derkal/Desktop/Downloaded_Files/cfitsio

//gcc -o make_fits_cube make_fits_cube.o -L. -lcfitsio -lm -L /home/derkal/Desktop/Downloaded_Files/cfitsio/lib/

int main()
{
  FILE *datain;
  char input_string[999];
  char fname[99];
  
  double pos[2];
  double NHI,NH2,Ngas,Zmet,CII_CNM,CII_WNM,CII_PDR,CII_CMB,Pressure,AVmag;
  int i,j,k,counter;

  fitsfile *fptr;
  int status=0;
  long naxes[4];
  long fpixel[4];
  long nelements;
  double *value=malloc(sizeof(double));
  double *optical_depth = malloc(sizeof(double)*VEL_RES);
  naxes[0]=AXES_LENGTH;
  naxes[1]=AXES_LENGTH;
  naxes[3]=1; //This is the polarization axis
  naxes[2]=VEL_RES;

  printf("Hello\n");

  sprintf(fname,"../../../RUN.3_25Mpc/a0.3200/by_halo_map/N2_9_face_on_velocity_map_halo_%d.bin",halo_counter);

  datain = fopen(fname,"rb");
  if(datain == NULL)
    {
      printf("Can't open the galaxy map\n");
      return -1;
    }

  sprintf(fname,"!halo_%d/CII_emission_halo_%d_cube.fits",halo_counter,halo_counter);

  fits_create_file(&fptr,fname,&status);
  if(status)
    {
      fits_report_error(stderr,status);
    }
  fits_create_img(fptr,-64,4,naxes,&status); //for double
  if(status)
    {
      fits_report_error(stderr,status);
    }

  fits_write_record(fptr,"BSCALE  =    1.00000000000E+00",&status);
  fits_write_record(fptr,"BZERO   =    0.00000000000E+00",&status);

  fits_write_record(fptr,"CDELT3  =    2.00000000000E+04",&status);
  fits_write_record(fptr,"CRPIX3  =    1.00000000000E+00",&status);
  fits_write_record(fptr,"CRVAL3  =    0.00000000000E+00",&status);
  fits_write_record(fptr,"CTYPE3  = 'VELO-LSR'",&status);
  
  if(status)
    {
      fits_report_error(stderr,status);
    }


  /*for(i=0;i<1000;i++)
    for(j=0;j<1000;j++)
      {
	fpixel[0]=i+1;
	fpixel[1]=j+1;
	*value = cos(5.0*(double)i+3.0*(double)j);
	fits_write_pix(fptr,TDOUBLE,fpixel,1,value,&status);
	}*/


  i=j=1;
  counter=0;
  fpixel[3]=1;

  printf("Opened the files\n");

  while(!feof(datain))
    {

      fread(pos,sizeof(double),2,datain);
      fread(optical_depth,sizeof(double),VEL_RES,datain);

      if(i>AXES_LENGTH)
	{
	  i=1;
	  j++;
	  counter++;
	  if(counter > (int)((double)AXES_LENGTH/10.0))
	    {
	      counter=0;
	      printf("Done with %.0f percent\n",(double)j/(double)AXES_LENGTH*100.0);
	    }
	}
      fpixel[0]=i;
      fpixel[1]=j;
      for(k=1;k<=VEL_RES;k++)
	{
	  fpixel[2]=k;
	  *value = 1.0e23*(1.0+z)/(4.0*M_PI*pow((1.0+z)*1.324e28*x,2.0))*optical_depth[k-1]*pow(l_cell*a_act,2.0)/(nu_158*(max_velocity-min_velocity)/(3.0e5*(double)VEL_RES));
	  fits_write_pix(fptr,TDOUBLE,fpixel,1,value,&status);
	}
      
      i++;
    }
  
  fclose(datain);
  
  fits_close_file(fptr,&status);
  if(status)
	{
	  fits_report_error(stderr,status);
	}
      
  free(value);
  free(optical_depth);

  return 0;
}
