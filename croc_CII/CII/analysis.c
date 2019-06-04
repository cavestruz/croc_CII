#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef RADIATIVE_TRANSFER
#include "analysis_dke.h"
#else
#include "analysis_dke_noRT.h"
#endif // RADIATIVE_TRANSFER

char *ngOutputDir = NULL;

const char *ngOutputFile(const char *filepath)
{
  int ret;
  static char str[999];
  char *name;

  if(ngOutputDir == NULL)
    {
      strcpy(str,filepath);
    }
  else
    {
      strcpy(str,ngOutputDir);
      strcat(str,"/");
      strcat(str,filepath);
    }

  name = strrchr(str,'/');
  if(name != NULL)
    {
      *name = 0;
      /*
      //  Create output directory
      */
      ret = system_mkdir(str);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	  return NULL;
	}
      *name = '/';
    }

  return str;
}


int main_analysis(int argc, char **argv)
{
  printf("Begun the main analysis has!\n");

  int ret;
  char str[999];

  //ngLoadHalos(NULL);

  sprintf(str,"a%6.4f",auni[min_level]);

  ret = system_mkdir(str);
  if(ret != 0)
    {
      cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
      return -1;
    }
  ngOutputDir = str;

#ifdef RADIATIVE_TRANSFER

  CII_map_with_galaxy_face_on_velocity_profiles();

#endif //RADIATIVE_TRANSFER
  return 0;
}
