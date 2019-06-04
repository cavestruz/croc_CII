/* Primary Settings */

#define HYDRO
#define GRAVITY
#define COOLING
#define STAR_FORMATION
#define COSMOLOGY
#define PARTICLES
#define REFINEMENT
#define RADIATIVE_TRANSFER
/* #define DUST_EVOLUTION */

/* Star formation settings */

#define ENRICHMENT
#define ENRICHMENT_SNIa
#define SF_RECIPE			<gk10-full> /* gk10-full if RT on, gk10-lite*/
#define SF_FEEDBACK                     <hart>/*<popM-thermal>*//*<hart>*/


#define DEBUG_MEMORY_USE
#define OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID /* Used in all of my runs */

/* For analysis including CO */

/*#define CO*/


/* for the 25Mpc runs on a huge*/

#define num_root_grid_refinements	8
#define num_refinement_levels		12
#define num_octs			60000000
#define num_particles		        320000000
#define num_star_particles              120000000
#define PARTICLE_HEADER_MAGIC           (0.0f)


/* for the 6Mpc RUN*/
/*
#define num_root_grid_refinements	6
#define num_refinement_levels		9
#define num_octs			6000000
#define num_particles		        9000000
#define num_star_particles              9000000
*/

/* For 6Mpc COSMO.BGHM */
/*
#define num_root_grid_refinements       6
#define num_refinement_levels           9
#define num_octs                        5000000
#define num_particles                   16000000
#define num_star_particles              16000000
*/
