#define RT_OUTPUT

#define RT_TRANSFER
#define RT_TRANSFER_METHOD RT_METHOD_OTVET
#define RT_UV
/*#define RT_UV_OLDSTYLE_3x1*/ /*for runs with the pre-ng branch RT_UV*/

#define RT_CHEMISTRY 
#define RT_LWBANDS
#define RT_DUST_EVOLUTION
#define RT_SIGNALSPEED_TO_C 
#define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_HAARDT_MADAU

#define RT_PARALLEL_NUM_OPENMP_BUFFERS 8
#define RT_PARALLEL_USE_MPI

/*#define BLASTWAVE_FEEDBACK*/ /* Needed for BLASTWAVE obviously */

/* FIXED ISM */
/*#define RT_FIXED_ISM*/
