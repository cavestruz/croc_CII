#From /lustre/sleitner/origin/analysis

TOP_DIR = /home/derkal/trunk_ng#
#
#  The rest of this file does not need to be changed
#

#  List your plugin file(s) here (if needed).
#  Adding /extra will cause all code from src/extra linked in as well.
#
PLUGIN = 


SRC_DIR = $(TOP_DIR)/src

compile: analysis

include $(TOP_DIR)/Make.config
#This one is apparently supposed to convert to the old file format
include $(TOP_DIR)/utils/universal_file_converter/Make.include 

analysis_dke.o: analysis_dke.c analysis_dke.h defs.h rt_defs.h $(ANL_OBJS)
	$(CC) $(CFLAGS) -I. $(INCLUDES) -I$(SRC_DIR) -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(SRC_DIR)/extra -c analysis_dke.c


analysis: .obj/.phony analysis.c analysis_dke.o analysis_dke.h $(ANL_OBJS)  defs.h rt_defs.h analysis.h
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(SRC_DIR)/extra $(INCLUDES) $(ANL_OBJS) analysis_dke.o analysis.c -o analysis $(LIBS) $(EXTRALIBS)
#	which mpirun

analysis_dke_noRT.o: analysis_dke_noRT.c analysis_dke_noRT.h defs.h $(ANL_OBJS)
	$(CC) $(CFLAGS) -I. $(INCLUDES) -I$(SRC_DIR) -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(SRC_DIR)/extra -c analysis_dke_noRT.c

analysis_noRT: .obj/.phony analysis.c analysis_dke_noRT.o analysis_dke_noRT.h $(ANL_OBJS)  defs.h analysis.h
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(SRC_DIR)/extra $(INCLUDES) $(ANL_OBJS) analysis_dke_noRT.o analysis.c -o analysis_noRT $(LIBS) $(EXTRALIBS)
#	which mpirun

CLEAN += *.o .obj analysis analysis_noRT
