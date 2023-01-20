F90	  =  mpif90
F90LINKER = $(F90)
LIB_DIR   = 
INC_DIR   = 
OPTFLAGS  = -Oss -noparallel
FFLAGS = $(OPTFLAGS)
FLIBS = 
F90LFLAGS=
#
TARGET = ../run/sol
default: $(TARGET)
OBJS =\
pfem_util.o \
solver_SR_3.o solver_CG_3.o \
solver33.o test1.o util.o pfem_init.o input_cntl.o input_grid.o define_file_name.o\
mat_con0.o mat_con1.o mat_ass_main.o mat_ass_bc.o pfem_finalize.o output_ucd.o recover_stress.o

$(TARGET):  $(OBJS)
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(OBJS) $(F90LFLAGS)
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
.f.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f
.f90.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f90
.SUFFIXES: .f90 .f
