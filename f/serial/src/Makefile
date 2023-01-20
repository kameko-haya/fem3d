F90	  =  f90
F90LINKER = $(F90)
OPTFLAGS  = -Oss -noparallel
FFLAGS = $(OPTFLAGS)
F90LFLAGS=
#
TARGET = ../run/sol
default: $(TARGET)
OBJS =\
pfem_util.o \
solver_CG_3.o \
solver33.o test1.o util.o input_cntl.o input_grid.o \
mat_con0.o mat_con1.o mat_ass_main.o mat_ass_bc.o output_ucd.o recover_stress.o

$(TARGET):  $(OBJS)
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(OBJS) $(F90LFLAGS)
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
.f.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f
.f90.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f90
.SUFFIXES: .f90 .f
