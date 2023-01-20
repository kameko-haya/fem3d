F90       =  gfortran
F90LINKER = $(F90)
OPTFLAGS  = -O2
FFLAGS = $(OPTFLAGS)
F90LFLAGS=
#
TARGET = ../run/sol
default: $(TARGET)
OBJS =\
m_util.o\
m_solver.o\
input_ctrl.o \
mat_con.o\
mat_ass.o\
main.o 

$(TARGET):  $(OBJS)
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(OBJS) $(F90LFLAGS)
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
.f.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f
.f90.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f90
.SUFFIXES: .f90 .f
