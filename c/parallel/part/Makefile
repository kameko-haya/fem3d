F77	= ifort 
F90	= ifort 
FLINKER   = $(F77)
F90LINKER = $(F90)
FLIB_PATH =
INC_DIR   =  
OPTFLAGS = -Oss -noparallel
FFLAGS = $(OPTFLAGS)
FLIBS  = $(FLIB_PATH) $(LIB_LIST)
#
#TARGET = hpcmw-solver-scalar-test
TARGET = ../mesh/part
default: $(TARGET)
OBJS =\
geofem_util.o \
partitioner.o analyzer.o geofem_tiger.o \
calc_edgcut.o cre_local_data.o define_file_name.o \
interface_nodes.o list.o create_metis_input.o\
neib_pe.o ordering.o paraset.o partition.o prep_next_partition.o\
reduce_edgcut.o rhoset.o proc_local.o local_data.o\
double_numbering.o output_ucd.o

$(TARGET):  $(OBJS)
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(OBJS) $(FLIBS)
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
.f.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f
.SUFFIXES: .f
