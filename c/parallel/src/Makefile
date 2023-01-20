CC	 =mpicc
LIB_DIR=
INC_DIR=  
OPTFLAGS= -Os -noparallel
LIBS = 
LFLAGS= -lm
#
TARGET = ../run/sol
default: $(TARGET)
OBJS =\
	test1_p.o\
	pfem_init.o\
	input_cntl_p.o\
	input_grid_p.o\
	define_file_name.o\
	mat_con0.o\
	mat_con1.o\
	mat_ass_main.o\
	mat_ass_bc.o\
	solve33_p.o\
	solver_CG_3_p.o\
	solver_SR_3.o\
	recover_stress_p.o\
	output_ucd_p.o\
	pfem_finalize.o\
	allocate.o\
	util.o

$(TARGET):  $(OBJS)
	$(CC) $(OPTFLAGS) -o $@ $(OBJS) $(LFLAGS)
.c.o:
	$(CC) $(OPTFLAGS) -c  $*.c
clean:
	/bin/rm -f *.o *~
