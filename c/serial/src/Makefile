CC 	    = cc
CCFLAGS     = -Os -noparallel

LD          = $(CC)
LDFLAGS     = 

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) -c $< $(CCFLAGS)

TARGET = ../run/sol

all: $(TARGET)

OBJS = \
	solver_CG_3.o \
	solve33.o \
	test1.o \
	util.o \
	input_cntl.o \
	input_grid.o \
	mat_con0.o \
	mat_con1.o \
	mat_ass_main.o \
	mat_ass_bc.o \
	output_ucd.o \
	recover_stress.o \
	allocate.o

$(TARGET): $(OBJS)
	$(LD) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod

