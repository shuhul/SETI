# Optimize
OPTIM =   -fbacktrace

# The system defines used in executing make
#FF = f77 -C
FF = gfortran  -g -fcheck=all 

# Glueing all the things together
FFLAGS   = $(OPTIM) 
#-Wall
LFLAGS   = $(OPTIM) 
#-Wall


#
# Names
#
COMMONS        = common_grid.h common_radialrt.h common_disk.h configure.h
SRCDISK        = diskevolfast.F 
OBJDISK        = diskevolfast.o
SRCRADIALRT    = radialrt2.F 
OBJRADIALRT    = radialrt2.o
SRCNRECIP      = nrecip.F 
OBJNRECIP      = nrecip.o
SRCDUST        = dusta.f90 
OBJDUST        = dusta.o
SRCDEVLV       = dustcollision.F
OBJDEVLV       = dustcollision.o
SRCALP         = alpha.f
OBJALP         = alpha.o

OBJ            = $(OBJDISK) $(OBJRADIALRT) $(OBJNRECIP) $(OBJDUST)   $(OBJDEVLV) $(OBJALP)


#################################################
#                   RULES                       #
#################################################

all:	      diskevolfast

diskevolfast:     $(OBJ) makefile
	      $(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

diskevolfast.o:   $(SRCDISK) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  diskevolfast.F -o $@

radialrt2.o:   $(SRCRADIALRT) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  radialrt2.F -o $@

nrecip.o:     $(SRCNRECIP) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  nrecip.F -o $@


dusta.o:       $(SRCDUST) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  dusta.f90 -o $@

dustcollision.o:   $(SRCDEVLV) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  dustcollision.F -o $@

alpha.o:       $(SRCALP) makefile $(COMMONS)
	      $(FF) -c $(FFLAGS)  alpha.f   -o $@


