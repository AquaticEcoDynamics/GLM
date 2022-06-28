###############################################################################
#                                                                             #
# Makefile for glm                                                            #
#                                                                             #
#  Part of GLM (General Lake Model)                                           #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2022 -  The University of Western Australia               #
#                                                                             #
#   GLM is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GLM is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

ifeq ($(shell uname),Linux)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),Darwin)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),FreeBSD)
  OSTYPE=$(shell uname -s)
else
  OSTYPE=$(shell uname -o)
endif

BUILDDATE=$(shell date -u +%Y%m%d-%H%MUTC)

ifeq ($(WITH_PLOTS),)
  WITH_PLOTS=true
  ifeq ($(WITH_XPLOTS),)
    WITH_XPLOTS=true
  endif
  ifeq ($(PLOTDIR),)
    PLOTDIR=../libplot
  endif
endif

ifeq ($(UTILDIR),)
  UTILDIR=../libutil
endif

ifeq ($(WITH_CHECKS),)
  ifeq ($(DEBUG),true)
    WITH_CHECKS=true
  else
    WITH_CHECKS=false
  endif
endif

srcdir=src
incdir=src
objdir=obj
moddir=mod

TARGETS=glm
DEFINES=
FINCLUDES=-I$(UTILDIR)/include
CINCLUDES=-I$(UTILDIR)/include
LIBS=-L$(UTILDIR)/lib -lutil
GLM_DEPS=$(UTILDIR)/lib/libutil.a
ifeq ($(WITH_PLOTS),true)
  DEFINES+=-DPLOTS
  ifeq ($(WITH_XPLOTS),true)
    DEFINES+=-DXPLOTS
  endif
  GLM_DEPS+=$(PLOTDIR)/lib/libplot.a
  CINCLUDES+=-I$(PLOTDIR)/include
endif

ifeq ($(OSTYPE),Darwin)
  ifeq ($(HOMEBREW),true)
     FINCLUDES+=-I/usr/local/include
     CINCLUDES+=-I/usr/local/include
     LIBS+=-L/usr/local/lib
  else
     FINCLUDES+=-I/opt/local/include
     CINCLUDES+=-I/opt/local/include
     LIBS+=-L/opt/local/lib
  endif
  #EXTRALINKFLAGS=-Wl,-no_compact_unwind
  EXTRALINKFLAGS=-Wl,-no_compact_unwind,-headerpad_max_install_names
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else ifeq ($(OSTYPE),Msys)
  CINCLUDES+=-I../win-3rd-party/x64-Release/include
  LIBS+=-L../win-3rd-party/x64-Release/lib
else ifeq ($(OSTYPE),FreeBSD)
  FINCLUDES+=-I../flang_extra/mod
  FINCLUDES+=-I/usr/local/include
  CINCLUDES+=-I/usr/local/include
  LIBS+=-L../flang_extra/ -lflang_extra -L/usr/local/lib
else
  CINCLUDES+=-I/usr/local/include
  EXTRALINKFLAGS=-Wl,--export-dynamic
  SHARED=-shared
  so_ext=so
endif

ifeq ($(FABM),true)
  ifeq ($(FABMDIR),)
    FABMDIR=../fabm-git
  endif
  DEFINES+=-DFABM

  FABMLIB=fabm
  ifeq ($(DEBUG),true)
    WITH_CHECKS=true
  endif

  FINCLUDES+=-I$(FABMDIR)/include -I$(FABMDIR)/src/drivers/glm -I$(FABMDIR)/build/modules
  FABMLIBS=-L$(FABMDIR)/build -l$(FABMLIB)

  ifeq ($(USE_DL),true)
    FABMTARGETS=libglm_wq_fabm.${so_ext}
  endif
endif

EXTFFLAGS=
ifeq ($(AED),true)
  DEFINES+=-DAED

  AEDWATDIR=../libaed-water
  FINCLUDES+=-I$(AEDWATDIR)/include -I$(AEDWATDIR)/mod
  AEDLIBS=-L$(AEDWATDIR)/lib -laed-water
  ifdef AEDBENDIR
    AEDLIBS+=-L$(AEDBENDIR)/lib -laed-benthic
  else
    EXTFFLAGS+=-DNO_BENTHIC
  endif
  ifdef AEDRIPDIR
    AEDLIBS+=-L$(AEDRIPDIR)/lib -laed-riparian
  else
    EXTFFLAGS+=-DNO_RIPARIAN
  endif
  ifdef AEDDMODIR
    AEDLIBS+=-L$(AEDDMODIR)/lib -laed-demo
  else
    EXTFFLAGS+=-DNO_DEMO
  endif
  ifdef AEDDEVDIR
    AEDLIBS+=-L$(AEDDEVDIR)/lib -laed-dev
  else
    EXTFFLAGS+=-DNO_DEV
  endif
  #AEDLIBS+=-ldl

  ifeq ($(USE_DL),true)
    AEDTARGETS=libglm_wq_aed.${so_ext}
  endif

  GLM_DEPS+=$(AEDWATDIR)/lib/libaed-water.a
endif

ifeq ($(AED2),true)
  DEFINES+=-DAED2

  ifeq ($(AED2DIR),)
    AED2DIR=../libaed2
  endif

  FINCLUDES+=-I$(AED2DIR)/include -I$(AED2DIR)/mod
  AED2LIBS=-L$(AED2DIR)/lib -laed2
  ifneq ("$(wildcard ${AED2PLS}/Makefile)","")
    AED2PLBS=-L${AED2PLS}/lib -laed2+
  endif

  ifeq ($(USE_DL),true)
    AED2TARGETS=libglm_wq_aed2.${so_ext}
  endif

  GLM_DEPS+=$(AED2DIR)/lib/libaed2.a
endif

FLIBS=
# Select specific compiler bits
ifeq ($(F90),ifort)
  LINK=$(CC)
  FINCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-warn all -module ${moddir} -static-intel -mp1 -stand f08 -warn nounused $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check bounds -check noarg_temp_created
  endif
  FFLAGS+=-real-size 64
  FLIBS+=-L/opt/intel/lib
  FLIBS+=-lifcore -lsvml -lifport
  FLIBS+=-limf -lintlc -liomp5  -lifport
  OMPFLAG=-qopenmp
  #EXTFFLAGS=-warn-no-unused-dummy-argument
else ifeq ($(F90),flang)
# LINK=$(FC) -fno-fortran-main
  LINK=$(CC)
  DEBUG_FFLAGS=-g -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-module ${moddir} $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
  FLIBS+=-L/usr/local/flang/lib
  FLIBS+=-lexecinfo -lpgmath -lomp
  FLIBS+=-lflang -lflangrti -lflangADT -lflangArgParser
else ifeq ($(F90),pgfortran)
  LINK=$(CC)
  DEBUG_FFLAGS=-g -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-module ${moddir} $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
  FLIBS+=-L/opt/pgi/linux86-64/18.10/lib
  FLIBS+=-lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902
  FLIBS+=-lpgftnrtl -lpgmp -lnuma -lpgmath -lpgc
else
  LINK=$(FC)
  DEBUG_FFLAGS=-g -fbacktrace -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-Wall -J ${moddir} -Wno-c-binding-type -ffree-line-length-none -std=f2008 $(DEFINES) $(FINCLUDES) -fall-intrinsics
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-fcheck=all,no-array-temps
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
  OMPFLAG=-fopenmp
  FLIBS+=-lgfortran -lgomp
  EXTFFLAGS+=-Wno-unused-dummy-argument
endif

ifneq ($(USE_DL),true)
  WQLIBS=$(AEDLIBS) $(FABMLIBS) $(AED2LIBS)
endif

ifeq ($(DEBUG),true)
  ifeq ($(OSTYPE),FreeBSD)
    DEBUG_CFLAGS=-g -fsanitize=address -DDEBUG=1
    LIBS+=-fsanitize=address -static-libsan
  else
    DEBUG_CFLAGS=-g -fbounds-check -DDEBUG=1
  endif
  OPT_CFLAGS=
  OPT_FFLAGS=
else
  DEBUG_FFLAGS=
  DEBUG_CFLAGS=
  # OPT_CFLAGS=-O4 -Ofast -frounding-math
  OPT_CFLAGS=-O3
  # OPT_CFLAGS=
  # OPT_FFLAGS=
endif

LIBS+=-lnetcdf
# If variable NETCDFLIB is not empty, use it to
# set path to the library
ifneq ($(NETCDFLIB),)
  LIBS+=-L$(NETCDFLIB)
endif


ifeq ($(PLOTDIR),)
  PLOTDIR=../../libplot
endif

RES=
RESP=
ifeq ($(WITH_PLOTS),true)
  #LIBS+=-L$(PLOTDIR)/lib -lplot -lgd -lpng -ljpeg -lm
  LIBS+=-L$(PLOTDIR)/lib -lplot -lgd -ljpeg -lm
  ifeq ($(WITH_XPLOTS),true)
    ifeq ($(OSTYPE),Darwin)
      LIBS+=-framework Cocoa
    else ifeq ($(OSTYPE),Msys)
      LIBS+=-lpng16 -lz -lcomdlg32 -lgdi32
      RES=${objdir}/glm_rc.o
      RESP=${objdir}/glm+_rc.o
    else
      LIBS+=-lX11
    endif
  endif
endif
ifeq ($(FENCE),true)
  LIBS+=-lefence
endif

CFLAGS=-Wall -I$(UTILDIR) -I$(PLOTDIR) $(CINCLUDES) $(DEFINES) $(DEBUG_CFLAGS) $(OPT_CFLAGS)
FFLAGS+=$(DEBUG_FFLAGS) $(OPT_FFLAGS)

OBJS=${objdir}/glm_globals.o \
     ${objdir}/glm_util.o \
     ${objdir}/glm_csv.o \
     ${objdir}/glm_mobl.o \
     ${objdir}/glm_mixu.o \
     ${objdir}/glm_wqual.o \
     ${objdir}/glm_layers.o \
     ${objdir}/glm_surface.o \
     ${objdir}/glm_input.o \
     ${objdir}/glm_plot.o \
     ${objdir}/glm_output.o \
     ${objdir}/glm_ncdf.o \
     ${objdir}/glm_lnum.o \
     ${objdir}/glm_init.o \
     ${objdir}/glm_flow.o \
     ${objdir}/glm_mixer.o \
     ${objdir}/glm_deep.o \
     ${objdir}/glm_stress.o \
     ${objdir}/glm_bird.o \
     ${objdir}/glm_model.o \
     ${objdir}/glm_types.o \
     ${objdir}/glm_const.o \
     ${objdir}/glm_debug.o \
     ${objdir}/glm_balance.o \
     ${objdir}/glm_main.o

ifeq ($(USE_DL),true)
  LIBS+=-ldl
  CFLAGS+=-DUSE_DL_LOADER=1
  FFLAGS+=-DUSE_DL_LOADER=1
  TARGETS+=$(AEDTARGETS) $(FABMTARGETS)
else
  ifeq ($(AED),true)
    OBJS+=${objdir}/glm_zones.o
  else ifeq ($(FABM),true)
    OBJS+=${objdir}/glm_zones.o
  endif
  ifeq ($(AED2),true)
    OBJS+=${objdir}/glm_aed2.o
  endif
  ifeq ($(AED),true)
    OBJS+=${objdir}/glm_aed.o \
          ${objdir}/aed_external.o
  endif
  ifeq ($(FABM),true)
  # OBJS+=${objdir}/glm_fabm.o ${objdir}/ode_solvers.o
    OBJS+=${objdir}/ode_solvers.o ${objdir}/glm_fabm.o
  endif
endif

all: $(TARGETS)

lib:
	@mkdir lib

${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

glm: ${objdir} ${moddir} $(OBJS) $(GLM_DEPS) $(RES)
	$(LINK) -o $@ $(EXTRALINKFLAGS) $(OBJS) $(RES) $(LIBS) $(WQLIBS) $(FLIBS)

glm+: ${objdir} ${moddir} $(OBJS) $(GLM_DEPS) $(RESP)
	$(LINK) -o $@ $(EXTRALINKFLAGS) $(OBJS) $(RESP) $(LIBS) $(WQLIBS) $(FLIBS)

clean: ${objdir} ${moddir}
	@touch ${objdir}/1.o ${moddir}/1.mod 1.t 1__genmod.f90 glm 1.${so_ext} glm_test_bird macos/glm.app macos/glm+.app
	@touch debian/.debhelper debian/files
	@touch debian/glm debian/glm.debhelper.log debian/glm.substvars
	@touch debian/glm+ debian/glm+.debhelper.log debian/glm+.substvars
	@/bin/rm ${moddir}/*.mod ${objdir}/*.o *.t *__genmod.f90 *.${so_ext} glm_test_bird
	@/bin/rm -rf debian/.debhelper debian/files
	@/bin/rm -rf debian/glm debian/glm.debhelper.log debian/glm.substvars
	@/bin/rm -rf debian/glm+ debian/glm+.debhelper.log debian/glm+.substvars
	@echo Made clean

distclean: clean
	@/bin/rm -rf ${objdir} ${moddir} glm glm+ macos/glm.app macos/glm+.app

${objdir}/%.o: ${srcdir}/%.F90 ${incdir}/glm.h
	$(FC) -fPIC $(FFLAGS) $(EXTRA_FFLAGS) -D_FORTRAN_SOURCE_ -c $< -o $@

${objdir}/glm_main.o: ${srcdir}/glm_main.c ${incdir}/glm.h
	$(CC) -fPIC -DBUILDDATE=\"${BUILDDATE}\" $(CFLAGS) $(EXTRA_FLAGS) -c $< -o $@

${objdir}/glm_rc.o: win/glm.rc
	windres $< -o $@

${objdir}/glm+_rc.o: win/glm+.rc
	windres $< -o $@

${objdir}/%.o: ${srcdir}/%.c ${incdir}/glm.h
	$(CC) -fPIC $(CFLAGS) $(EXTRA_FLAGS) -c $< -o $@

%.${so_ext}:
	$(LD) ${SHARED} $(LDFLAGS) \
                        -o $@ $^ -L/opt/intel/lib/intel64/ $(LIBS)

#                        -E -Bdynamic -undefined suppress -o $@ $^ -L/opt/intel/lib/intel64/ $(LIBS)

# Build rules

libglm_wq_aed.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_aed.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(AEDLIBS)

libglm_wq_aed+.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_aed.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(AEDPLBS)

libglm_wq_fabm.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_fabm.o ${objdir}/ode_solvers.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(FABMLIBS)

# special needs dependancies

${objdir}/aed_external.o: ../libaed-water/src/aed_external.F90
	$(FC) -fPIC $(FFLAGS) $(EXTFFLAGS) -D_FORTRAN_SOURCE_ $(OMPFLAG) -c $< -o $@

${objdir}/glm_globals.o: ${srcdir}/glm_globals.c ${incdir}/glm_globals.h ${incdir}/glm.h
${objdir}/glm_plugin.o: ${srcdir}/glm_plugin.c ${incdir}/glm_plugin.h ${incdir}/glm.h
${objdir}/glm_mixer.o: ${srcdir}/glm_mixer.c ${incdir}/glm_mixer.h ${incdir}/glm.h ${srcdir}/glm_debug.h
