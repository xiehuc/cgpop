#
# File:  sgi.gnu
#
# The commenting in this file is intended for occasional maintainers who 
# have better uses for their time than learning "make", "awk", etc.  There 
# will someday be a file which is a cookbook in Q&A style: "How do I do X?" 
# is followed by something like "Go to file Y and add Z to line NNN."
#
FC = f90
LD = f90
CC = cc
Cp = /usr/bin/cp
Cpp = /lib/cpp -P -traditional
#AWK = /fs/tools/sparc.sunos-5.5.1/gawk-3.0.2/bin/awk
AWK = /usr/xpg4/bin/awk

ABI =
COMMDIR = serial
 
# These have been loaded as a module so no values necessary
#NETCDFINC = -I/usr/local/include -I/usr/local/lib64/r4i4
#NETCDFLIB = -L/usr/local/lib64/r4i4 
NETCDFINC = -I/fs/tools/sparc.solaris-7/netcdf-3.5.1-beta/lib
NETCDFLIB = -L/fs/tools/sparc.solaris-7/netcdf-3.5.1-beta/lib 
NETCDFMOD = -M/fs/tools/sparc.solaris-7/netcdf-3.5.1-beta/lib

#  Enable MPI library for parallel code, yes/no.

MPI = no

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

# Enable hardware trace library 
PAPIDIR = /fs/projects/css/dennis/PAPI
TRACE = -D_HTRACE
TLIBS = -L/fs/projects/css/dennis/validate/lib -lhtrace $(PAPIDIR)/lib/libpapi.a -lcpc 

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled

Cpp_opts =  $(DCOUPL) $(TRACE)

Cpp_opts := $(Cpp_opts) -DPOSIX $(NETCDFINC)
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI)

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -xO5 -fast -xarch=v8plusa -xchip=native -dalign
else
  CFLAGS := $(CFLAGS) -g
endif
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------
 
FBASE = $(ABI) $(NETCDFMOD) -moddir=$(ObjDepDir)

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) -DEBUG:trap_uninitialized=ON
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -fast -e -xarch=v8plusa -xchip=native -dalign -stackvar
else
# FFLAGS := $(FBASE) -g 
  FFLAGS := $(FBASE) -g -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) $(FBASE) -fast 
 
#LIBS = $(NETCDFLIB) -lnetcdf -lX11
LIBS = $(NETCDFLIB) -lnetcdf $(TLIBS)
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) -lfpe
endif
 
#LDLIBS = $(TARGETLIB) $(LIBRARIES) $(LIBS)
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------
#
#                           Explicit Rules for Compilation Problems
#
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#
#                           Implicit Rules for Compilation
#
#----------------------------------------------------------------------------
 
# Cancel the implicit gmake rules for compiling
%.o : %.f
%.o : %.f90
%.o : %.c

%.o: %.f
	@echo Solaris SERIAL Compiling with implicit rule $<
	@$(FC) $(FFLAGS) -c $<
	@if test -f *.mod; then mv -f *.mod $(ObjDepDir); fi
 
%.o: %.f90
	@echo Solaris SERIAL Compiling with implicit rule $<
	$(FC) $(FFLAGS) -c $<
	@if test -f *.mod; then mv -f *.mod $(ObjDepDir); fi
 
%.o: %.c
	@echo Solaris SERIAL Compiling with implicit rule $<
	$(CC) $(Cpp_opts) $(CFLAGS) -c $<

#----------------------------------------------------------------------------
#
#                           Implicit Rules for Dependencies
#
#----------------------------------------------------------------------------
 
ifeq ($(OPTIMIZE),yes)
  DEPSUF = .do
else
  DEPSUF = .d
endif

# Cancel the implicit gmake rules for preprocessing

%.c : %.C
%.o : %.C

%.f90 : %.F90
%.o : %.F90

%.f : %.F
%.o : %.F

%.h : %.H
%.o : %.H

# Preprocessing  dependencies are generated for Fortran (.F, F90) and C files
$(SrcDepDir)/%$(DEPSUF): %.F
	@echo 'Solaris SERIAL Making depends for preprocessing' $<
	$(Cpp) $(Cpp_opts) $< >$(TOP)/compile/$*.f
	@echo '$(*).f: $(basename $<)$(suffix $<)' > $(SrcDepDir)/$(@F)

$(SrcDepDir)/%$(DEPSUF): %.F90
	@echo 'Solaris SERIAL Making depends for preprocessing' $<
	$(Cpp) $(Cpp_opts) $< >$(TOP)/compile/$*.f90
	echo '$(*).f90: $(basename $<)$(suffix $<)' > $(SrcDepDir)/$(@F)

$(SrcDepDir)/%$(DEPSUF): %.C
	@echo 'Solaris SERIAL Making depends for preprocessing' $<
#  For some reason, our current Cpp options are incorrect for C files.
#  Therefore, let the C compiler take care of #ifdef's, etc.  Just copy.
#	@$(Cpp) $(Cpp_opts) $< >$(TOP)/compile/$*.c
	@$(Cp) $< $(TOP)/compile/$*.c
	@echo '$(*).c: $(basename $<)$(suffix $<)' > $(SrcDepDir)/$(@F)

# Compiling dependencies are generated for all normal .f files
$(ObjDepDir)/%$(DEPSUF): %.f
	@if test -f $(TOP)/compile/$*.f;  \
       then : ; \
       else $(Cp) $< $(TOP)/compile/$*.f; fi
	@echo 'Solaris SERIAL Making depends for compiling' $<
	@$(AWK) -f $(TOP)/fdepends.awk -v NAME=$(basename $<) -v ObjDepDir=$(ObjDepDir) -v SUF=$(suffix $<) -v DEPSUF=$(DEPSUF) $< > $(ObjDepDir)/$(@F)

# Compiling dependencies are generated for all normal .f90 files
$(ObjDepDir)/%$(DEPSUF): %.f90
	@if test -f $(TOP)/compile/$*.f90;  \
       then : ; \
       else $(Cp) $< $(TOP)/compile/$*.f90; fi
	@echo 'Solaris SERIAL Making depends for compiling' $<
	@$(AWK) -f $(TOP)/fdepends.awk -v NAME=$(basename $<) -v ObjDepDir=$(ObjDepDir) -v SUF=$(suffix $<) -v DEPSUF=$(DEPSUF) $< > $(ObjDepDir)/$(@F)

# Compiling dependencies are also generated for all .c files, but 
# locally included .h files are not treated.  None exist at this 
# time.  The two .c files include only system .h files with names 
# delimited by angle brackets, "<...>"; these are not, and should 
# not, be analyzed.  If the c programming associated with this code 
# gets complicated enough to warrant it, the file "cdepends.awk" 
# will need to test for includes delimited by quotes.
$(ObjDepDir)/%$(DEPSUF): %.c
	@if test -f $(TOP)/compile/$*.c;  \
       then : ; \
       else $(Cp) $< $(TOP)/compile/$*.c; fi
	@echo 'Solaris SERIAL Making depends for compiling' $<
	@echo '$(*).o $(ObjDepDir)/$(*)$(DEPSUF): $(basename $<)$(suffix $<)' > $(ObjDepDir)/$(@F)
