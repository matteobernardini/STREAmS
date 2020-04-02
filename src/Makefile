#----------------------------------------------------------------------
# SUPPORTED COMPILE: pgi-cuda, pgi-cuda-sync, pgi, intel, gnu, ibmxl, cray-cuda
# SUPPORTED MODES  : opt, debug [default=debug]
# SUPPORTED PRECISIONS (PREC) : single, double [default=double]
#----------------------------------------------------------------------
COMPILE = "gnu"
MODE    = "opt"
PREC    = "double"
#----------------------------------------------------------------------
ifeq (${PREC},"single") # mpif90 shipped with PGI compiler
    PRECFLAGS = -DSINGLE_PRECISION
endif
#----------------------------------------------------------------------

ifeq (${COMPILE},"cray-cuda") # mpif90 shipped with PGI compiler
    # Pitz-Daint instructions (2020)
    # module load pgi
    # module swap PrgEnv-cray PrgEnv-pgi
    CUDA_HOME = $(CRAY_CUDATOOLKIT_DIR)
    INCLUDE   = -I $(MPICH_DIR)/include 
    FCLIBS    = -L$(MPICH_DIR)/lib -lmpich 
    FC        = pgf90
    LDFLAGS   = -L$(CUDA_HOME)/lib64
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 -Mcuda=ccall -DUSE_CUDA -DCUDA_ASYNC -DNOCUDAAWAREMPI
    else
        OLEVEL = -Mcuda=ccall -DUSE_CUDA -DCUDA_ASYNC -g -DNOCUDAAWAREMPI
    endif
endif

ifeq (${COMPILE},"pgi-cuda") # mpif90 shipped with PGI compiler
    FC = mpif90 
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 -Mcuda=ccall,cuda10.0 -DUSE_CUDA -DCUDA_ASYNC
    else
        OLEVEL = -Mcuda=ccall,cuda10.0 -DUSE_CUDA -DCUDA_ASYNC
    endif
endif

ifeq (${COMPILE},"pgi-cuda-sync") # mpif90 shipped with PGI compiler
    FC = mpif90 
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 -Mcuda=ccall,cuda10.0 -DUSE_CUDA
    else
        OLEVEL = -Mcuda=ccall,cuda10.0 -DUSE_CUDA
    endif
endif

ifeq (${COMPILE},"pgi") # mpif90 shipped with PGI compiler or other MPI
    FC = mpif90 
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 
    else
        OLEVEL = -O2 -Mbounds -Ktrap=fp -g
    endif
endif

ifeq (${COMPILE},"intel") # mpiifort from IntelMPI
    FC = mpiifort
    ifeq (${MODE},"opt") 
        OLEVEL = -O3
    else
        OLEVEL = -O2 -check bounds,uninit -g -fpe0 -traceback
    endif
endif

ifeq (${COMPILE},"gnu") # mpif90 from MPI library
    FC = mpif90 
    ifeq (${MODE},"opt") 
        OLEVEL = -O3
    else
        OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core -ffpe-trap=zero,invalid,overflow
    endif
endif

ifeq (${COMPILE},"ibmxl") # mpif90 from MPI library
    FC = mpif90
    ifeq (${MODE},"opt") 
        OLEVEL = -O3
    else
        OLEVEL = -O2 -g
    endif
endif

RM = rm -f

#----------------------------------------------------------------------
# Base code
#----------------------------------------------------------------------
MAIN = main.o
OBJ_FILES = alloc.o bcdf.o bcextr.o bcfree.o bc.o bcrelax.o bcshk.o bcswap.o \
    bcwall.o bcwall_pl.o bcwall_staggered.o check_input.o computeav.o computemetrics.o \
    computeresidual.o constants.o df_par.o digital_filtering_util.o euler_i.o euler_j.o \
    euler_k.o finalize.o generategrid.o generateinflowrand.o generatewmean_channel.o \
    generatewmean.o init_channel.o init.o initurb.o init_windtunnel.o osw.o pgrad.o \
    prims.o printres.o readdf.o readgrid.o readinp.o readrst.o readstat.o rk.o sensor.o \
    setup.o solver.o startmpi.o stats.o step.o target_reystress.o updateghost.o utility.o \
    visflx.o writedf.o writefield.o writefieldvtk.o writegridplot3d.o writerst.o \
    writestatbl.o writestatchann.o writestat.o writestatzbl.o

OBJ_MODS = mod_streams.o mod_sys.o
OBJS = $(OBJ_MODS) $(OBJ_FILES) $(MAIN)
EXEC = streams
#----------------------------------------------------------------------
# Additional package
#----------------------------------------------------------------------
-include makefile.inc
#----------------------------------------------------------------------

$(EXEC): $(OBJS) Makefile
	$(FC) $(PRECFLAGS) $(OLEVEL) -o $@.exe $(OBJS) $(FCLIBS) $(LDFLAGS)

%.o:%.F90 Makefile
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(EXEC).exe $(OBJS) $(OBJ_MODS) $(GRID) *.mod
