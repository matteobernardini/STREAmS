```
!=============================================================
!
! ███████╗████████╗██████╗ ███████╗ █████╗ ███╗   ███╗███████╗
! ██╔════╝╚══██╔══╝██╔══██╗██╔════╝██╔══██╗████╗ ████║██╔════╝
! ███████╗   ██║   ██████╔╝█████╗  ███████║██╔████╔██║███████╗
! ╚════██║   ██║   ██╔══██╗██╔══╝  ██╔══██║██║╚██╔╝██║╚════██║
! ███████║   ██║   ██║  ██║███████╗██║  ██║██║ ╚═╝ ██║███████║
! ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝
!
! Supersonic TuRbulEnt Accelerated navier stokes Solver
!
! input file description
!
!=============================================================
```

`STREAmS` performs Direct Numerical Simulations of compressible turbulent flows in Cartesian geometry
solving the fully compressible Navier-Stokes equations. Currently, three canonical wall-bounded flows can be simulated:
* compressible turbulent channel flow
* compressible zero-pressure-gradient turbulent boundary layer
* supersonic oblique shock wave/turbulent boundary-layer interaction

STREAmS can be used in local environments but its main expected usage is with HPC architectures.

# Compiling

STREAmS requires (1) a Fortran compiler and (2) an MPI library. For the GPU CUDA version, the PGI compiler
is required (tested using PGI 19.4 or more recent compilers). To perform the compilation a Makefile with 
some predefined settings is available.
It is expected that the user only modifies the first three assignments to select the compilation type.
For example:

```
COMPILE = "gnu"
MODE    = "opt"
PREC    = "double"
```

`COMPILE` currently supports these choices:

* `pgi-cuda`: PGI compiler Cuda Fortran asyncrhonous version with MPI library (tested with OpenMPI provided by PGI)
* `pgi-cuda-sync`: PGI compiler Cuda Fortran synchornous version with MPI library (tested with OpenMPI provided by PGI)
* `pgi`: PGI compiler CPU version with OpenMPI library (tested with OpenMPI provided by PGI)
* `intel`: Intel compiler with MPI library (tested with IntelMPI library)
* `gnu`: gnu compiler with MPI library (tested with OpenMPI)
* `ibmxl`: XL IBM compiler with MPI library (tested with IBM MPI library)
* `cray-cuda`: PGI compiler with Cray mpich library without support of CUDA-Aware MPI (currently oriented to Pitz-Daint cluster)

`MODE` can be one of:

* `opt`: optimized compilation, for standard production runs
* `debug`: debugging compilation, with run-time checks enabled and backtracing when available 

`PREC` defines the precision of the float numbers:

* `double`: double precision floating point numbers
* `single`: single precision floating point numbers

Additional Makefile customizations can be made manually modifying the Makefile.

# Running

Once compiled, STREAmS can be usually executed using a standard MPI launcher, e.g. `mpirun`.
In order to execute, in addition to the executable, in the running folder you need:
* `input.dat`: file defining the physical and numerical setup of the simulation, to be customized
according to the desired needs. For an explaination of the input.dat parameters, read below.
Some examples of input.dat files are available in the `examples` folder.
* `database_bl.dat`: only required for the `boundary layer` and `shock-boundary layer interaction` flows
(this file does not have to be modified by the user). The file is available in the examples folder.

To run a simulation, type, e.g.:
```
mpirun -np 8 ./streams
```
or (for SLURM jobs)
```
srun ./streams
```

For CUDA versions in cluster environments, you must distribute MPI processes according to the number of GPUs 
available for each node. For CINECA Marconi-100 cluster -- 4 GPUS per node --  a possible submission script using
8 GPUs is:

```
#!/bin/bash
#SBATCH -N 2 
#SBATCH --tasks-per-node 4
#SBATCH --mem=64G 
#SBATCH --partition=debug 
#SBATCH --time=00:30:00 
#SBATCH --gres=gpu:4
#SBATCH --partition=debug 
module load profile/global pgi
srun ./streams
```

For CSCS Pitz-Daint cluster -- 1 GPU per node -- a possible submission script using 8 GPUs is:
```
#!/bin/bash
#SBATCH -N 8
#SBATCH --tasks-per-node 1 
#SBATCH --mem=15G 
#SBATCH --partition=debug 
#SBATCH --time=00:30:00 
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu 
module swap PrgEnv-cray PrgEnv-pgi
srun ./streams
```

# Preparing input.dat file

`iflow` defines the type of flow. 0 = channel, 1 = boundary layer, 2 = shock/boundary-layer interaction

`rlx rly rlz` real numbers defining the size of the Cartesian domain along the three directions 
(x=streamwise, y=wall-normal, z=spanwise)

`nxmax  nymax  nzmax` integer numbers defining the number of grid nodes along each direction
 
 `nymaxwr   rlywr    dyp_target` specify the wall-normal mesh features. dyp_target is the desired spacing at the wall in inner units.
 nymaxwr denotes the number of grid points in the wall-resolved region, ranging from y=0 (wall) up to y=rlywr, where a sinh mapping
 is applied. Both nymaxwr and rlywr can be specified when iflow = 1 (turbulent boundary layer), then a geometric progession is applied from
 y=rlywr up to y = rly. When iflow = 2 (shock/boundary-layer interaction), nymaxwr must be specified but rlywr is automatically computed.
 In this case a constant spacing in the wall-normal direction is applied from y=rlywr up to y=rly.

 `nblocks(1)  nblocks(3)` define the MPI decomposition along x (streamwise) and z (spanwise). These numbers 
 must be consistent to nxmax, nymax, and nzmax. In particular the following divisions must have zero remainder:
 nxmax/nblocks(1), nzmax/nzmax. Furthermore, for iflow>0 cases also nxmax/nblocks(3) and nymax/nblocks(1)
 cannot have remainders.  The product  nblocks(1)\*nblocks(3) must equal the total number of launched MPI processes.

 `tresduc  ximp  deflec` tresduc is the weno threshold (weno is active if the shock sensor value exceeds tresduc),
 ximp is the shock abscissa and deflec is the shock angle. ximp and deflec are used only if iflow > 0.

 `idiski ncyc   cfl  nstep nprint` idiski selects the start type (0=init, 1=restart), ncyc the total number
 of iterations, cfl is the stability limit, nstep the interval of iterations between dt evaluation,  and nprint 
 the interval of iteration between output residual print

 `rm      re (friction)  twall/taw` Mach number, $Re_\tau$ number, and the ratio between wall temperature and
 adiabatic wall temperature.

 `istat  nstat` interval of iterations between statistics evaluation, number of abscissas for 
 statistics extractions. nstat meaningful only if iflow > 0 

 `xstat` abscissas of statistic. Meaningful only if iflow > 0
  20. 30. 40. 50. 60. 70. 80. 87.

 `dtsave dtsave_restart plot3d vtk` dt_save is the time interval between field output, dtsave_restart is the time interval
 between output for simulation restart, plot3d is a integer which, if larger than 0, activates the plot3d format output,
 vtk is the corresponding integer to activate vtk format output

  `irand` if < 0 produce not reproducible random sequences based on current time, if >=0 produces random sequences which
are reproducible across different runs with the same configuration

Look at the input files in the examples folder to start with typical cases

# Understanding the output files

 Plot3d binary output files -- e.g. plot3dgrid.xyz and field_0001.q -- can be read using Tecplot or Paraview

 VTK Rectilinear grid files -- e.g. field_0001.vtr -- can be read using Paraview.

 Other statistics files are automatically produced according to the input configuration.
 
 ## Channel flow output files
### Residual file `output_streams.dat`
 When running channel flow cases (`iflow`=0) the residuals file (output_streams.dat) contains seven columns:
 1) number of cycles,
 2) elapsed time,
 3) streamwise momentum residual,
 4) pressure gradient,
 5) bulk density (conserved to machine accuracy),
 6) bulk flow velocity (conserved to machine accuracy),
 7) bulk temperature.
 For instance, plotting the fourth column vs. the second allows the user to check the time history of the pressure gradient.



