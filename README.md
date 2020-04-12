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

 `idiski ncyc   cfl  nstep nprint` idiski selects the start type (0=init, 1=restart, 2=collect statistics at runtime), ncyc the total number
 of iterations, cfl is the stability limit, nstep the interval of iterations between dt evaluation,  and nprint 
 the interval of iteration between output residual print

 `rm      re (friction)  twall/taw` Mach number, <img src="svgs/1f71bd9db75247d7e8fed8df71d9e9f8.svg?invert_in_darkmode&sanitize=true" align=middle width=27.55203pt height=22.38192pt/> number, and the ratio between wall temperature and
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
### Mean flow statistics `channstat.prof`
    The file `channstat.prof` contains the mean channel flow statistics.
    This file is printed at the end of each run and it contains the mean flow statistics averaged
    in the homogeneous spatial directions and in time (statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1).
    The file `channstat.prof` contains 15 columns:
1) <img src="svgs/08b335933c05df68e8439a765e868cdf.svg?invert_in_darkmode&sanitize=true" align=middle width=26.243415pt height=24.56553pt/> is the wall-normal coordinate, normalized with the channel half width
2) <img src="svgs/9cbaafb1656e7b37560d2705ed2653e1.svg?invert_in_darkmode&sanitize=true" align=middle width=72.450675pt height=26.12412pt/> the wall-normal coordinate in viscous units
3) <img src="svgs/c29b0a6ba7b896aed4d1ff94a957d475.svg?invert_in_darkmode&sanitize=true" align=middle width=82.21653pt height=28.25757pt/> the wall-normal coordinate transformed according to Trettel & Larsson in viscous units
4) <img src="svgs/e1210f0fa24ba17ea9fd4da41c310dab.svg?invert_in_darkmode&sanitize=true" align=middle width=32.733195pt height=24.56553pt/> the mean streamwise velocity averaged according to Favre, normalized with the bulk flow velocity
5) <img src="svgs/edd39792d631fff9fc9340cbcd03af76.svg?invert_in_darkmode&sanitize=true" align=middle width=34.33848pt height=24.56553pt/> the mean streamwise velocity averaged according to Favre, normalized with the friction velocity
6) <img src="svgs/5168dd0c22b7ee67557ebdc76a52a925.svg?invert_in_darkmode&sanitize=true" align=middle width=48.956655pt height=24.56553pt/> the mean streamwise velocity transformed according to van Driest, normalized with the friction velocity
7) <img src="svgs/2fb4a896c47e22ceb0b1ad62e552a0bf.svg?invert_in_darkmode&sanitize=true" align=middle width=44.6589pt height=24.56553pt/> the mean streamwise velocity transformed according to Trettel & Larsson in viscous units
8) <img src="svgs/6d150b0d92a9b0f9f8ce9bf06e99af7c.svg?invert_in_darkmode&sanitize=true" align=middle width=34.968945pt height=24.56553pt/> the mean density profile, normalized with the mean wall density
9) <img src="svgs/9c48fce70048f82f658109eef23b7b71.svg?invert_in_darkmode&sanitize=true" align=middle width=69.25347pt height=34.08471pt/> the Favre streamwise Reynolds stress, normalized with the wall-shear stress
10) <img src="svgs/7c2b51045dc88620a88a255165e20708.svg?invert_in_darkmode&sanitize=true" align=middle width=67.54869pt height=34.08471pt/> the Favre wall-normal Reynolds stress, normalized with the wall-shear stress
11) <img src="svgs/f243b45c92e4fdc93b7e8280bcebb799.svg?invert_in_darkmode&sanitize=true" align=middle width=74.85456pt height=36.5409pt/> the Favre spanwise Reynolds stress, normalized with the wall-shear stress
12) <img src="svgs/b1254d1845c343e24eb7465dad8339d0.svg?invert_in_darkmode&sanitize=true" align=middle width=68.40108pt height=34.08471pt/> the Favre shear Reynolds stress, normalized with the wall-shear stress
13) <img src="svgs/41ea061742760dfbfe5de2c01f1c8bc1.svg?invert_in_darkmode&sanitize=true" align=middle width=39.43071pt height=27.72594pt/>  the mean temperature profile, normalized with the wall temperature
14) <img src="svgs/6673acc30acf7858774f24361b9bcecf.svg?invert_in_darkmode&sanitize=true" align=middle width=46.78971pt height=30.4656pt/> the density fluctuations normalized with the mean wall density
15) <img src="svgs/6b014bcc811232a0bb81b6b3ab828c69.svg?invert_in_darkmode&sanitize=true" align=middle width=51.25164pt height=30.4656pt/> The temperature fluctuations, normalized with the wall temperature
## Boundary layer output files
### Residual file `output_streams.dat`
 When running boundary layer cases (`iflow`=1) the residuals file (output_streams.dat) contains three columns:
 1) number of cycles, 
 2) elapsed time, 
 3) streamwise momentum residual, 
### Wavplot files 
  The files wavplot_xxx_yyy.dat are ASCII files in Tecplot data format.
  The numbers xxx and yyy refer to the Cartesian MPI block to which the file belong.
  The files are printed at the end of each run. Statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1.
  The files contain the spanwise Favre averaged velocity components as a function
  of the streamwise and wall-normal spatial coordinates. The file contains four variables:
1) <img src="svgs/a7de9ea243c799172d54f6c73ca17aa2.svg?invert_in_darkmode&sanitize=true" align=middle width=31.355115pt height=24.56553pt/> streamwise coordinate normalized with the inlet boundary layer thickness
2) <img src="svgs/e3a9a01189dcd6e5f33eeafb6d00056a.svg?invert_in_darkmode&sanitize=true" align=middle width=30.61443pt height=24.56553pt/> wall-normal coordinate normalized with the inlet boundary layer thickness
3) <img src="svgs/09a4eb9ca66b0dde022422feeb68d606.svg?invert_in_darkmode&sanitize=true" align=middle width=33.50193pt height=24.56553pt/> the Favre averaged streamwise velocity component, normalized with the free-stream velocity
4) <img src="svgs/b8d101a623e31826e689ce804017c38b.svg?invert_in_darkmode&sanitize=true" align=middle width=32.64954pt height=24.56553pt/> the Favre averaged wall-normal velocity  component, normalized with the free-stream velocity
5) <img src="svgs/368d663ad6f8d1f9c03a3b43ddc750bc.svg?invert_in_darkmode&sanitize=true" align=middle width=36.30264pt height=24.56553pt/> the Favre averaged spanwise velocity  component, normalized with the free-stream velocity
### cf files
  The files cf_xxx.dat are ASCII files and the number xxx refers to the Cartesian streamwise MPI block to which the file belong.
  These files is printed at the end of each run. Statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1.
  The file contains 13 columns:
1) <img src="svgs/a7de9ea243c799172d54f6c73ca17aa2.svg?invert_in_darkmode&sanitize=true" align=middle width=31.355115pt height=24.56553pt/> streamwise coordinate normalized with the inlet boundary layer thickness
2) <img src="svgs/95b92a3b01f370498645efb59e52ff36.svg?invert_in_darkmode&sanitize=true" align=middle width=19.376115pt height=22.38192pt/> skin-friction coefficient
3) <img src="svgs/b91858acbc5918f19aa579c1c71a454a.svg?invert_in_darkmode&sanitize=true" align=middle width=19.19247pt height=22.74591pt/> friction Reynolds number 
4) <img src="svgs/7b9a0316a2fcd7f01cfd556eedf72e96.svg?invert_in_darkmode&sanitize=true" align=middle width=14.94405pt height=22.38192pt/> Compressible shape factor
5) <img src="svgs/b76f125299e13d61d8c9b80929221c5f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.195295pt height=22.38192pt/> Incompressible shape factor 
6) <img src="svgs/018562a8c71022027090ab072944bb8e.svg?invert_in_darkmode&sanitize=true" align=middle width=43.20096pt height=24.56553pt/> boundary layer thickness, normalized with the inlet boundary layer thickness
7) <img src="svgs/27e556cf3caa0673ac49a8f0de3c73ca.svg?invert_in_darkmode&sanitize=true" align=middle width=8.1430305pt height=22.74591pt/> Compressible displacement thickness
8) <img src="svgs/19889dd126af8c9b2d5d30c60cb5606d.svg?invert_in_darkmode&sanitize=true" align=middle width=26.26998pt height=22.74591pt/> Compressible momentum thickness
9) <img src="svgs/89144b2a5e5b3f3a3663daf61886c11b.svg?invert_in_darkmode&sanitize=true" align=middle width=16.739745pt height=14.10255pt/> friction velocity
10) <img src="svgs/a90189de2f00cdfb33219eb692f700a7.svg?invert_in_darkmode&sanitize=true" align=middle width=45.36114pt height=22.38192pt/> Reynolds number based on the incompressible momentum thickness 
11) <img src="svgs/e2c20c0828264306966bbd3444f28680.svg?invert_in_darkmode&sanitize=true" align=middle width=37.96287pt height=22.38192pt/> Incompressible friction coefficient, according to van Driest II transofrmation
12) <img src="svgs/a90189de2f00cdfb33219eb692f700a7.svg?invert_in_darkmode&sanitize=true" align=middle width=45.36114pt height=22.38192pt/> Reynolds number based on the compressible momentum thickness
13) <img src="svgs/ac9c81a75856142b321298cf3bd14b97.svg?invert_in_darkmode&sanitize=true" align=middle width=62.57625pt height=42.84918pt/> Mean pressure rms at the wall normalized with the wall-shear stress 
### stat files
   The files stat_nnn.dat are ASCII containing the mean boundary layer statistics.
   The files are printed at the end of each run. Statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1.
   The number nnn indicates the global mesh index in the streamwise direction at which statistics are printed.
   The files contains 10 columns:
1) <img src="svgs/f709689c86153cfc47fc13bd71263621.svg?invert_in_darkmode&sanitize=true" align=middle width=69.11487pt height=24.56553pt/>, wall-normal coordinate, normalized with the local boundary layer thickness
2) <img src="svgs/413e5e559443639a2283994f9ccf0131.svg?invert_in_darkmode&sanitize=true" align=middle width=101.102265pt height=26.12412pt/> wall-normal coordinate, normalized with the viscous lenght scale
3) <img src="svgs/c003ef34c0f5583624f88cf9611fa556.svg?invert_in_darkmode&sanitize=true" align=middle width=76.46166pt height=26.12412pt/> the Favre averaged streamwise velocity, normalized with the friction velocity
4) <img src="svgs/5168dd0c22b7ee67557ebdc76a52a925.svg?invert_in_darkmode&sanitize=true" align=middle width=48.956655pt height=24.56553pt/> the streamwise velocity transformed according to van Driest
5) <img src="svgs/e387aa38fb0526c04f95bdf2108c30dd.svg?invert_in_darkmode&sanitize=true" align=middle width=33.612315pt height=26.12412pt/> the density scaled streamwise velocity rms
6) <img src="svgs/8c8bd17eadff3457884dec9db2dc1881.svg?invert_in_darkmode&sanitize=true" align=middle width=32.175495pt height=26.12412pt/> the density scaled wall-normal velocity rms, in viscous units
7) <img src="svgs/e8dfa6bfdb753c7e780b5fa7aee92971.svg?invert_in_darkmode&sanitize=true" align=middle width=35.96175pt height=26.12412pt/> the density scaled spanwise velocity rms, in viscous units  
8) <img src="svgs/7764c497204b291fd5e3218cb50d2a2a.svg?invert_in_darkmode&sanitize=true" align=middle width=27.95694pt height=26.12412pt/>      the density scaled Reynolds shear stress, in viscous units
9) <img src="svgs/b9943858c2a1621a8a2449fdec1bf1a9.svg?invert_in_darkmode&sanitize=true" align=middle width=52.296585pt height=29.42478pt/> the square root of the mean density, normalized with the wall density
10) The pressure rms, normalized with the square root of the wall-shear stress. 
