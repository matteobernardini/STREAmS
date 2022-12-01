**A new version of the solver STREAmS-2.0 is now available https://github.com/STREAmS-CFD/STREAmS-2. 
We strongly suggest to use the new release.**












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
!=============================================================
```

`STREAmS` performs Direct Numerical Simulations of compressible turbulent flows in Cartesian geometry
solving the unsteady, fully compressible Navier-Stokes equations for a perfect gas.
Currently, three canonical wall-bounded flows can be simulated:
* compressible turbulent channel flow
* compressible zero-pressure-gradient turbulent boundary layer
* supersonic oblique shock-wave/turbulent boundary-layer interaction

STREAmS can be used on both local clusters and massively parallel HPC architectures, including those based on
Graphical Processing Units (GPUs). 

# References

Bernardini, M., Modesti, D., Salvadore, F., & Pirozzoli, S. (2021). STREAmS: A high-fidelity accelerated solver for direct numerical simulation of compressible turbulent flows. Computer Physics Communications, 263, 107906. https://doi.org/10.1016/j.cpc.2021.107906

# Compiling

STREAmS requires (1) a Fortran compiler and (2) an MPI library. For the GPU CUDA version, the NVIDIA compiler
is required (tested using PGI 19.4 or more recent compilers). A Makefile with 
predefined settings is available to facilitate compiling.
Different compiling modes can be selected by changing the three variables:

```
COMPILE 
MODE    
PREC    
```

`COMPILE` currently supports these choices:

* `pgi-cuda`: PGI compiler Cuda Fortran asynchronous version with MPI library (tested with OpenMPI provided by PGI)
* `pgi-cuda-sync`: PGI compiler Cuda Fortran synchronous version with MPI library (tested with OpenMPI provided by PGI)
* `pgi`: PGI compiler CPU version with OpenMPI library (tested with OpenMPI provided by PGI)
* `intel`: Intel compiler with MPI library (tested with IntelMPI library)
* `gnu`: gnu compiler with MPI library (tested with OpenMPI and MPICH)
* `ibmxl`: XL IBM compiler with MPI library (tested with IBM MPI library)
* `cray-cuda`: PGI compiler with Cray mpich library without support of CUDA-Aware MPI (currently oriented to Pitz-Daint cluster)

`MODE` can be one of:

* `opt`: optimized compilation, for standard production runs
* `debug`: debugging compilation, with run-time checks enabled and backtracing when available 

`PREC` defines the precision of the float numbers:

* `double`: double precision floating point numbers
* `single`: single precision floating point numbers

# Running

STREAmS can be usually executed using a standard MPI launcher, e.g. `mpirun`.
In addition to the executable, in the running folder you need:
* `input.dat`: file defining the physical and numerical setup of the simulation, to be customized
according to the desired needs. A detailed description of the input.dat file is given below.
Some examples of input.dat files are available in the `examples` folder.
* `database_bl.dat`: only required for initialization of `boundary layer` and `shock-boundary layer interaction` flows
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

`flow_type` defines the type of flow. 0 = channel flow, 1 = boundary layer, 2 = shock/boundary-layer interaction

`Lx Ly Lz` real numbers defining the size of the Cartesian domain along the three coordinate directions 
(x=streamwise, y=wall-normal, z=spanwise)

`Nx  Ny  Nz` integer numbers defining the number of grid nodes along each direction
 
 `Ny_wr     Ly_wr      dy+_w  jbgrid` specify the wall-normal mesh features. In all cases`dy+_w` is the desired spacing at the wall in inner units, <img src="svgs/f860ead74a4dff89d33b26468c06bbf7.svg?invert_in_darkmode&sanitize=true" align=middle width=32.43933pt height=26.17758pt/>.
 When `flow_type >0`, `Ny_wr` denotes the number of grid points in the wall-resolved region, ranging from y=0 (wall) up to y=`Ly_wr`, where a sinh mapping is applied. Both `Ny_wr` and `Ly_wr` can be specified when iflow = 1 (turbulent boundary layer), then a geometric progression is applied from y=`Ly_wr` up to y = Ly. When iflow = 2 (shock/boundary-layer interaction), `Ny_wr` must be specified but `Ly_wr` is automatically computed.

 `ng  visc_ord  ep_ord  weno_par`. Parameters to control the order of accuracy of the discretization. `ng` specifies the number of ghost
nodes. `visc_ord` represents the order of accuracy for the computation of the viscous terms (must be even <=6). `ep_ord` represents the order of accuracy for the computation of the convective terms in the smooth flow regions (central scheme, must be even <=6). `weno_par` selects the order of accuracy (order = 2*`weno_par`-1) for the computation of the convective terms in the shocked flow regions (WENO recontruction, must be <=3).

 `MPI_x_split  MPI_z_split` define the MPI decomposition along x (streamwise) and z (spanwise). These numbers 
 must be consistent to `Nx`, `Ny`, and `Nz`. In particular the following divisions must have zero remainder:
 `Nx/MPI_x_split`, `Nz/MPI_z_split`. Moreover, for `flow_type`>0 cases also `Nx/MPI_z_split` and `Ny/MPI_x_split`
 must have zero remainder.  The product `MPI_x_split * MPI_z_split` must equal the total number of launched MPI processes.

 `sensor_threshold   xshock_imp   deflec_shock`. `sensor_threshold` is the WENO threshold (WENO is active if the shock sensor value exceeds `sensor_threshold`), `xshock_imp` is the shock abscissa and `deflec_shock` is the shock angle (used only when `flow_type` > 0).

 `restart   num_iter   cfl   dt_control  print_control  io_type`. `restart` selects the start type (0=init, 1=restart run and statistics, 2=restart and continue statistics), `num_iter` the total number of iterations, `cfl` is the CFL number. The time step is re-evaluated every `dt_control` iterations. The residual file is printed every `print_control` iterations and `io_type` selects the type of I/O for restart (io_type=0 no I/O, io_type=1 serial I/O, io_type=2 MPI I/O).

 `Mach      Reynolds (friction)  temp_ratio   visc_type   Tref (dimensional)   turb_inflow`. `Mach` is the bulk Mach number for channel flow (defined with Twall) and the freestream Mach number for boundary layer and SBLI cases. `Reynolds` is the target friction Reynolds number. `temp_ratio` is the ratio between the wall temperature and the adiabatic wall temperature for boundary layer and SBLI, whereas for channel flow represents the ratio between bulk and wall temperature.  `visc_type` selects the viscosity law (visc_type=1 power law, visc_type=2 Sutherland law, with reference dimensional temperature specified by `Tref`). For boundary layers and SBLI, `turb_inflow` selects the turbulent inflow type (if > 0. implies digital filtering and its value is used to control/reduce temperature fluctuations, if < 0. implies recycling-rescaling method with constant spanwise shifting and its value denotes the location of the recycling station).  

 `stat_control  xstat_num`. Cumulative flow statistics are evaluate every `stat_control` iterations. `xstat_num` is the number of streamwise stations at which flow statistics are extracted and is only meaningful for boundary layer and SBLI cases (`flow_type` > 0). 

 `xstat_list` streamwise locations of boundary layer flow statistics (`flow_type` > 0)

 `dtsave dtsave_restart enable_plot3d enable_vtk`. `dt_save` is the time interval between field output, `dtsave_restart` is the time interval
 between output of restart files, `enable_plot3d`>0 activates the plot3d format output,
 `enable_vtk`>0 activates vtk format output.

  `rand_type` if < 0 produces not reproducible random sequences based on current time, if >=0 produces random sequences which
are reproducible across different runs with the same configuration.

Typical input files for canonical flow cases are available.

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
1) <img src="svgs/08b335933c05df68e8439a765e868cdf.svg?invert_in_darkmode&sanitize=true" align=middle width=26.33961pt height=24.6576pt/> is the wall-normal coordinate, normalized with the channel half width
2) <img src="svgs/9cbaafb1656e7b37560d2705ed2653e1.svg?invert_in_darkmode&sanitize=true" align=middle width=72.642735pt height=26.17758pt/> the wall-normal coordinate in viscous units
3) <img src="svgs/c29b0a6ba7b896aed4d1ff94a957d475.svg?invert_in_darkmode&sanitize=true" align=middle width=82.40859pt height=28.3107pt/> the wall-normal coordinate transformed according to Trettel & Larsson in viscous units
4) <img src="svgs/e1210f0fa24ba17ea9fd4da41c310dab.svg?invert_in_darkmode&sanitize=true" align=middle width=32.820645pt height=24.6576pt/> the mean streamwise velocity averaged according to Favre, normalized with the bulk flow velocity
5) <img src="svgs/edd39792d631fff9fc9340cbcd03af76.svg?invert_in_darkmode&sanitize=true" align=middle width=34.432035pt height=24.6576pt/> the mean streamwise velocity averaged according to Favre, normalized with the friction velocity
6) <img src="svgs/5168dd0c22b7ee67557ebdc76a52a925.svg?invert_in_darkmode&sanitize=true" align=middle width=49.08519pt height=24.6576pt/> the mean streamwise velocity transformed according to van Driest, normalized with the friction velocity
7) <img src="svgs/2fb4a896c47e22ceb0b1ad62e552a0bf.svg?invert_in_darkmode&sanitize=true" align=middle width=44.7876pt height=24.6576pt/> the mean streamwise velocity transformed according to Trettel & Larsson in viscous units
8) <img src="svgs/6d150b0d92a9b0f9f8ce9bf06e99af7c.svg?invert_in_darkmode&sanitize=true" align=middle width=35.03643pt height=24.6576pt/> the mean density profile, normalized with the mean wall density
9) <img src="svgs/9c48fce70048f82f658109eef23b7b71.svg?invert_in_darkmode&sanitize=true" align=middle width=69.34785pt height=34.18602pt/> the Favre streamwise Reynolds stress, normalized with the wall-shear stress
10) <img src="svgs/7c2b51045dc88620a88a255165e20708.svg?invert_in_darkmode&sanitize=true" align=middle width=67.642905pt height=34.18602pt/> the Favre wall-normal Reynolds stress, normalized with the wall-shear stress
11) <img src="svgs/f243b45c92e4fdc93b7e8280bcebb799.svg?invert_in_darkmode&sanitize=true" align=middle width=74.94894pt height=36.65145pt/> the Favre spanwise Reynolds stress, normalized with the wall-shear stress
12) <img src="svgs/b1254d1845c343e24eb7465dad8339d0.svg?invert_in_darkmode&sanitize=true" align=middle width=68.495295pt height=34.18602pt/> the Favre shear Reynolds stress, normalized with the wall-shear stress
13) <img src="svgs/41ea061742760dfbfe5de2c01f1c8bc1.svg?invert_in_darkmode&sanitize=true" align=middle width=39.534pt height=27.72594pt/>  the mean temperature profile, normalized with the wall temperature
14) <img src="svgs/6673acc30acf7858774f24361b9bcecf.svg?invert_in_darkmode&sanitize=true" align=middle width=46.857195pt height=30.4656pt/> the density fluctuations normalized with the mean wall density
15) <img src="svgs/6b014bcc811232a0bb81b6b3ab828c69.svg?invert_in_darkmode&sanitize=true" align=middle width=51.35493pt height=30.4656pt/> The temperature fluctuations, normalized with the wall temperature
## Boundary layer output files
### Residual file `output_streams.dat`
 When running boundary layer cases (`iflow`=1 or 2) the residuals file (output_streams.dat) contains three columns:
 1) number of cycles, 
 2) elapsed time, 
 3) streamwise momentum residual, 
### cf files
  The files cf_xxx.dat are ASCII files and the number xxx refers to the Cartesian streamwise MPI block to which the file belong.
  These files is printed at the end of each run. Statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1.
  The file contains 13 columns:
1) <img src="svgs/a7de9ea243c799172d54f6c73ca17aa2.svg?invert_in_darkmode&sanitize=true" align=middle width=31.47276pt height=24.6576pt/> streamwise coordinate normalized with the inlet boundary layer thickness
2) <img src="svgs/95b92a3b01f370498645efb59e52ff36.svg?invert_in_darkmode&sanitize=true" align=middle width=19.448715pt height=22.46574pt/> skin-friction coefficient
3) <img src="svgs/b91858acbc5918f19aa579c1c71a454a.svg?invert_in_darkmode&sanitize=true" align=middle width=19.26441pt height=22.83138pt/> friction Reynolds number 
4) <img src="svgs/7b9a0316a2fcd7f01cfd556eedf72e96.svg?invert_in_darkmode&sanitize=true" align=middle width=14.999985pt height=22.46574pt/> Compressible shape factor
5) <img src="svgs/b76f125299e13d61d8c9b80929221c5f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.316075pt height=22.46574pt/> Incompressible shape factor 
6) <img src="svgs/018562a8c71022027090ab072944bb8e.svg?invert_in_darkmode&sanitize=true" align=middle width=43.310685pt height=24.6576pt/> boundary layer thickness, normalized with the inlet boundary layer thickness
7) <img src="svgs/27e556cf3caa0673ac49a8f0de3c73ca.svg?invert_in_darkmode&sanitize=true" align=middle width=8.1735885pt height=22.83138pt/> Compressible displacement thickness
8) <img src="svgs/19889dd126af8c9b2d5d30c60cb5606d.svg?invert_in_darkmode&sanitize=true" align=middle width=26.368485pt height=22.83138pt/> Compressible momentum thickness
9) <img src="svgs/89144b2a5e5b3f3a3663daf61886c11b.svg?invert_in_darkmode&sanitize=true" align=middle width=16.80261pt height=14.15535pt/> friction velocity
10) Reynolds number based on the incompressible momentum thickness 
11) Incompressible friction coefficient, according to van Driest II transofrmation
12) Reynolds number based on the compressible momentum thickness
13) Mean pressure rms at the wall normalized with the wall-shear stress 
### stat files
   The files stat_nnn.dat are ASCII containing the mean boundary layer statistics.
   The files are printed at the end of each run. Statistics are progressively updated in time at each restart if idiski=2, or collected from scratch if idisk=1.
   The number nnn indicates the global mesh index in the streamwise direction at which statistics are printed.
   The files contains 10 columns:
1) <img src="svgs/06aaa96f13e0f03a6cde0227940dee9e.svg?invert_in_darkmode&sanitize=true" align=middle width=37.27944pt height=24.6576pt/>, wall-normal coordinate, normalized with the local boundary layer thickness
2) <img src="svgs/b93877c5651e7a90ae95501151ef0443.svg?invert_in_darkmode&sanitize=true" align=middle width=72.642735pt height=26.17758pt/> wall-normal coordinate, normalized with the viscous lenght scale
3) <img src="svgs/c003ef34c0f5583624f88cf9611fa556.svg?invert_in_darkmode&sanitize=true" align=middle width=76.67319pt height=26.17758pt/> the Favre averaged streamwise velocity, normalized with the friction velocity
4) <img src="svgs/5168dd0c22b7ee67557ebdc76a52a925.svg?invert_in_darkmode&sanitize=true" align=middle width=49.08519pt height=24.6576pt/> the streamwise velocity transformed according to van Driest
5) <img src="svgs/e387aa38fb0526c04f95bdf2108c30dd.svg?invert_in_darkmode&sanitize=true" align=middle width=33.737055pt height=26.17758pt/> the density scaled streamwise velocity rms
6) <img src="svgs/8c8bd17eadff3457884dec9db2dc1881.svg?invert_in_darkmode&sanitize=true" align=middle width=32.29479pt height=26.17758pt/> the density scaled wall-normal velocity rms, in viscous units
7) <img src="svgs/e8dfa6bfdb753c7e780b5fa7aee92971.svg?invert_in_darkmode&sanitize=true" align=middle width=36.095235pt height=26.17758pt/> the density scaled spanwise velocity rms, in viscous units  
8) <img src="svgs/7764c497204b291fd5e3218cb50d2a2a.svg?invert_in_darkmode&sanitize=true" align=middle width=28.05957pt height=26.17758pt/>      the density scaled Reynolds shear stress, in viscous units
9) <img src="svgs/b9943858c2a1621a8a2449fdec1bf1a9.svg?invert_in_darkmode&sanitize=true" align=middle width=52.296585pt height=29.42511pt/> the square root of the mean density, normalized with the wall density
10) The pressure rms, normalized with the wall-shear stress. 
