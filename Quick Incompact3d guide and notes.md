This is a quick guide on to how to use **Incompact3d** and to know some details of its functioning. The present guide and the `modified` version of the code (baseline `v4.0`) are unofficial. Please refer to the official website of the code: [https://www.incompact3d.com](https://www.incompact3d.com/) and the major references [1],[2],[3],[4].

### Compiling
It is suggested to create a *build* directory in the same folder of the solver. 

<div style="text-align:center; font-family: monospace; background-color: #f2f2f2; padding: 10px;">mkdir build</div>

In the *build* directory run 
<div style="text-align:center; font-family: monospace; background-color: #f2f2f2; padding: 10px;">cmake ../</div>

After that, the *makefile* will be created. If a specific version of the Fortran compiler is needed, the option to be added is the following (with the example of `gfortran`):
<div style="text-align:center; font-family: monospace; background-color: #f2f2f2; padding: 10px;">cmake -DCMAKE_Fortran_COMPILER=gfortran ../</div>

It is now possibly to compile the binary file in the build directory with the following command:
<div style="text-align:center; font-family: monospace; background-color: #f2f2f2; padding: 10px;">make -j n</div>

where `n` is the number of processors that will be used to compile the program.
### Basic functioning
In this section, the main features of the code that can be useful to run simulations are reported.

In Incompact3d:

1) Non-dimensional Navier-Stokes equations are solved.
2) The $Re$ specified in the input file is used to compute the kinematic viscosity as $$\nu = \frac{1}{Re}$$so we are assuming unitary reference length and velocity scales. Reference scales depend then on the specific flow case.

	**Channel flow**
	For a channel flow, the default condition of simulation is to enforce a constant flow rate (CFR), while constant pressure gradient (CPG) can be enabled if necessary. Both conditions are rescaled with the centerline Reynolds number of a laminar Poiseuille flow,  $Re_0 = \frac{U_0 h}{\nu}$. A relation is available in order to estimate the related friction Reynolds number $Re_\tau$:
	$$Re_\tau \approx 0.116Re_0^{0.88}$$
	With CPG option enabled, the same relation is used to impose the same non-dimensionalization:
	$$Re_0 \approx \left( \frac{Re_\tau}{0.116} \right)^{1/0.88}$$
	In CFR conditions, being $U_0$ constant and unitary, this means that the bulk velocity $U_B$ must be kept to the value of $2/3$. This can be demonstrated being the laminar profile parabolic.  
	
	Finally, in order to estimate the related bulk Reynolds number $Re_B = \frac{2 U_B h}{\nu}$, the following relation can be employed (Pope, "Turbulent Flows"):
		$$ Re_\tau \approx 0.09Re_B^{0.88} \Rightarrow Re_B \approx \left( \frac{Re_\tau}{0.09} \right)^{1/0.88}$$
		
1) Pressure field is made non-dimensional with the freestream specific kinetic energy content $q_\infty = \frac{1}{2} \rho u_\infty^2$. The code performs operation on pressure field with a factor $\Delta t$, but it is rescaled before saving the snapshots with `rescale_pressure` subroutine.
2) To evaluate the stretching parameter for the mesh $\beta$, some useful scripts can be found in the folder `1-pre_processing/Stretching mesh`. The default code of Incompact3d is `stretching_parameter_channel.f90` and it is used to setup channel flow simulations. For temporal TBL simulations, the Python script `mesh_evaluation_ttbl.py` was developed. Low $\beta$ values correspond to a strong stretching of the mesh, while high $\beta$ values correspond to almost uniform mesh. 
3) Variables are saved on $y_p$ points, that are the faces of the grid elements. On the other hand, $y_{pi}$ points are the centers of the grid elements ($i:$ internal).
4) Boundary values of velocity are specified through the $b_{ijk}$ variables, where $i$ is the wall-normal direction of the boundary, $j$ is the direction of the specific velocity component and $k$ specifies if we are considering the bottom or the top walls (e.g. $b_{yx1}$ refers to the $x$ velocity component, specified at the bottom boundary with normal direction $y$).
5) Boundary conditions are specified in the input file `input.i3d` with the following variables: `nclx1, nclxn, ncly1, nclyn, nclz1, nclzn`, that specify the normal direction to the boundary and if we are considering the first or the last element along the specific direction. Values that can be adopted are: 0, for *periodic BC*, 1 for *free-slip BC* and 2 for *Dirichlet BC* (so imposed velocity value). Boundary conditions can be different along the same dimension, so different combinations can be enforced. 
6) Boundary conditions for scalar fields have the same functioning of the velocity BCs. They are called: `nclxS1, nclxSn, nclyS1, nclySn, nclzS1, nclzSn`.
7) The hyperviscous option for the second order derivative is a way to increase the numerical dissipation of the standard 6th order accurate scheme (that is sub-dissipative). In this manner, it is possible to increase the modified wavenumber at high wavenumbers, thus increasing the *dissipation error* $E_{diss}$, similarly to what is observed in high-order upwind schemes. This approach allows to prevent *wiggles* and thus to improve stability, even at high cell Reynolds number (or Péclet) $Pé \approx 200$. The two parameters available can be used also to have control on the numerical dissipation in the context of ILES simulations. 
8) Implicit time integration is available for the diffusive terms in $y$ direction. This approach is called **semi-implicit time integration**. From preliminary analyses, it appears that it allows to drop the restriction due to the stability parameter $S<1$ of fully explicit time integration schemes. The stability parameter $S$ can be calculated in the three directions
$$ S_{i} = \frac{u_i^2 \Delta t}{2 \nu}$$
		considering $i=x,y,z$.

11) Avoid to cancel the latest `restart.info` file, in order to being able to correctly calculate the current time unit (since the current time unit before restart is read from `restart.info` file itself).
12) The total shear velocity is calculated as: 
$$u_\tau = \sqrt{\nu \frac{\partial U_{\parallel,w}}{\partial y}}$$
	$U_\parallel = \sqrt{U^2 + W^2}$ and similarly for the $x$ and $z$ components: $u_{\tau,x} = \sqrt{\nu |\frac{\partial U_w}{\partial y}|}$ and $u_{\tau,z} = \sqrt{\nu |\frac{\partial W_w}{\partial y}|}$.
13) The friction coefficient $c_f$ is calculated as:
$$c_f = 2 \left( \frac{u_\tau}{U_{ref}} \right)^2$$
	where $U_{ref}$ is the reference velocity according to the specific flow case ($U_{w}$ for a TTBL and $U_B$ for a channel). 
###### References
1) [[High-order compact schemes for incompressible ﬂows_ A simple and efﬁcient method with quasi-spectral accuracy - Laizet & Lamballais - 2009.pdf|High-order compact schemes for incompressible ﬂows: A simple and efﬁcient method with quasi-spectral accuracy - Laizet & Lamballais - 2009]]
2) [[Incompact3d - A powerful tool to tackle turbulence problems with up to O(10**5 ) computational cores - Laizet & Li - 2011.pdf|Incompact3d - A powerful tool to tackle turbulence problems with up to O(10^5 ) computational cores - Laizet & Li - 2011]]
3) [[Xcompact3D - An open-source framework for solving turbulence problems on a Cartesian mesh - Bartholomew et al. - 2020.pdf|Xcompact3D - An open-source framework for solving turbulence problems on a Cartesian mesh - Bartholomew et al. - 2020]]
4) [[Straightforward high-order numerical dissipation via the viscous term for DNS and LES - Lamballais et al. - 2011.pdf|Straightforward high-order numerical dissipation via the viscous term for DNS and LES - Lamballais et al. - 2011]]
5) [[Turbulent Flows - Pope.pdf|Turbulent Flows - Pope]]

