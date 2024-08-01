# PhantomRevealer

A Julia interface for analyzing dump files from the Phantom Smoothed Particle Hydrodynamics code.

Most of the analysis is based on SPH interpolation. Check out Price2010 for further information.

## Installation

### 1. Install Julia

Julia is required before installing this package. You can check this website: [Julia Downloads](https://julialang.org/downloads/) for further information.

After installing Julia, type 

```bash
julia
```
to start the julia REPL. 

### 2. Install package
In the Julia REPL, activate the Pkg module to install the package:

```julia
using Pkg
```

Next, install this package directly from the Git repository:
```julia
Pkg.add(url="https://github.com/weishansu011017/PhantomRevealer.git")
```

## Usage

### Read Phantom dumpfiles

To read the dumpfiles to the Julia environment, importing this package by typing

~~~julia
using PhantomRevealer
~~~

Next, loading your dumpfiles

~~~julia
filepath :: String = "/path/to/your/dumpfile"
prdf_list :: Vector{PhantomRevealerDataFrame} = read_phantom(filepath, "all")
~~~

You will get the vector that contains all of the particles in a Sarracen-alike dataframes structure, which is named as `PhantomRevealerDataFrame`,  with the separation of different type of particles. 

Note that the final element in this vector will always be the data of sinks particles (if exist.)

The structure of `PhantomRevealerDataFrame` is shown below:

~~~julia
struct PhantomRevealerDataFrame <: PhantomRevealerDataStructures
    dfdata::DataFrame
    params::Dict
end
~~~

The table of particles is stored in the `dfdata` field, and the other information of this simulation, including mass of particles, time...etc, can be found in the `params` field.

For example, to get the unit of mass in the code unit, type

~~~julia
umass :: Float64 = prdf_list[1].params["umass"]
~~~

to get the unit mass in cgs.

Just like sarracen, it is easy to add more quantity for each particles. **PhantomRevealer** provide some inbulid method. For instance, adding the position and velocity in cylindrical coordinate by entering

~~~julia
add_cylindrical!(prdf_list[1])
~~~

Note that **any fucntion/method name that includes `!` in PhantomRevealer mean that the function would change the inner state of first argument/input**. Becareful if you don't want to change the value, or just copy the variable by

~~~julia
Copy_data = deepcopy(prdf_list[1])
~~~



Beside, the coordinate transformation that shifts the origin to the positional of specific sinks particles is provied

~~~julia
sinks_data :: PhantomRevealerDataFrame = prdf_list[end]
COM2star!(prdf_list, sinks_data, 1)
~~~

This statement would change the whole coordinate into the new coordinate with the sinks particle in first row locating at the origin. The current coordinate can be checked by

~~~julia
prdf_list[1].params["COM_coordinate"]
prdf_list[1].params["Origin_sink_id"]
prdf_list[1].params["Origin_sink_mass"]
~~~

if the coordinate is in COM point of view. These value would be in 

~~~julia
>> [0.0, 0.0, 0.0]
>> -1
>> NaN
~~~

otherwise, they will give their coorespoding result.

### K-dimensional tree (KD-Tree) structure support

In some of the calculation, we would like to filter some of the particles to reduce the calculation. However, taking the geometry distance calculation directly will lead to a $\mathcal{O}(n^2)$ time complexity. Therefore, **PhantomRevealer** provides a method to filter your `PhantomRevealerDataFrame` which is based on K-dimensional tree structure that is provided by `Neighborhood.jl`(https://github.com/JuliaNeighbors/Neighborhood.jl). To start with, generate the kd-tree by following

~~~julia
kdtree3d :: KDTree = Generate_KDtree(data, dim=3)
kdtree2d :: KDTree = Generate_KDtree(data, dim=2)
~~~

Next, filtering the data by entering

~~~julia
target :: Vector = [10.0, π/4, 0.0]   # In cylindrical coordinate (s, ϕ, z) for exmaple.
radius :: Float64 = 2.5

filtered_data :: PhantomRevealerDataFrame = KDtree_filter(data, kdtree3d, target, radius, "polar" )  # You can change the final input to "cart" if the `target` input is in cartisian coordinate.
~~~

### Grid generator

We would like to look at the quantities at an arbitrary point in a grid; therefore, a grid generator is necessary. In **PhantomRevealer**, the grid structure is designed as 

~~~julia
struct gridbackend
    grid::Array
    axes::Vector{LinRange}
    dimension::Vector{Int64}
end
~~~

where `grid` field can store the value of quantities on each point, and the coordinate of each element is defined by the axes in `axes` field. The `dimension` field provies the size of this array.

To generate this grid automatically, entering

~~~julia
imin :: Vector{Float64} = [10.0, 0.0]     # The minimum value of each axes
imax :: Vector{Float64} = [120.0, 2π]			# The maximum value of each axes
in :: Vector{Float64} = [181, 301]				# The number of separation of each axes
type :: Type = Float64

test_grid :: gridbackend = generate_empty_grid(imin, imax, in ,type)
~~~

For a disk-shape grid, a useful function is provied

~~~julia
grid :: gridbackend = disk_2d_grid_generator(imin, iman, in)
~~~

the `grid` is equivalent to `test_grid`

### SPH interpolation

In SPH, an arbitrary quantity $\mathbf{A}(\mathbf{r})$ can be interpolated by the following formula.
```math
\mathbf{A}(\mathbf{r}) \approx \sum_b \frac{\mathbf{A}_b}{\rho_b} W(|\mathbf{r}-\mathbf{r}_b|; h)
```
Where $W$ is the kernel function, and the $h$ is the smoothed radius. **PhantomRevealer** provides 6 kernel functions: $M_4$ B-Spline, $M_5$ B-Spline, $M_6$ B-Spline, Wendland $C_2$, Wendland $C_4$ and Wendland $C_6$. The smoothed radius $h$​ for each point can be determine by three ways: taking the average of all of the particles `"mean"`, choosing the smoothed radius of closest particles to the target`"closest"`, calculating by interpolating the nearest particles `"intep"` .

### File extraction

The interpolation result needs to be extracted. A structure is designed to saving these result, extracting to a **HDF5** format.

~~~julia
struct Analysis_result_buffer <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,gridbackend}
    axes::Vector{LinRange}
    column_names::Dict{Int,String}
    params::Dict{String,Any}
end

mutable struct Analysis_result <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,Array{Float64}}
    axes::Dict{Int,Vector{Float64}}
    column_names::Dict{Int,String}
    params::Dict{String,Any}
end
~~~

You can transfer the `buffer` structure into the normal structure with some data testament by

~~~julia
result_buffer :: Analysis_result_buffer(SOMEPARAMETER...)  # Check out the section "Example 1" for further invesgation!
result :: Analysis_result = buffer2output(result_buffer)
~~~

After generateing the `Analysis_result` , the data can be written out by

~~~julia
Write_HDF5(filepath, result, "TEST")    # The `filepath` is used for extracting the step of simulation to prevent conflicting among different output e.g "disc_00100" -> "00100".
~~~

The `data_prefix` is the `PREFIX` part in the result string `PREFIX_00XXX.h5`

The output dump files can be loaded by `Read_HDF5()`.

~~~julia
result = Read_HDF5(filepath)
~~~



#### Example : Disk interpolation

Our goal is to interpolating a disk-shape (annulus) grid from a dumpfile `disc_00000` which including gaseous and dusty particles. We should also calculate the mid-plane average of density `rho` and velocity `vs, vϕ` . Here is a example to do it.

~~~julia
function Disk_Faceon_interpolation(filepath :: String)
    # ------------------------------PARAMETER SETTING------------------------------
    # parameters of radial axis
    smin :: Float64 = 10.0
    smax :: Float64 = 120.0
    sn :: Int64 = 141
  
    # parameters of azimuthal axis
    ϕmin :: Float64 = 0.0
    ϕmax :: Float64 = 2π
    ϕn :: Int64 = 301
  
    # Other parameters
    z_separation :: Int64 = 5												# The number of point selecting of vertical direction for taking average.
    midH_frac :: Float64 = 0.5 												# The ratio between the mid-plane disk scale height and the gaseous disk scale height
    column_names :: Vector = ["e"]									    	# The quantities that would be interpolate except for surface density `Sigma`.
    mid_column_names :: Vector = ["rho","vs", "vϕ"]                         # The quantities that would be interpolate in the midplane.
    center_sinks_id :: Int64 = 1											# The id of sink at the middle of disk for analysis.
    smoothed_kernel :: Function = M6_spline
    h_mode :: String = "closest"
    
    # Output setting
    File_prefix :: String = "Faceon"
    # -----------------------------------------------------------------------------
      # Packaging parameters
    sparams :: Tuple{Float64,Float64,Int} = (smin, smax, sn)
    ϕparams :: Tuple{Float64,Float64,Int} = (ϕmin, ϕmax, ϕn)
    columns_order :: Vector = ["Sigma", "∇Sigmas", "∇Sigmaϕ", column_names..., (mid_column_names.*"m")...] # construct a ordered column names (Those quantities with taking mid-plane average will have a suffix "m")
    
    # Load file
    prdf_list :: Vector{PhantomRevealerDataFrame} = read_phantom(filepath, "all")
    COM2star!(prdf_list, prdf_list[end], center_sinks_id)
    datag :: PhantomRevealerDataFrame = prdf_list[1]
    datad :: PhantomRevealerDataFrame = prdf_list[2]
    sinks_data :: PhantomRevealerDataFrame = prdf_list[3]
    
    # Add extra quantity for interpolation 
    add_cylindrical!(datag)
    add_cylindrical!(datad)
    add_eccentricity!(datag)
    add_eccentricity!(datad)
    
    # Make the `params` field
    time :: Float64 = get_time(datag)
    params :: Dict{String, Any} = Analysis_params_recording(datag, File_prefix)
    params["GasDiskMass"] = get_disk_mass(datag, sinks_data, smax, center_sinks_id)
    params["DustDiskMass"] = get_disk_mass(datad, sinks_data, smax, center_sinks_id)
  
    # Calculate the scale height of gaseous disk
    H_g = Disk_scale_height_analysis(datag, sparams) 
    
    # Interpolation
    grids_gas :: Dict{String, gridbackend} = Disk_2D_FaceOn_Grid_analysis(datag, sparams, ϕparams, column_names, mid_column_names, H_g, midH_frac, z_separation, smoothed_kernel, h_mode)
    grids_dust :: Dict{String, gridbackend} = Disk_2D_FaceOn_Grid_analysis(datad, sparams, ϕparams, column_names, mid_column_names, H_g, midH_frac, z_separation, smoothed_kernel, h_mode)
    
    # Combine these dictionaries of grids with suffix
    final_dict = create_grids_dict(["g","d"], [grids_gas, grids_dust])
    
    # Packaging the result
    Result_buffer :: Analysis_result_buffer = Analysis_result_buffer(time, final_dict, columns_order ,params)
    Result :: Analysis_result = buffer2output(Result_buffer)
    
    # Write out the result to HDF5
    Write_HDF5(filepath, Result, File_prefix)
    @info "-------------------------------------------------------"
end
~~~

### Built-in plotting backend

There has a built-in plotting backend in **PhantomRevealer** which is based on `matplotlib`. To use it, initializing the backend by

~~~julia
prplt = initialize_pyplot_backend()
~~~

The basic structure of backend is shown as the following

~~~python
class figure_ax:
    '''
    The plotting object in the code.
    
    Field in Class:
        fig: The figure/Canvas.
        ax: The axes of plotting/colorbar
        proj: The projection of plotting
    '''
    def __init__(self,proj=None):
        '''
        Two kinds of projection
        "None": Normal cartisian plotting
        "polar": Polar plotting
        '''
        self.fig: mfg.Figure = None
        self.ax: maxe._axes.Axes = None
        self.ncols = None
        self.nrows = None
        self.proj = proj
~~~

Follow the guideline to use the structure 

1. Generate an object

   ~~~julia
   fax = prplt.figure_ax()
   ~~~

2. Setup the desired figure

   ~~~julia
   fax.setup_fig(2,2,(10,6))
   ~~~

3. Make the plot at each axes

   ~~~julia
   fax.ax[0].plot(x,y,*plotting_params)
   cont = fax.ax[1].pcolor(x,y,z,*plotting_params)
   ~~~

   

4. Setup the colorbar for each map (if necessary). Allowing giving the colorbar a specific label so that the new colorbar can be drawn on the old ax rather then the new one.

   ~~~julia
   colorbar = fax.setup_colorbar(cont,"LABEL_OF_COLORBAR")
   ~~~

5. Save the figure

   ~~~julia
   save_fig("figure_1.png",450)
   ~~~

6. Close the figure

   ~~~julia
   fax.close_fig()
   ~~~

   Or, two subclasses for plotting is also provided 

   ~~~julia
   prplt.cart_plots()
   prplt.polar_plots()
   ~~~

   Check out the source code under `src/` to get further information.

## Relative links

Phantom homepage: [Phantom SPH](https://phantomsph.github.io)

Sarracen documentation: [Sarracen Documentation](https://sarracen.readthedocs.io/en/latest/)

## References

[1]:[Smoothed Particle Hydrodynamics and Magnetohydrodynamics ](https://ui.adsabs.harvard.edu/abs/2012JCoPh.231..759P/abstract) (Daniel J. Price, *Journal of Computational Physics, Volume 231, Issue 3, p. 759-794.* 2012)

