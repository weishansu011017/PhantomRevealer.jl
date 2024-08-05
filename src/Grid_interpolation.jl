"""
The grid SPH interpolation
    by Wei-Shan Su,
    June 27, 2024

Generate a `gridbackend`, calculate the value for each column.

# Progression
Step 1. Generate a 3d KDTree to speed up the searching of particles

Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

Step 3. Calculate the result for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"
"""

"""
    Disk_3D_Grid_analysis(data::PhantomRevealerDataFrame ,s_params::Tuple{Float64,Float64,Int} ,ϕ_params :: Tuple{Float64,Float64,Int} ,z_params::Tuple{Float64,Float64,Int},column_names::Vector{String}, smoothed_kernal:: Function = M5_spline,h_mode::String="closest")
Calculate the SPH interpolation on a Edge-on grid that is described as a cylindrical coordinate (s,ϕ,z) for a disk.

Note: The grident values of arbitary quantities haven't supported yet!

The density `rho` and its grident vector `∇rho` would be calculated automatically. The `∇ρ` would be in cylindrical coordinate (∇ρs, ∇ρϕ, ∇ρz)

h_mode: 
"closest": Choose h of the closest particles.
"mean": use the mean value of h in the data.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
- `s_params :: Tuple{Float64,Float64,Int}`: The radial parameters with [smin, smax, sn]
- `ϕ_params :: Tuple{Float64,Float64,Int}`: The azimuthal parameters with [ϕmin, ϕmax, ϕn]
- `z_params :: Tuple{Float64,Float64,Int}`: The height parameters with [zmin, zmax, zn]
- `column_names :: Vector{String}`: The quantities that would be interpolated.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="closest"`: The mode for finding a proper smoothed radius. (Allowed value: "closest", "mean")

# Returns
- `Dict{String, gridbackend}`: The dictionary that contains all of the result by the form `gridbackend`.

# Examples
```julia
# Preparation
prdf_list :: Vector = read_phantom("dumpfile_00000", "all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
add_cylindrical!(data)
add_eccentricity!(data)
sparams :: Tuple{Float64,Float64,Int} = (10.0,100.0,91)
ϕparams :: Tuple{Float64,Float64,Int} = (0.0,2π,12)
zparams :: Tuple{Float64,Float64,Int} = (0.0,30.0,151)
smoothed_kernal :: Function = M6_spline
hmode :: String = "closest"

H_func = Disk_scale_height_analysis(data, sparams)

column_names :: Vector = [ "vs", "vϕ", "vz", "e"]
result1 :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(data, sparams,ϕparams, zparams , column_name, smoothed_kernal,hmode)
println(keys(result1)) # Print out ["Sigma", "∇Sigmas", "∇Sigmaϕ", "e", "rhom", "vsm", "vϕm"]
```
"""
function Disk_3D_Grid_analysis(
    data::PhantomRevealerDataFrame,
    s_params::Tuple{Float64,Float64,Int},
    ϕ_params::Tuple{Float64,Float64,Int},
    z_params::Tuple{Float64,Float64,Int},
    column_names::Vector{String},
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "closest"
)
    function wrap_dens(data::PhantomRevealerDataFrame, point::Array)::Float64
        return density(data, point, smoothed_kernal, h_mode, "polar")
    end
    function wrap_graddens(data::PhantomRevealerDataFrame, point::Array)::Vector
        return gradient_density(data, point, smoothed_kernal, h_mode, "polar")
    end
    function wrap_quant(data::PhantomRevealerDataFrame, point::Array)::Dict{String,Float64}
        return quantity_intepolate(
            data,
            point,
            column_names,
            smoothed_kernal,
            h_mode,
            "polar",
        )
    end
    @info "Start 3D disk grid analysis."
    # Add necessary quantities
    add_necessary_quantity!(data)

    # Checking data before interpolation
    if (data.params["Origin_sink_id"] == -1)
        error("IntepolateError: Wrong origin located!")
    end

    for column_name in column_names
        if !(hasproperty(data.dfdata, column_name))
            error("IntepolateError: Missing column name $column_name !")
        end
    end

    # Generate kd tree in 3D space
    kdtree3d = Generate_KDtree(data, 3)

    # Generate Edge-on grid 
    imin::Vector = [s_params[1], ϕ_params[1], z_params[1]]
    imax::Vector = [s_params[2], ϕ_params[2], z_params[2]]
    in::Vector = [s_params[3], ϕ_params[3], z_params[3]]
    empty_gridbackend::gridbackend = disk_3d_grid_generator(imin, imax, in)

    # Generate the coordinate array for the grid interpolation
    gridv::Array{Vector{Float64}} = generate_coordinate_grid(empty_gridbackend)

    # Preparation of result dictionary
    Result_dict = Dict{String,gridbackend}()
    Result_dict["rho"] = deepcopy(empty_gridbackend)
    Result_dict["∇rhos"] = deepcopy(empty_gridbackend)
    Result_dict["∇rhoϕ"] = deepcopy(empty_gridbackend)
    Result_dict["∇rhoz"] = deepcopy(empty_gridbackend)
    for column_name in column_names
        (column_name == "rho") && continue
        Result_dict[column_name] = deepcopy(empty_gridbackend)
    end

    # Prepare a roughly truncate radius for KD-tree filtering.
    roughly_truncated_radius::Float64 =
        get_truncated_radius(data, -1.0f0, 0.5, smoothed_kernal)

    # Iteration
    @threads for i in eachindex(gridv)
        target = gridv[i]
        kdtf_data = KDtree_filter(data, kdtree3d, target, roughly_truncated_radius, "polar") # New data that has been filtered.
        Result_dict["rho"].grid[i] = wrap_dens(kdtf_data, target)
        ∇dens = wrap_graddens(kdtf_data, target)
        Result_dict["∇rhos"].grid[i] = ∇dens[1]
        Result_dict["∇rhoϕ"].grid[i] = ∇dens[2]
        Result_dict["∇rhoz"].grid[i] = ∇dens[3]
        quantity_interpolation_dict::Dict{String,Float64} = wrap_quant(kdtf_data, target)
        if all(key -> haskey(Result_dict, key), keys(quantity_interpolation_dict))
            for key in keys(quantity_interpolation_dict)
                Result_dict[key].grid[i] = quantity_interpolation_dict[key]
            end
        else
            error("IntepolateError: Missing column name!")
        end
    end

    @info "End 3D disk grid analysis."
    return Result_dict
end

"""
    Disk_scale_height_analysis(edgeon_data_3D :: Dict{String, gridbackend})
Calculate the scale height of disk from existing edgeon data.

# Parameters
- `edgeon_data_3D :: Dict{String, gridbackend}`: The result from a edge-on analysis.

# Returns
- `Interpolations.Extrapolation`: The interpolation object of scale height as the function of `s`
"""
function Disk_scale_height_analysis(edgeon_data_3D :: Dict{String, gridbackend})
    rho_grid = grid_reduction(edgeon_data_3D["rho"], 2)
    s_separate = rho_grid.dimension[1]
    z_separate = rho_grid.dimension[2]
    s_array = rho_grid.axes[1]
    z_array = collect(rho_grid.axes[2])
    H_array = zeros(Float64, s_separate)
    @threads for i = 1:s_separate
        rho_array = rho_grid.grid[i, :]
        rho_sum = cumsum(rho_array) ./ maximum(rho_array)
        rho_cdf = rho_sum ./ maximum(rho_sum)
        for j = 1:z_separate
            cdfrho = rho_cdf[j]
            if cdfrho > 0.99999999
                H_array[i] = z_array[j]
                break
            end
        end
    end
    result = CubicSplineInterpolation(s_array, H_array, extrapolation_bc = Line())
    return result
end

"""
    Disk_scale_height_analysis(data::PhantomRevealerDataFrame, s_params::Tuple{Float64,Float64,Int}, ϕ_params :: Tuple{Float64,Float64,Int} = (0.0,2π,8) ,z_params::Tuple{Float64,Float64,Int} = (0.0, 28.0, 70),smoothed_kernal:: Function = M5_spline,h_mode::String="closest")
Calculate the scale height of disk.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
- `s_params :: Tuple{Float64,Float64,Int}`: The radial parameters with [smin, smax, sn]
- `ϕ_params :: Tuple{Float64,Float64,Int}`: The azimuthal parameters with [ϕmin, ϕmax, ϕn]
- `z_params :: Tuple{Float64,Float64,Int}`: The height parameters with [zmin, zmax, zn]
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="closest"`: The mode for finding a proper smoothed radius. (Allowed value: "closest", "mean")

# Returns
- `Interpolations.Extrapolation`: The interpolation object of scale height as the function of `s`
"""
function Disk_scale_height_analysis(
    data::PhantomRevealerDataFrame,
    s_params::Tuple{Float64,Float64,Int},
    ϕ_params::Tuple{Float64,Float64,Int} = (0.0, 2π, 8),
    z_params::Tuple{Float64,Float64,Int} = (0.0, 28.0, 70),
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "closest"
)
    edgeon_data_3D = Disk_3D_Grid_analysis(
        data,
        s_params,
        ϕ_params,
        z_params,
        Vector{String}(),
        smoothed_kernal,
        h_mode,
    )
    result = Disk_scale_height_analysis(edgeon_data_3D)
    return result
end

"""
    Disk_2D_FaceOn_Grid_analysis(data::PhantomRevealerDataFrame ,s_params::Tuple{Float64,Float64,Int} ,ϕ_params :: Tuple{Float64,Float64,Int}, 
    column_names::Vector{String}, mid_column_names::Vector{String}, 
    H_func::Interpolations.Extrapolation, midH_frac::Float64=0.5, midz_seperation::Int = 5,
    smoothed_kernal:: Function = M5_spline,h_mode::String="closest")
Calculate the SPH interpolation on a Face-on grid that is described as a polar coordinate (s,ϕ) for a disk.

Note: The grident values of arbitary quantities haven't supported yet!

The surface density `Sigma` and its grident vector `∇Sigma` would be calculated automatically. The `∇Sigma` would be in polar/cylindrical coordinate (∇Sigmas, ∇Sigmaϕ)

For those quantities that calculate by taking average in mid-plane would have a suffix `m`. e.g rhom, vsm,...

h_mode: 
"closest": Choose h of the closest particles.
"mean": use the mean value of h in the data.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
- `s_params :: Tuple{Float64,Float64,Int}`: The radial parameters with [smin, smax, sn]
- `ϕ_params :: Tuple{Float64,Float64,Int}`: The azimuthal parameters with [ϕmin, ϕmax, ϕn]
- `column_names :: Vector{String}`: The quantities that would be interpolated.
- `mid_column_names :: Vector{String}`: The quantities that would be evaluated by taking the midplane average.
- `H_func :: Interpolations.Extrapolation`: The interpolation object of scale height as the function of `s`
- `midH_frac :: Float64=0.5`: Fraction between The disk scale height and mid plane scale height i.e. Hmid/Hg
- `midz_seperation :: Int = 5`: The number of grid along the z axis for taking the average of midplane
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="closest"`: The mode for finding a proper smoothed radius. (Allowed value: "closest", "mean")

# Returns
- `Dict{String, gridbackend}`: The dictionary that contains all of the result by the form `gridbackend`.

# Examples
```julia
# Preparation
prdf_list :: Vector = read_phantom("dumpfile_00000", "all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
add_cylindrical!(data)
add_eccentricity!(data)
sparams :: Tuple{Float64,Float64,Int} = (10.0,100.0,91)
ϕparams :: Tuple{Float64,Float64,Int} = (0.0,2π,12)
midH_frac = 0.5
z_separate :: Int = 5
smoothed_kernal :: Function = M6_spline
hmode :: String = "closest"

H_func = Disk_scale_height_analysis(data, sparams)

column_names :: Vector = ["e"]
mid_column_names :: Vector = ["rho", "vs", "vϕ"]
result1 :: Dict{String, gridbackend} = Disk_2D_FaceOn_Grid_analysis(data, sparams, ϕparams, column_names, mid_column_names, H_func, midH_frac, z_separate, smoothed_kernal, hmode)
println(keys(result1)) # Print out ["Sigma", "∇Sigmas", "∇Sigmaϕ", "e", "rhom", "vsm", "vϕm"]
```
"""
function Disk_2D_FaceOn_Grid_analysis(
    data::PhantomRevealerDataFrame,
    s_params::Tuple{Float64,Float64,Int},
    ϕ_params::Tuple{Float64,Float64,Int},
    column_names::Vector{String},
    mid_column_names::Vector{String},
    H_func::Interpolations.Extrapolation,
    midH_frac::Float64 = 0.5,
    midz_seperation::Int = 5,
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "closest"
)
    function wrap_surf_dens(data::PhantomRevealerDataFrame, point::Array)::Float64
        return surface_density(data, point, smoothed_kernal, h_mode, "polar")
    end
    function wrap_grad_surf_dens(
        data::PhantomRevealerDataFrame,
        point::Array,
    )::Vector{Float64}
        return gradient_surface_density(data, point, smoothed_kernal, h_mode, "polar")
    end
    function wrap_quant(data::PhantomRevealerDataFrame, point::Array)::Dict{String,Float64}
        return quantity_intepolate(
            data,
            point,
            mid_column_names,
            smoothed_kernal,
            h_mode,
            "polar",
        )
    end
    function wrap_quant2D(
        data::PhantomRevealerDataFrame,
        point::Array,
    )::Dict{String,Float64}
        return quantity_intepolate_2D(
            data,
            point,
            column_names,
            smoothed_kernal,
            h_mode,
            "polar",
        )
    end
    function buffer_average_taker(buffer::Vector{Dict})
        buffer_keys = keys(buffer[1])
        buffer_length::Int64 = length(buffer)
        result_dict = Dict{String,Float64}()
        for key in buffer_keys
            value::Float64 = 0.0
            for i in eachindex(buffer)
                value += buffer[i][key]
            end
            value /= buffer_length
            result_dict[key] = value
        end
        return result_dict
    end
    @info "Start 2D disk grid analysis."
    # Add necessary quantities
    add_necessary_quantity!(data)

    # Checking data before interpolation
    if (data.params["Origin_sink_id"] == -1)
        error("IntepolateError: Wrong origin located!")
    end

    for column_name in column_names
        if !(hasproperty(data.dfdata, column_name))
            error("IntepolateError: Missing column name $column_name !")
        end
    end
    for column_name in mid_column_names
        if !(hasproperty(data.dfdata, column_name))
            error("IntepolateError: Missing column name $column_name !")
        end
    end

    # Generate kd tree
    kdtree3d = Generate_KDtree(data, 3)
    kdtree2d = Generate_KDtree(data, 2)

    # Generate Face-on grid 
    imin::Vector = [s_params[1], ϕ_params[1]]
    imax::Vector = [s_params[2], ϕ_params[2]]
    in::Vector = [s_params[3], ϕ_params[3]]
    empty_gridbackend::gridbackend = disk_2d_grid_generator(imin, imax, in)

    # Generate the coordinate array for the grid interpolation
    gridv::Array{Vector{Float64}} = generate_coordinate_grid(empty_gridbackend)

    # Preparation of result dictionary
    Result_dict = Dict{String,gridbackend}()
    Result_dict["Sigma"] = deepcopy(empty_gridbackend)
    Result_dict["∇Sigmas"] = deepcopy(empty_gridbackend)
    Result_dict["∇Sigmaϕ"] = deepcopy(empty_gridbackend)
    for column_name in column_names
        (column_name == "Sigma") && continue
        Result_dict[column_name] = deepcopy(empty_gridbackend)
    end
    imid_column_names = mid_column_names .* "m"
    for column_name in imid_column_names
        Result_dict[column_name] = deepcopy(empty_gridbackend)
    end

    # Prepare a roughly truncate radius for KD-tree filtering.
    roughly_truncated_radius::Float64 =
        get_truncated_radius(data, -1.0f0, 0.5, smoothed_kernal)

    # Iteration
    @threads for i in eachindex(gridv)
        target = gridv[i]
        # 2D intepolation
        kdtf_data2d =
            KDtree_filter(data, kdtree2d, target, roughly_truncated_radius, "polar") # New data that has been filtered.
        Result_dict["Sigma"].grid[i] = wrap_surf_dens(kdtf_data2d, target)
        ∇dens = wrap_grad_surf_dens(kdtf_data2d, target)
        Result_dict["∇Sigmas"].grid[i] = ∇dens[1]
        Result_dict["∇Sigmaϕ"].grid[i] = ∇dens[2]
        quantity_interpolation_dict::Dict{String,Float64} = wrap_quant2D(kdtf_data2d, target)
        if all(key -> haskey(Result_dict, key), keys(quantity_interpolation_dict))
            for key in keys(quantity_interpolation_dict)
                Result_dict[key].grid[i] = quantity_interpolation_dict[key]
            end
        else
            error("IntepolateError: Missing column name!")
        end
        # 3D midplane intepolation
        midH = H_func(target[1]) * midH_frac
        inteval = (-midH, midH)
        z_array = LinRange(inteval..., midz_seperation)
        buffer_array = Vector{Dict}(undef, midz_seperation)
        for k in eachindex(z_array)
            target3D = [target..., z_array[k]]
            kdtf_data3d = KDtree_filter(data, kdtree3d, target3D, roughly_truncated_radius, "polar") # New data that has been filtered.
            buffer_array[k] = wrap_quant(kdtf_data3d, target3D)
        end
        mid_interpolation_dict = buffer_average_taker(buffer_array)
        for j in eachindex(mid_column_names)
            Result_dict[imid_column_names[j]].grid[i] = mid_interpolation_dict[mid_column_names[j]]
        end
    end
    @info "End 2D disk grid analysis."
    return Result_dict
end
