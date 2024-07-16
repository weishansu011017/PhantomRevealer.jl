"""
The single point SPH interpolation
    by Wei-Shan Su,
    June 27, 2024


# Structure:
    ## Preparation
        add_necessary_quantity()
    ## 3D interpolation
        ### density
        ### gradient density
        ### 3D quantity interpolation
    ## 2D interpolation
        ### surface density
        ### gradient surface density
        ### 2D quantity interpolation
"""
function add_necessary_quantity!(data::PhantomRevealerDataFrame)
    if !(hasproperty(data.dfdata, "rho"))
        add_rho!(data)
    end
    if !(hasproperty(data.dfdata, "Sigma"))
        add_Sigma!(data)
    end
end

function _standard_estimate_h_intepolate(dfdata::DataFrame, rnorm::Vector)
    """Assume data has been transfered to the reference_point-based coordinate"""
    fallback::Float32 = 0.0
    h_intepolate::Float32 = 0.0
    total_weight::Float32 = 0.0
    j = 0
    for particle in eachrow(dfdata)
        j += 1
        distance = rnorm[j]
        h_particle = particle.h
        weight = 1.0 / distance
        h_intepolate += weight * h_particle
        total_weight += weight
    end
    return total_weight > 0 ? h_intepolate / total_weight : fallback
end

function _easy_estimate_h_intepolate(dfdata::DataFrame, rnorm::Vector)
    """
    Give a specific smoothed radius to calculate
    Assume data has been transfered to the reference_point-based coordinate
    """
    if (isempty(rnorm))
        return NaN32
    end
    min_index = argmin(rnorm)
    return dfdata[min_index, "h"]
end

"""
    estimate_h_intepolate(data::PhantomRevealerDataFrame, rnorm::Vector,mode::String ="intep")
Provide a smoothed radius for interpolation.

mode: 
"intep": Intepolate h by local
"closest": Choose h of the closest particles.
"mean": use the mean value of h

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
- `rnorm :: Vector`: The distance between reference point and every particles. 
- `mode :: String ="intep"`: The mode of interpolation. (Allowed value: "intep", "closest", "mean")

# Returns
- `Float32`: The smoothed radius h.
"""
function estimate_h_intepolate(
    data::PhantomRevealerDataFrame,
    rnorm::Vector,
    mode::String = "intep",
)
    dfdata = data.dfdata
    if !(haskey(data.params, "h_mean"))
        add_mean_h!(data)
    end
    if (mode == "mean")
        return data.params["h_mean"]
    elseif (mode == "closest")
        return _easy_estimate_h_intepolate(dfdata, rnorm)
    elseif (mode == "intep")
        return _standard_estimate_h_intepolate(dfdata, rnorm)
    else
        error("IntepolateError: Invaild mode of h calculattion.")
    end
end

"""
    density(data::PhantomRevealerDataFrame, reference_point::Vector, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the density ρ at the given reference point by SPH interpolation

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Float64`: The value of density at the reference point.

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
reference_point :: Vector = [25.0, π/2, 0.0]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
dens = density(data, reference_point, smoothed_kernal, "closest", "polar")
```
"""
function density(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    reference_point is in "3D"
    "cart" = cartitian
    "polar" = cylindrical
    """
    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    rnorm = get_rnorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, rnorm, h_mode) # _easy_estimate_h_intepolate(dfdata, rnorm,)
    if h_intepolate == 0.0
        density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        filtered_rnorm = filter(r -> r <= truncate_radius, rnorm)
        density = sum(
            particle_mass .*
            Smoothed_kernel_function.(smoothed_kernal, h_intepolate, filtered_rnorm, 3),
        )
    end
    return density
end

"""
    gradient_density(data::PhantomRevealerDataFrame, reference_point::Vector, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the grident of density ∇ρ at the given reference point by SPH interpolation

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Vector`: The gradient vector of density at the reference point. The coordinate would be as same as the input(depends on "coordinate_flag").

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
reference_point :: Vector = [25.0, π/2, 0.0]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
∇dens = gradient_density(data, reference_point, smoothed_kernal, "closest", "polar")
```
"""
function gradient_density(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    rnorm, xyzref = get_r_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, rnorm, h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    grad_density::Vector = Vector{Float64}(undef, 3)
    if h_intepolate == 0.0
        fill!(grad_density, NaN)
        return grad_density
    else
        truncate_radius = truncate_multiplier * h_intepolate
        mask_rnorm = rnorm .< truncate_radius
        xyz_filtered = xyzref[mask_rnorm, :]
        buffer_array = zeros(Float64, size(xyz_filtered))
        for i = 1:size(xyz_filtered)[1]
            buffer_array[i, :] =
                particle_mass .* Smoothed_greident_kernel_function(
                    smoothed_kernal,
                    h_intepolate,
                    xyz_filtered[i, :],
                )
        end
        for j in eachindex(grad_density)
            grad_density[j] = sum(buffer_array[:, j])
        end

        if coordinate_flag == "polar"
            cylindrical_grad_density = _cart2cylin(grad_density)
            return cylindrical_grad_density
        else
            return grad_density
        end
    end
end

"""
    quantity_intepolate(data::PhantomRevealerDataFrame, reference_point::Vector, column_names::Vector{String}, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the value of every requested quantities at the given reference point by SPH interpolation.
Those points whose neighborhood has no particles around it would be label as `NaN`

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `column_names :: Vector{String}`: The quantities that would be interpolated.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Dict{String, Float64}`: A dictionary that contains the value of quantities at the reference point.

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
column_names :: Vector{String} = ["vr", "vϕ", "vz", "e", "St"]
reference_point :: Vector = [25.0, π/2, 0.0]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
quantities_dict :: Dict{Float64} = quantity_intepolate(data, reference_point, column_names, smoothed_kernal, "closest", "polar")
```
"""
function quantity_intepolate(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    column_names::Vector{String},
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    "cart" = cartitian
    "polar" = cylindrical
    """
    if !(hasproperty(data.dfdata, "rho"))
        add_rho!(data)
    end
    for column_name in column_names
        if !(hasproperty(data.dfdata, column_name))
            error("IntepolateError: No matching column '$(column_name)'.")
        end
    end
    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    divided_column = "rho"
    quantity_result = Dict{String,Float64}()
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    rnorm = get_rnorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, rnorm, h_mode)
    dfdata = data.dfdata

    if h_intepolate == 0.0
        for column_name in column_names
            quantity_result[column_name] = NaN
        end
    else
        truncate_radius = truncate_multiplier * h_intepolate
        indices = findall(x -> x <= truncate_radius, rnorm)
        filtered_dfdata = dfdata[indices, :]
        filtered_rnorm = filter(r -> r <= truncate_radius, rnorm)
        for column_name in column_names
            filtered_dfdata[!, column_name] ./= filtered_dfdata[!, divided_column]
            quantity_result[column_name] = sum(
                particle_mass .* (filtered_dfdata[!, column_name]) .* (
                    Smoothed_kernel_function.(
                        smoothed_kernal,
                        h_intepolate,
                        filtered_rnorm,
                        3,
                    )
                ),
            )
        end
    end
    return quantity_result
end

"""
    surface_density(data::PhantomRevealerDataFrame, reference_point::Vector, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the surface density Σ at the given reference point by SPH interpolation

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Float64`: The value of surface density at the reference point.

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
reference_point :: Vector = [25.0, π/2]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
surf_dens :: Vector = surface_density(data, reference_point, smoothed_kernal, "closest", "polar")
```
"""
function surface_density(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm = get_snorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, snorm, h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    if h_intepolate == 0.0
        surface_density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        filtered_snorm = filter(r -> r <= truncate_radius, snorm)
        surface_density = sum(
            particle_mass .*
            Smoothed_kernel_function.(smoothed_kernal, h_intepolate, filtered_snorm, 2),
        )
    end
    return surface_density
end

"""
    gradient_surface_density(data::PhantomRevealerDataFrame, reference_point::Vector, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the gredient of surface density ∇Σ at the given reference point by SPH interpolation

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Vector`: The gradient vector of surface density at the reference point. The coordinate would be as same as the input(depends on "coordinate_flag").

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
reference_point :: Vector = [25.0, π/2]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
grad_surf_dens :: Vector = grad_surface_density(data, reference_point, smoothed_kernal, "closest", "polar")
```
"""
function gradient_surface_density(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm, xyref = get_s_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, snorm, h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    grad_surface_density::Vector = Vector{Float64}(undef, 2)
    if h_intepolate == 0.0
        fill!(grad_surface_density, NaN)
        return grad_surface_density
    else
        truncate_radius = truncate_multiplier * h_intepolate
        mask_snorm = snorm .< truncate_radius
        xy_filtered = xyref[mask_snorm, :]
        buffer_array = zeros(Float64, size(xy_filtered))
        for i = 1:size(xy_filtered)[1]
            buffer_array[i, :] =
                particle_mass .* Smoothed_greident_kernel_function(
                    smoothed_kernal,
                    h_intepolate,
                    xy_filtered[i, :],
                )
        end
        for j in eachindex(grad_surface_density)
            grad_surface_density[j] = sum(buffer_array[:, j])
        end
        return grad_surface_density
    end
end

"""
    quantity_intepolate_2D(data::PhantomRevealerDataFrame, reference_point::Vector, column_names::Vector{String}, smoothed_kernal:: Function = M5_spline,h_mode::String="intep", coordinate_flag::String = "cart")
Calculate the value of every requested quantities at the given reference point by SPH interpolation ON THE 2D PLANE.
Those points whose neighborhood has no particles around it would be label as `NaN`

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame`
- `reference_point :: Vector`: The reference point for interpolation.
- `column_names :: Vector{String}`: The quantities that would be interpolated.
- `smoothed_kernal :: Function = M5_spline`: The Kernel function for interpolation.
- `h_mode :: String="intep"`: The mode for finding a proper smoothed radius. (Allowed value: "intep", "closest", "mean")
- `coordinate_flag :: String = "cart"`: The coordinate system that is used for given the target. Allowed value: ("cart", "polar") 

# Returns
- `Dict{String, Float64}`: A dictionary that contains the value of quantities at the reference point ON THE 2D PLANE.

# Example
```julia
prdf_list = read_phantom(dumpfile_00000,"all")
COM2star!(prdf_list, prdf_list[end],1)
data :: PhantomRevealerDataFrame = prdf_list[1]
column_names :: Vector{String} = ["vr", "vϕ", "e", "St"]
reference_point :: Vector = [25.0, π/2]  # Choosing the cylindrical coordinate
smoothed_kernal :: Function = M6_spline
quantities_dict :: Dict{Float64} = quantity_intepolate2D(data, reference_point, column_names, smoothed_kernal, "closest", "polar")
```
"""
function quantity_intepolate_2D(
    data::PhantomRevealerDataFrame,
    reference_point::Vector,
    column_names::Vector{String},
    smoothed_kernal::Function = M5_spline,
    h_mode::String = "intep",
    coordinate_flag::String = "cart",
)
    """
    Here recommended to use a single type of particle.
    coordinate_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if !(hasproperty(data.dfdata, "Sigma"))
        add_Sigma!(data)
    end
    for column_name in column_names
        if !(hasproperty(data.dfdata, column_name))
            error("IntepolateError: No matching column '$(column_name)'.")
        end
    end

    if coordinate_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    if get_dim(data) == 3
        divided_column = "Sigma"
    elseif get_dim(data) == 2
        divided_column = "rho"
    else
        error("IntepolateError: The dimension of data is not supported to this function!")
    end

    quantity_result = Dict{String,Float64}()
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm = get_snorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, snorm, h_mode)
    dfdata = data.dfdata


    if h_intepolate == 0.0
        for column_name in column_name
            quantity_result[column_name] = NaN
        end
    else
        truncate_radius = truncate_multiplier * h_intepolate
        indices = findall(x -> x <= truncate_radius, snorm)
        filtered_dfdata = dfdata[indices, :]
        filtered_snorm = filter(r -> r <= truncate_radius, snorm)
        for column_name in column_names
            filtered_dfdata[!, column_name] ./= filtered_dfdata[!, divided_column]
            quantity_result[column_name] = sum(
                particle_mass .* (filtered_dfdata[!, column_name]) .* (
                    Smoothed_kernel_function.(
                        smoothed_kernal,
                        h_intepolate,
                        filtered_snorm,
                        2,
                    )
                ),
            )
        end
    end
    return quantity_result
end