"""
Spiral detection in a Face-on data by using the ridge detection technique from Lindeberg(1996)(doi=10.1023/A:1008097225773)
    by Wei-Shan Su,
    Augest 30, 2024
"""

function image_gradient(image::Array{Float64,2}, axes::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}, Tuple{LinRange{Float64}, LinRange{Float64}}}=nothing)
    function compute_gradient_autodiff(func, x, y)
        g = ForwardDiff.gradient(z -> func(z[1], z[2]), [x, y])
        return g[1], g[2]
    end
    image_shape = size(image)
    if isnothing(axes)
        x = collect(LinRange(1,image_shape[1],image_shape[1]))
        y = collect(LinRange(1,image_shape[2],image_shape[2]))
        axes = (x,y)
    else
        axes_length = (length(axes[1]),length(axes[2]))
        if image_shape != axes_length
            error("ValueError: The length  of axes $axes_length should be equal to the shape of image $image_shape")
        end
    end
    image_func :: Interpolations.Extrapolation = LinearInterpolation(axes,image)

    grad_x = zeros(image_shape)
    grad_y = zeros(image_shape)

    for i in 1:image_shape[1]
        for j in 1:image_shape[2]
            x, y = axes[1][i], axes[2][j]
            grad_x[i, j], grad_y[i, j] = compute_gradient_autodiff(image_func, x, y)
        end
    end

    return grad_x, grad_y
end

function spiral_detection(Disk2Ddata :: Analysis_result,slim :: Tuple{Float64,Float64},image_subfix::String="_g", gaussian_sigma::Float64=1.0, draw::Bool = true;fax::Union{PyObject, Nothing} =nothing)
    if Disk2Ddata.params["Analysis_type"] != "Faceon_disk"
        error("InputError: The Analysis type of data needs to be `Faceon_disk`!")
    end
    radius_array = Disk2Ddata.axes[1]
    azimuth_array = Disk2Ddata.axes[2]
    slim_index = zeros(Int64,2)
    for i in eachindex(slim_index)
        slim_index[i] = value2closestvalueindex(radius_array,slim[i])
    end
    radius_array = radius_array[slim_index[1]:slim_index[2]]
    column_names = Disk2Ddata.column_names
    Sigma_index :: Union{Nothing, Int64} = nothing
    ∇Sigmas_index :: Union{Nothing, Int64} = nothing
    ∇Sigmaϕ_index :: Union{Nothing, Int64} = nothing
    for key in keys(column_names)
        column_name = column_names[key]
        if occursin("Sigma$image_subfix", column_name)
            Sigma_index = key
        elseif occursin("∇Sigmas$image_subfix", column_name)
            ∇Sigmas_index = key
        elseif occursin("∇Sigmaϕ$image_subfix", column_name)
            ∇Sigmaϕ_index = key
        end
    end
    if isnothing(Sigma_index)
        error("MissingColumn: Cannot find the column `Sigma`!")
    elseif isnothing(∇Sigmas_index)
        error("MissingColumn: Cannot find the column `∇Sigmas`!")
    elseif isnothing(∇Sigmaϕ_index)
        error("MissingColumn: Cannot find the column `∇Sigmaϕ`!")
    end
    RawSigma = Disk2Ddata.data_dict[Sigma_index][slim_index[1]:slim_index[2],:]
    RawSigmaave = mean(RawSigma)
    RawSigmastd = std(RawSigma)
    Sigma = imfilter(RawSigma, Kernel.gaussian(gaussian_sigma),Pad(:symmetric))
    Raw∇Sigmassph = Disk2Ddata.data_dict[∇Sigmas_index][slim_index[1]:slim_index[2],:]
    Raw∇Sigmaϕsph = Disk2Ddata.data_dict[∇Sigmaϕ_index][slim_index[1]:slim_index[2],:]
    ∇Sigmassph = imfilter(Raw∇Sigmassph, Kernel.gaussian(8.0),Pad(:symmetric))
    ∇Sigmaϕsph = imfilter(Raw∇Sigmaϕsph, Kernel.gaussian(8.0),Pad(:symmetric))
    ∇Sigmas, ∇Sigmaϕ = image_gradient(Sigma,(radius_array,azimuth_array))
    ∇∇Sigmass, ∇∇Sigmasϕ = image_gradient(∇Sigmas,(radius_array,azimuth_array))
    ∇∇Sigmaϕs, ∇∇Sigmaϕϕ = image_gradient(∇Sigmaϕ,(radius_array,azimuth_array))

    collected_points_radius :: Vector{Float64} = Vector{Float64}(undef, 0)
    collected_points_azimuth :: Vector{Float64} = Vector{Float64}(undef, 0)

    for (i, j) in product(eachindex(radius_array), eachindex(azimuth_array))
        if RawSigma[i,j] <= RawSigmaave + 1*RawSigmastd
            continue
        end
        raidus = radius_array[i]
        ∇SigSPHx = ∇Sigmassph[i,j]
        ∇SigSPHy = raidus * ∇Sigmaϕsph[i,j]
        Lxx = ∇∇Sigmass[i,j]
        Lxy = ∇∇Sigmasϕ[i,j]
        Lyy = ∇∇Sigmaϕϕ[i,j]
        Ldiagsub = Lxx-Lyy
        rocosβ = sqrt((1+((Ldiagsub)/sqrt(Ldiagsub^2 + 4*Lxy^2)))/2) 
        rosinβ = sign(Lxy)*sqrt((1-((Ldiagsub)/sqrt(Ldiagsub^2 + 4*Lxy^2)))/2)
        ∇SigSPHp = rosinβ*∇SigSPHx - rocosβ*∇SigSPHy
        Lpp = Lxx*rosinβ^2 - 2*Lxy*rosinβ*rocosβ + Lyy*rocosβ^2
        Lqq = Lxx*rocosβ^2 + 2*Lxy*rosinβ*rocosβ + Lyy*rosinβ^2
        if abs(∇SigSPHp)<=1e-14 && Lpp<=0.0 && abs(Lpp) >= abs(Lqq)
            push!(collected_points_radius, radius_array[i])
            push!(collected_points_azimuth, azimuth_array[j])
        end
    end
    collected_points = transpose(hcat(collected_points_azimuth,collected_points_radius))
    if draw
        fax = Faceon_polar_plot(fax=fax,Disk2Ddata,Sigma_index)
        for k in eachindex(collected_points_radius)
            fax.ax.scatter(collected_points[:,k]...,c="k",marker=".")
        end
        fax.draw_fig()
    end
    return fax, collected_points
end