"""
    The tool kits for analysis/visualizing the dumpfile from PhantomRevealer analysis
        by Wei-Shan Su,
        July 18, 2024
"""
initialization_modules()
using Statistics
const TRANSFER_DICT = Dict{String, LaTeXString}(
    "∇" => L"$\nabla$",
    "ϕ" => L"$\phi$",
    "θ" => L"$\theta$",
    "ρ" => L"$\rho$",
    "Σ" => L"$\Sigma$",
)
"""
    replace_trans_LaTeXStr(str::String)
Replace all the weird charter to the LaTeX format in `TRANSFER_DICT`, and change the string into latexstring.

# Parameters
- `str :: String`: The string that would be transform.

# Returns
- `LaTeXString`: The transformation results.
"""
function replace_trans_LaTeXStr(str::String)
    transfer_dict = TRANSFER_DICT
    for target in keys(transfer_dict)
        result = transfer_dict[target]
        str = LaTeXString(replace(str, target=>result))
    end
    return str
end

"""
    Faceon_polar_plot(Disk2Ddata :: Analysis_result, array_index :: Int64,Log_flag::Bool=false,minzero::Bool=false, vlim :: Vector = [],colormap::String="plasma", figzise :: Tuple = (12,8))
Draw the face-on polar plot from the Faceon data.

#Parameters
- `Disk2Ddata :: Analysis_result`: The analysis result from PhantomRevealer
- `array_index :: Int64`: The column index of array.
- `Log_flag :: Bool = true`: Change the colorbar to log scale
- `vlim :: Union{Nothing,Vector} = nothing`: The colorbar range
- `colormap :: String = "plasma"`: The colormap.
- `figzise :: Tuple = (10,6)` The size of figure.

# Keyword argument
- `fax :: Union{PyObject, Nothing} = nothing`: The existing object of canvas. If `nothing` will generate a new object.

# Returns
- `PyCall.PyObject <pyplot_backend.polar_plot>`: The object of plotting.
"""
function Faceon_polar_plot(Disk2Ddata :: Analysis_result, array_index :: Int64,Log_flag::Bool=false, vlim :: Union{Nothing,Vector} = nothing,colormap::String="plasma", figzise :: Tuple = (8,6);fax::Union{PyObject, Nothing} =nothing)
    if Disk2Ddata.params["Analysis_type"] != "Faceon_disk"
        error("InputError: The Analysis type of data needs to be `Faceon_disk`!")
    end
    transfer_cgs!(Disk2Ddata)
    s = Disk2Ddata.axes[1]
    slabel = latexstring(L"$r$ [au]")
    ϕ = Disk2Ddata.axes[2]
    ϕlabel = latexstring(L"$\phi$")
    z = Disk2Ddata.data_dict[array_index]
    z_unit = Disk2Ddata.params["column_units"][array_index]
    time = Disk2Ddata.time
    label_left = latexstring(L"$t = ",Int64(round(time)), L"$")
    label_right = replace_trans_LaTeXStr(Disk2Ddata.column_names[array_index])

    if ϕ[end] != 2π
        ϕ = vcat(ϕ, 2π)
        z = hcat(z, z[:, 1])
    end

    prplt = initialize_pyplot_backend()
    if isnothing(fax)
        fax = prplt.polar_plot(s, ϕ, slabel, ϕlabel)
        fax.setup_fig(1,1,figzise)
    else
        fax.reset_fig()
    end
    fax.pcolor_draw([z], [colormap],[z_unit],[Log_flag],[label_left],vlim,false)
    fax.ax.text(0.9,0.95,label_right,transform=fax.ax.transAxes,c="white", fontsize=14, verticalalignment="top", bbox=fax.props)
    fax.draw_fig()
    return fax
end

"""
    Check_array_quantity(data :: Analysis_result, array_index :: Int64)
Print out the statistical properties of array with given array index.

# Parameters
- `data :: Analysis_result`: The analysis result from PhantomRevealer.
- `array_index :: Int64`: The column index of array.
"""
function Check_array_quantities(data :: Analysis_result, array_index :: Int64)
    try
        _ = data.params["_cgs"]
    catch err
        In_cgs = false
    else
        In_cgs = true
    end
    nanmax(itr) = maximum(x for x in itr if !isnan(x))
    nanmin(itr) = minimum(x for x in itr if !isnan(x))
    array = data.data_dict[array_index]
    column_name = data.column_names[array_index]
    shape = size(array)
    max = nanmax(array)
    min = nanmin(array)
    average = mean(array)
    med = median(array)
    STD = std(array)
    println("--------------Properties of array $column_name--------------")
    println("Size: $shape")
    println("In cgs unit?: $In_cgs")
    println("Maximum: $max")
    println("Minimum: $min")
    println("Average: $average")
    println("Median: $med")
    println("STD: $STD")
    println("----------------------------------------------------------------")
end

function spiral_detection(Disk2Ddata :: Analysis_result)
    function inverse_logarithmic_spiral(r, k_0, ϕmax, rmax)
        spiral = mod((log(r/rmax)/k_0)+ϕmax, 2π)
    end
    function k2pitch(k)
        beta = atan(abs(k))*(180/π)
    end
    error("This method have not done yet!")
end

function TESTHI()
    println("HI")
end