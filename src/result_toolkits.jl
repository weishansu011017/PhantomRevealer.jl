"""
    The tool kits for analysis/visualizing the dumpfile from PhantomRevealer analysis
        by Wei-Shan Su,
        July 18, 2024
"""
initialization_modules()

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
- `array_index :: Int64`: The column index of plotting array.
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
    if Disk2Ddata.params["Analysis_type"] != "Faceon"
        error("InputError: The Analysis type of data needs to be `Faceon`!")
    end
    transfer_cgs!(Disk2Ddata)
    s = Disk2Ddata.axis[1]
    slabel = latexstring(L"$r$ [au]")
    ϕ = Disk2Ddata.axis[2]
    ϕlabel = latexstring(L"$\phi$")
    z = Disk2Ddata.data_dict[array_index]
    z_unit = Disk2Ddata.params["column_units"][array_index]
    time = Disk2Ddata.time
    label_left = latexstring(L"$t = ",Int64(round(time)), L"$")
    label_right = replace_trans_LaTeXStr(Disk2Ddata.column_names[array_index])

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

function TESTHI()
    println("HI")
end