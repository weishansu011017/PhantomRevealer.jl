"""
    The tool kits for analysis/visualizing the dumpfile from PhantomRevealer analysis
        by Wei-Shan Su,
        July 18, 2024
"""
module_initialization()

"""
    Faceon_polar_plot(Disk2Ddata :: Analysis_result, array_index :: Int64,Log_flag::Bool=false,minzero::Bool=false, vlim :: Vector = [],colormap::String="plasma", figzise :: Tuple = (12,8))
Draw the face-on polar plot from the Faceon data.

#Parameters
- `Disk2Ddata :: Analysis_result`: The analysis result from PhantomRevealer
- `array_index :: Int64`: The column index of plotting array.
- `Log_flag :: Bool = true`: Change the colorbar to log scale
- `minzero::Bool=false`: Fix the minimum of colorbar to 0 if the automatical calculation of minimum is less then zero.
- `vlim :: Vector = []`: The colorbar range
- `colormap :: String = "plasma"`: The colormap.
- `figzise :: Tuple = (12,8)` The size of figure.

# Returns
- `PyCall.PyObject <pyplot_backend.polar_plot>`: The object of plotting.
"""
function Faceon_polar_plot(Disk2Ddata :: Analysis_result, array_index :: Int64,Log_flag::Bool=false,minzero::Bool=false, vlim :: Vector = [],colormap::String="plasma", figzise :: Tuple = (12,8))
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
    label_left = latexstring(L"$t = ",Int64(round(time)))
    label_right = Disk2Ddata.column_names[array_index]

    prplt = initialize_pyplot_backend()
    fax = prplt.polar_plot(s, ϕ, slabel, ϕlabel)
    fax.setup_fig(1,1,figzise)
    if !isempty(vlim)
        fax.pcolor_draw([z], [colormap],[z_unit],[Log_flag],[label_left],[minzero],[vlim])
    else
        fax.pcolor_draw([z], [colormap],[z_unit],[Log_flag],[label_left],[minzero])
    end
    fax.ax.test(0.7,0.95,label_right,transform=ax.transAxes,c="white", fontsize=14, verticalalignment="top", bbox=fax.props)
    return fax
end
