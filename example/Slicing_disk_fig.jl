using PhantomRevealer
initialize_modules()

"""
Generate the image from the dump files generated by `Slicing_disk.jl`
    by Wei-Shan Su,
    July 16, 2024
"""

function Constructing_figure(filepath::String)
    # ------------------------------PARAMETER SETTING------------------------------
    # Figure setting
    File_prefix :: String = "Slice" 
    Figure_format :: String = "png"
    figsize :: Tuple = (12,8)
    dpi = 450
    slabel = latexstring(L"$r$ [au]")
    zlabel = latexstring(L"$z$ [au]")
    colormap_gas :: String = "RdYlGn"
    colormap_dust :: String = "RdYlGn"
    clim_gas :: Vector = [1e-17,4e-14]
    clim_dust :: Vector = [1e-19,1e-15]
    streamline_density :: Float64 = 3.0
    streamline_color :: String = "black"
    # -----------------------------------------------------------------------------
    # Packaging parameters
    colormaps = [colormap_gas,colormap_dust]
    clims = [clim_gas, clim_dust]
    filename = splitext(filepath)[1]
    number_data = extract_number(filename)

    # Verifing data
    data :: Analysis_result = Read_HDF5(filepath)
    if data.params["Analysis_type"] != "Slicing_disk"
        error("InputError: The Analysis type of data needs to be `Slicing_disk`!")
    end
    transfer_cgs!(data)

    # Finding column_index
    target_columns :: Vector{String} = ["rho_g","rho_d","vs_g","vs_d","vz_g","vz_d"]
    target_column_index :: Vector{Int64} = zeros(Int64, length(target_columns))
    for key in keys(data.column_names)
        value = data.column_names[key]
        for (i,target) in enumerate(target_columns)
            if occursin(target, value)
                target_column_index[i] = key
            end
        end
    end
    rhog_index, rhod_index, vsg_index, vsd_index, vzg_index, vzd_index = target_column_index

    # Setup plotting target
    timestamp = data.time
    s = data.axes[1]
    z = data.axes[3]
    rhog_label = data.params["column_units"][rhog_index]
    rhod_label = data.params["column_units"][rhod_index]
    clabels :: Vector = [rhog_label, rhod_label]

    reduced_array = Vector{Array{Float64}}(undef, length(target_columns))
    for (i,index) in enumerate(target_column_index)
        reduced_array[i] = dropdims(mean(data.data_dict[index], dims = 2), dims = 2)'
    end

    rhog, rhod, vsg, vsd, vzg, vzd = reduced_array

    # Preparing plotting backend
    prplt = initialize_pyplot_backend()
    fax = prplt.cart_plot(s, z, slabel, zlabel)
    fax.setup_fig(2,1,figsize)

    fax.pcolor_draw([rhog,rhod], colormaps, clabels,[1,1],[latexstring(L"$t = ",Int64(round(timestamp)), L"$ yr")], clims, false)
    fax.streamplot_draw([vsg,vsd],[vzg,vzd],streamline_color, streamline_density, false)
    fax.set_ylabel(0)
    fax.set_ylabel(1)
    fax.set_xlabel(1)
    fax.save_fig("$(File_prefix)_$(number_data).$(Figure_format)",dpi)
    fax.close_fig()
end

function main()
    # Commendline variable setting
    if length(ARGS) < 1
        println("Usage: julia Slicing_disk_fig.jl <filename>")
        exit(1)
    end

    files = ARGS             

    First_logging()

    for file in files
        @info "File: $file"
        @time_and_print begin
            Constructing_figure(file)
        end 
    end

    @info "\nEnd analysis!"
end

main()