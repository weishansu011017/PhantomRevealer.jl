using PhantomRevealer
initialize_modules()

"""
Slice the disk for checking the edge-on vertical structure.
    by Wei-Shan Su,
    July 16, 2024
"""

function Slicing_disk(file::String)
    @info "-------------------------------------------------------"
    # ------------------------------PARAMETER SETTING------------------------------
    Analysis_tag :: String = "Slicing_disk"
    # General setting
    smoothed_kernel :: Function = M6_spline                              # Allowed function: M4_spline, M5_spline, M6_spline, C2_Wendland, C4_Wendland, C6_Wendland
    h_mode :: String = "closest"                                         # Allowed mode: "mean", "closest"
    DiskMass_OuterRadius :: Float64 = 150.0                               # The outer radius of disk while estimating the mass of disk

    # Output setting
    File_prefix :: String = "Slice"
    HDF5 :: Bool = true                                                  # Extract the final result as HDF5 format
    figure :: Bool = false                                               # Extract the final result as figure

    # Disk generating setting (Base on cylindrical coordinate (s,ϕ,z))
    Origin_sinks_id = 1                                                  # The id of sink which is located at the middle of disk.
    ## s
    smin :: Float64 = 10.0                                               # The minimum radius of disk.
    smax :: Float64  = 140.0                                             # The maximum radius of disk.
    sn :: Int64 = 131                                                    # The number of separation alone the radial direction on the disk.

    ## ϕ
    ϕmin :: Float64  = 0.0                                               # The lower bound of angle for analysis.
    ϕmax :: Float64  = 2π                                                # The upper bound of angle for analysis.
    ϕn :: Int64 = 24                                                     # The number of separation alone the azimuthal direction on the disk.
    
    ## z
    zmin :: Float64  = -33.0                                             # The lower bound of height from the z = 0 plane.
    zmax :: Float64  = 33.0                                              # The upper bound of height from the z = 0 plane.
    zn :: Int64 = 80                                                     # The number of separation alone the vertical direction on the disk.
    
    column_names :: Vector{String} = ["vs", "vz"]                        # The column that would be interpolated

    # Figure setting (Valid only if figure == ture)
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
    # Setup info
    initial_logging(get_analysis_info(file))

    # Packaging the parameters
    sparams :: Tuple{Float64,Float64,Int} = (smin, smax, sn)
    ϕparams :: Tuple{Float64,Float64,Int} = (ϕmin, ϕmax, ϕn)
    zparams :: Tuple{Float64,Float64,Int} = (zmin, zmax, zn)
    columns_order :: Vector = ["rho", "∇rhos", "∇rhoϕ", column_names...] # construct a ordered column names 
    
    # Read file
    prdf_list :: Vector{PhantomRevealerDataFrame} = read_phantom(file,"all")
    sinks_data :: PhantomRevealerDataFrame = prdf_list[end]
    COM2star!(prdf_list, sinks_data ,Origin_sinks_id)
    datag :: PhantomRevealerDataFrame = prdf_list[1]
    datad :: PhantomRevealerDataFrame = prdf_list[2]
    
    # Get time stamp of the dumpfile
    time = get_time(datag)

    # Get params
    params = Analysis_params_recording(datag, Analysis_tag)
    diskg_mass = get_disk_mass(datag, sinks_data, DiskMass_OuterRadius , Origin_sinks_id)
    diskd_mass = get_disk_mass(datad, sinks_data, DiskMass_OuterRadius , Origin_sinks_id)
    params["MassGaseousDisk"] = diskg_mass
    params["MassDustyDisk"] = diskd_mass

    # Add the cylindrical parameters
    add_cylindrical!(datag)
    add_cylindrical!(datad)

    # Main_analysis
    grids_dictg :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datag, sparams, ϕparams, zparams, column_names, smoothed_kernel, h_mode)
    grids_dictd :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datad, sparams, ϕparams, zparams, column_names, smoothed_kernel, h_mode)

    # Packaging the grids dictionary
    final_dict :: OrderedDict = create_grids_dict(["g","d"], [grids_dictg, grids_dictd])

    # Packaging the result
    Result_buffer :: Analysis_result_buffer = Analysis_result_buffer(time, final_dict, columns_order,params)
    Result :: Analysis_result = buffer2output(Result_buffer)

    # Write the file to HDF5
    if HDF5
        Write_HDF5(file, Result, File_prefix)
    end

    # Construct the figure
    if figure
        # Packaging parameters
        colormaps = [colormap_gas,colormap_dust]
        clims = [clim_gas, clim_dust]
        number_data = extract_number(file)

        transfer_cgs!(Result)

        # Finding column_index
        target_columns :: Vector{String} = ["rho_g","rho_d","vs_g","vs_d","vz_g","vz_d"]
        target_column_index :: Vector{Int64} = zeros(Int64, length(target_columns))
        for key in keys(Result.column_names)
            value = Result.column_names[key]
            for (i,target) in enumerate(target_columns)
                if occursin(target, value)
                    target_column_index[i] = key
                end
            end
        end
        rhog_index, rhod_index, vsg_index, vsd_index, vzg_index, vzd_index = target_column_index

        # Setup plotting target
        timestamp = Result.time
        s = Result.axes[1]
        z = Result.axes[3]
        rhog_label = Result.params["column_units"][rhog_index]
        rhod_label = Result.params["column_units"][rhod_index]
        clabels :: Vector = [rhog_label, rhod_label]

        reduced_array = Vector{Array{Float64}}(undef, length(target_columns))
        for (i,index) in enumerate(target_column_index)
            reduced_array[i] = grid_reduction(Result.data_dict[index],2)'
        end

        rhog, rhod, vsg, vsd, vzg, vzd = reduced_array

        # Preparing plotting backend
        prplt = initialize_pyplot_backend()
        fax = prplt.cart_plot(s, z, slabel, zlabel)
        fax.__class__.anato_text_position = [0.01,0.96]
        fax.setup_fig(2,1,figsize)

        fax.pcolor_draw([rhog,rhod], colormaps, clabels,[1,1],[latexstring(L"$t = ",Int64(round(timestamp)), L"$ yr")], clims, false)
        fax.streamplot_draw([vsg,vsd],[vzg,vzd],streamline_color, streamline_density, false)
        fax.set_ylabel(0)
        fax.set_ylabel(1)
        fax.set_xlabel(1)
        fax.save_fig("$(File_prefix)_$(number_data).$(Figure_format)",dpi)
        fax.close_fig()
    end
    
    @info "-------------------------------------------------------"
end

function main()
    # Commendline variable setting
    if length(ARGS) < 1
        println("Usage: julia Slicing_disk.jl <filename>")
        exit(1)
    end

    files = ARGS             

    First_logging()

    for file in files
        @info "File: $file"
        @time_and_print begin
            Slicing_disk(file)
        end 
    end

    @info "\nEnd analysis!"
end

main()