using PhantomRevealer
initialize_modules()

function main_analysis(filepath :: String)
    # ------------------------------PARAMETER SETTING------------------------------
    Analysis_tag_faceon :: String = "Faceon_disk"
    Analysis_tag_slice :: String = "Slicing_disk"
    # Disk generating setting (Base on cylindrical coordinate (s,ϕ,z))
    Origin_sinks_id = 1                                                     # The id of sink which is located at the middle of disk.
    ## s
    smin :: Float64 = 10.0                                                  # The minimum radius of disk.
    smax :: Float64  = 140.0                                                # The maximum radius of disk.
    sn :: Int64 = 131                                                       # The number of separation alone the radial direction on the disk.

    ## ϕ
    ϕmin :: Float64  = 0.0                                                  # The lower bound of angle for analysis.
    ϕmax :: Float64  = 2π                                                   # The upper bound of angle for analysis.
    ϕn :: Int64 = 24                                                        # The number of separation alone the azimuthal direction on the disk.
    
    ## z
    zmin :: Float64  = -33.0                                                # The lower bound of height from the z = 0 plane.
    zmax :: Float64  = 33.0                                                 # The upper bound of height from the z = 0 plane.
    zn :: Int64 = 80                                                        # The number of separation alone the vertical direction on the disk.
    
    column_names_faceon :: Vector{String} = ["e"]                    # The column that would be interpolated
    mid_column_names_faceon :: Vector{String} = ["rho","vs", "vϕ"]
    column_names_slice :: Vector{String} = ["vs", "vz"]                    # The column that would be interpolated
  
    # Other parameters
    z_separation :: Int64 = 5												# The number of point selecting of vertical direction for taking average.
    midH_frac :: Float64 = 0.5 												# The ratio between the mid-plane disk scale height and the gaseous disk scale height
    smoothed_kernel :: Function = M6_spline
    h_mode :: String = "closest"
    
    # Output setting (Faceon)
    File_prefix_faceon :: String = "Faceon"

    # Output setting (Edgeon)
    File_prefix_Slice :: String = "Slice"
    Save_figure :: Bool = true   
    Figure_format :: String = "png"

    # Figure setting
    figsize :: Tuple = (12,8)
    dpi = 450
    slabel = latexstring(L"$r$ [au]")
    zlabel = latexstring(L"$z$ [au]")
    colormap_gas :: String = "RdYlGn"
    colormap_dust :: String = "RdYlGn"
    clim_gas :: Vector = [1e-11,1e-7]
    clim_dust :: Vector = [1e-13,5e-8]
    streamline_density :: Float64 = 3.0
    streamline_color :: String = "black"
    # -----------------------------------------------------------------------------
    # Packaging parameters
    sparams :: Tuple{Float64,Float64,Int} = (smin, smax, sn)
    ϕparams :: Tuple{Float64,Float64,Int} = (ϕmin, ϕmax, ϕn)
    zparams :: Tuple{Float64,Float64,Int} = (zmin, zmax, zn)
    columns_order_faceon :: Vector = ["Sigma", "∇Sigmas", "∇Sigmaϕ", column_names_faceon..., (mid_column_names_faceon.*"m")...] # construct a ordered column names (Those quantities with taking mid-plane average will have a suffix "m")
    columns_order_slice :: Vector = ["rho", "∇rhos", "∇rhoϕ", column_names_slice...] # construct a ordered column names 
    colormaps = [colormap_gas, colormap_dust]
    clims = [clim_gas, clim_dust]

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

    # Get time stamp of the dumpfile
    time = get_time(datag)

    # Make the `params` field
    params_faceon :: Dict{String, Any} = Analysis_params_recording(datag, Analysis_tag_faceon)
    params_slice :: Dict{String, Any} = Analysis_params_recording(datag, Analysis_tag_slice)
    params_faceon["GasDiskMass"] = get_disk_mass(datag, sinks_data, smax, center_sinks_id)
    params_faceon["DustDiskMass"] = get_disk_mass(datad, sinks_data, smax, center_sinks_id)
    params_slice["GasDiskMass"] = get_disk_mass(datag, sinks_data, smax, center_sinks_id)
    params_slice["DustDiskMass"] = get_disk_mass(datad, sinks_data, smax, center_sinks_id)

    # Main analysis Slice
    grids_gas_slice :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datag, sparams, ϕparams, zparams, column_names_slice, smoothed_kernel, h_mode)
    grids_dust_slice :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datad, sparams, ϕparams, zparams, column_names_slice, smoothed_kernel, h_mode)
    
    # Calculate the scale height of gaseous disk
    H_g = Disk_scale_height_analysis(grids_dictg_slice) 

    # Main analysis Faceon
    grids_gas_faceon :: Dict{String, gridbackend} = Disk_2D_FaceOn_Grid_analysis(datag, sparams, ϕparams, column_names, mid_column_names, H_g, midH_frac, z_separation, smoothed_kernel, h_mode)
    grids_dust_faceon :: Dict{String, gridbackend} = Disk_2D_FaceOn_Grid_analysis(datad, sparams, ϕparams, column_names, mid_column_names, H_g, midH_frac, z_separation, smoothed_kernel, h_mode)

    # Combine these dictionaries of grids with suffix
    final_dict_slice = create_grids_dict(["g","d"], [grids_gas_slice, grids_dust_slice])
    final_dict_faceon = create_grids_dict(["g","d"], [grids_gas_faceon, grids_dust_faceon])

    Result_buffer_slice :: Analysis_result_buffer = Analysis_result_buffer(time, final_dict_slice, columns_order_slice ,params_slice)
    Result_buffer_faceon :: Analysis_result_buffer = Analysis_result_buffer(time, final_dict_faceon, columns_order_faceon ,params_faceon)

    # Packaging the result
    Result_slice :: Analysis_result = buffer2output(Result_buffer_slice)
    Result_faceon :: Analysis_result = buffer2output(Result_buffer_faceon)

    # Write out the results to HDF5
    Write_HDF5(filepath, Result_slice, File_prefix_slice)
    Write_HDF5(filepath, Result_faceon, File_prefix_faceon)
    @info "-------------------------------------------------------"
end

function main()
    # Commendline variable setting
    if length(ARGS) < 1
        println("Usage: julia Faceon_disk.jl <filename>")
        exit(1)
    end

    files = ARGS             

    First_logging()

    for file in files
        println("File: $file")
        @time_and_print begin
            main_analysis(file)
        end 
    end

    @info "\nEnd analysis!"
end

main()