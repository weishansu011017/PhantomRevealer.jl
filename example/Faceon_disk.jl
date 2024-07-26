using PhantomRevealer

"""
Analysis the data in a face-on annulus grid.
    by Wei-Shan Su,
    July 17, 2024
"""

function Disk_Faceon_interpolation(filepath :: String)
    # ------------------------------PARAMETER SETTING------------------------------
    Analysis_tag :: String = "Faceon_disk"
    # parameters of radial axis
    smin :: Float64 = 10.0
    smax :: Float64 = 120.0
    sn :: Int64 = 141
    scale_height_sn :: Int64 = 75
  
    # parameters of azimuthal axis
    ϕmin :: Float64 = 0.0
    ϕmax :: Float64 = 2π
    ϕn :: Int64 = 301
  
    # Other parameters
    z_separation :: Int64 = 5												# The number of point selecting of vertical direction for taking average.
    midH_frac :: Float64 = 0.3 												# The ratio between the mid-plane disk scale height and the gaseous disk scale height
    column_names :: Vector = ["e"]									    	# The quantities that would be interpolate except for surface density `Sigma`.
    mid_column_names :: Vector = ["rho","vs", "vϕ"]                         # The quantities that would be interpolate in the midplane.
    Origin_sinks_id :: Int64 = 1											# The id of sink at the middle of disk for analysis.
    smoothed_kernel :: Function = M6_spline
    h_mode :: String = "closest"
    
    # Output setting
    File_prefix :: String = "Faceon"
    # -----------------------------------------------------------------------------
    # Packaging parameters
    sparams :: Tuple{Float64,Float64,Int} = (smin, smax, sn)
    scale_height_params :: Tuple{Float64,Float64,Int} = (smin, smax, scale_height_sn)
    ϕparams :: Tuple{Float64,Float64,Int} = (ϕmin, ϕmax, ϕn)
    columns_order :: Vector = ["Sigma", "∇Sigmas", "∇Sigmaϕ", column_names..., (mid_column_names.*"m")...] # construct a ordered column names (Those quantities with taking mid-plane average will have a suffix "m")
    
    # Load file
    prdf_list :: Vector{PhantomRevealerDataFrame} = read_phantom(filepath, "all")
    COM2star!(prdf_list, prdf_list[end], Origin_sinks_id)
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
    params :: Dict{String, Any} = Analysis_params_recording(datag, Analysis_tag)
    params["GasDiskMass"] = get_disk_mass(datag, sinks_data, smax, Origin_sinks_id)
    params["DustDiskMass"] = get_disk_mass(datad, sinks_data, smax, Origin_sinks_id)
    params["HmidOverH_g"] = midH_frac
  
    # Calculate the scale height of gaseous disk
    H_g = Disk_scale_height_analysis(datag, scale_height_params) 
    params["Hg(r)"] = H_g(LinRange(sparams...))

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

function main()
    # Commendline variable setting
    if length(ARGS) < 1
        println("Usage: julia Faceon_disk.jl <filename>")
        exit(1)
    end

    files = ARGS             

    First_logging()

    for file in files
        @info "File: $file"
        @time_and_print begin
            Disk_Faceon_interpolation(file)
        end 
    end

    @info "\nEnd analysis!"
end

main()