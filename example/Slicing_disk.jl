using PhantomRevealer
using LaTeXStrings
using DataStructures
using PyCall

"""
Slice the disk for checking the edge-on vertical structure.
    Made by Wei-Shan Su, 2024
"""

function Slicing_disk(file::String)
    Analysis_tag :: String = "Slicing_disk"
    # General setting
    Smoothed_kernel_function :: Function = M6_spline                     # Allowed function: M4_spline, M5_spline, M6_spline, C2_Wendland, C4_Wendland, C6_Wendland
    h_mode :: String = "closest"                                         # Allowed mode: "mean", "closest", "intep"
    Rotate_xy_plane :: Bool = false                                      # Rotate the whole coordinate to the coordinate with z' axis paralleling to the direction of angular_momentum_vector of primary disc
    Extracting_dumpfile :: Bool = false                                  # Extract the analysis result.
    Save_figure :: Bool = true                                           # Saving the result as figure.

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

    # Output setting

    File_header :: String = "Slice"
    Figure_format :: String = "png"

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

    # Setup info
    initial_logging(get_analysis_info(file))

    # Packaging the parameters
    sparams :: Tuple{Float64,Float64,Int} = (smin, smax, sn)
    ϕparams :: Tuple{Float64,Float64,Int} = (ϕmin, ϕmax, ϕn)
    zparams :: Tuple{Float64,Float64,Int} = (zmin, zmax, zn)
    colormaps = [colormap_gas, colormap_dust]
    clims = [clim_gas, clim_dust]
    
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
    diskg_mass = get_disk_mass(datag, sinks_data, smax, Origin_sinks_id)
    diskd_mass = get_disk_mass(datad, sinks_data, smax, Origin_sinks_id)
    params["MassGaseousDisk"] = diskg_mass
    params["MassDustyDisk"] = diskd_mass

    # Add the cylindrical parameters
    add_cylindrical!(datag)
    add_cylindrical!(datad)

    # Main_analysis
    grids_dictg :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datag, sparams, ϕparams, zparams, column_names, Smoothed_kernel_function, h_mode)
    grids_dictd :: Dict{String, gridbackend} = Disk_3D_Grid_analysis(datad, sparams, ϕparams, zparams, column_names, Smoothed_kernel_function, h_mode)

    # Packaging the grids dictionary
    final_dict :: OrderedDict = create_grids_dict(["g","d"], [grids_dictg, grids_dictd])

    # Packaging the result
    Result_buffer :: Analysis_result_buffer = Analysis_result_buffer(time, final_dict, ["rho","∇rhos","∇rhoϕ","vs", "vz"],params)
    Result :: Analysis_result = buffer2output(Result_buffer)

    if Extracting_dumpfile
        Write_HDF5(file, Result, File_header)
    end

    if Save_figure
        # Initialize built-in plotting backend
        prplt = initialize_pyplot_backend()

        transfer_cgs!(Result)
        timestamp = Result.time
        s = Result.axis[1]
        z = Result.axis[3]
        rhog_label = Result.params["column_units"][2]
        rhod_label = Result.params["column_units"][3]
        clabels :: Vector = [rhog_label, rhod_label]
        fax = prplt.cart_plot(s,z, slabel, zlabel)
        fax.setup_fig(2,1,figsize)

        rhog = grid_reduction(grids_dictg["rho"],2).grid'
        rhod = grid_reduction(grids_dictd["rho"],2).grid'

        vsg = grid_reduction(grids_dictg["vs"],2).grid'
        vsd = grid_reduction(grids_dictd["vs"],2).grid'

        vzg = grid_reduction(grids_dictg["vz"],2).grid'
        vzd = grid_reduction(grids_dictd["vz"],2).grid'

        fax.pcolor_draw([rhog,rhod], colormaps, clabels,[1,1],[latexstring(L"$t = ",Int64(round(timestamp)), L"$ yr")],[1,1], clims, false)
        fax.streamplot_draw([vsg,vsd],[vzg,vzd],streamline_color, streamline_density, false)

        fax.set_ylabel(0)
        fax.set_ylabel(1)
        fax.set_xlabel(1)

        fax.save_fig("$(File_header)_$(extract_number(file)).$(Figure_format)",dpi)
    end

    @info "-------------------------------------------------------"
end

function main()
    # Commendline variable setting
    if length(ARGS) < 1
        println("Usage: julia Slicing_disk <filename>")
        exit(1)
    end

    files = ARGS             

    First_logging()

    for file in files
        println("File: $file")
        @time_and_print begin
            Slicing_disk(file)
        end 
    end

    @info "\nEnd analysis!"
end

main()