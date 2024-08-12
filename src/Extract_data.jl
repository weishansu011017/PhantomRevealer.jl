"""
Analysis result extraction.
    by Wei-Shan Su,
    July 11, 2024
"""

"""
    Analysis_result_buffer
The structure that contains all of the analysis result but not prepare for extraction yet. 

# Fields
- `time :: Float64`: The timestamp of simulation.
- `data_dict :: Dict{Int, gridbackend}`: The dictionary that contains all of the data.
- `axes :: Vector{LinRange}`: The axes of grid.
- `column_names :: Dict{Int,String}`: The column name of each data.
- `params :: Dict{String, Any}`: Other information of analysis/simulation.
"""
struct Analysis_result_buffer <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,gridbackend}
    axes::Vector{LinRange}
    column_names::Dict{Int,String}
    params::Dict{String,Any}
end

"""
    Analysis_result
The mutable structure that contains all of the analysis result for extraction.  

# Fields
- `time :: Float64`: The timestamp of simulation.
- `data_dict :: Dict{Int, Array{Float64}}`: The dictionary that contains all of the data.
- `axes :: Dict{Int, Vector{Float64}}`: The axes of grid.
- `column_names :: Dict{Int,String}`: The column name of each data.
- `params :: Dict{String, Any}`: Other information of analysis/simulation.
"""
mutable struct Analysis_result <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,Array{Float64}}
    axes::Dict{Int,Vector{Float64}}
    column_names::Dict{Int,String}
    params::Dict{String,Any}
end

"""
    extract_number(str::String)
Extract the lase group of number in a string.

# Parameters
- `str :: String`: The string that contains number.

# Returns
- `String`: The string of number.
"""
function extract_number(str::String)
    matches = findall(r"\d+", str)
    return !isempty(matches) ? str[matches[end].start:matches[end].stop] : ""
end

"""
    convert_field(value)
Convert the field in the `buffer` for extraction.

# Parameters
- `value`: The value in the `buffer`.

# Returns
The value in the converted type.
"""
function convert_field(value)
    if typeof(value) <: Dict{Int,gridbackend}
        converted_value = Dict{Int,Array{Float64}}()
        for key in keys(value)
            converted_value[key] = value[key].grid
        end
        return converted_value
    elseif typeof(value) <: Vector{LinRange}
        converted_value = Dict{Int,Vector{Float64}}()
        for i in eachindex(value)
            converted_value[i] = collect(value[i])
        end
        return converted_value
    else
        return value
    end
end

"""
    create_column_names(keys_order::Vector{String}, suffixes::Vector)
Create the array of column names in a interleaved order.

# Parameters
- `keys_order :: Vector{String}`: The sorted keys for all column.
- `suffixes :: Vector`: The suffixes of keys.

# Returns
- `Vector`: The array of column names.

# Example
```julia
keys = ["rho", "vs", "vϕ", "vz", "e"]
suffixes = ["g", "d"]
column_names = create_column_names(keys, suffixes)
println(column_names)  # Print out: ["phi", "rho_g", "rho_d", "vs_g", "vs_d", "vϕ_g", "vϕ_d", "vz_g", "vz_d", "e_g", "e_d"]
```
"""
function create_column_names(keys_order::Vector{String}, suffixes::Vector)
    column_names = ["phi"]
    for key in keys_order
        for suffix in suffixes
            push!(column_names, string(key, "_", suffix))
        end
    end
    return column_names
end

"""
    create_column_dict(column_names::Array{String})
Create the dictionary of column names in the `[0X     NAME]`-alike format.

# Parameters
- `column_names :: Array{String}`: The column name array.

# Returns
- `Dict{Int, String}`: The column name dictionary.

# Example
```julia
keys = ["rho", "vs", "vϕ", "vz", "e"]
suffixes = ["g", "d"]
column_names = create_column_names(keys, suffixes)
println(column_names)  # Print out: ["phi", "rho_g", "rho_d", "vs_g", "vs_d", "vϕ_g", "vϕ_d", "vz_g", "vz_d", "e_g", "e_d"]
column_dict = create_column_dict(column_names)
println(values(column_dict))  # Print out(non-ordered): ["[05        vs_d]", "[08        vz_g]", "[01         phi]", "[06        vϕ_g]", "[11         e_d]", "[09        vz_d]", "[03       rho_d]", "[07        vϕ_d]", "[04        vs_g]", "[02       rho_g]", "[10         e_g]"]
```
"""
function create_column_dict(column_names::Array{String})
    column_format = Dict{Int,String}()

    total_length = 16
    index_length = 2

    for (i, name) in enumerate(column_names)
        index_str = lpad(i, index_length, '0')
        space_padding = total_length - length(index_str) - length(name) - 3
        formatted_name = "[" * index_str * " " * repeat(" ", space_padding) * name * "]"
        column_format[i] = formatted_name
    end
    return column_format
end

"""
    create_grids_dict(suffixes :: Vector{String}, grids :: Vector{Dict{String, gridbackend}})
Create a dictionary that contains all of the data with different kinds of SPH particles.

# Parameters
- `suffixes :: Vector{String}`: The key of dict for the corresponding grid.
- `grids :: Vector{Dict{String, gridbackend}}`: The dictionary of data.

# Returns
- `OrderedDict{String, Dict{String, gridbackend}}`: The dictionary that contains all of the data with different kinds of SPH particles.
"""
function create_grids_dict(
    suffixes::Vector{String},
    grids::Vector{Dict{String,gridbackend}},
)
    if length(suffixes) != length(grids)
        error("DictGenerateError: Mismatching size of array `suffixes` and `grids`!")
    end
    grids_dict = OrderedDict{String,Dict{String,gridbackend}}()
    for i in eachindex(suffixes)
        grids_dict[suffixes[i]] = grids[i]
    end
    return grids_dict
end

"""
    Analysis_result_buffer(time :: Float64, grids_dict :: Dict{String, Dict{String, gridbackend}}, key_values::Vector{String})
Construct a `Analysis_result_buffer` structure with different suffixes

# Parameters
- `time :: Float64`: The timestamp of simulation.
- `grids_dict :: AbstractDict{String, Dict{String, gridbackend}}`: The dictionary that contains the analysis results in the format of `gridbackend` which is stored in the dictionaries. The keys would be the suffix of column names. 
- `keys_array :: Vector{String}`: The name of quantities for each column.

# Returns
- `Analysis_result_buffer`: The analysis data.
"""
function Analysis_result_buffer(
    time::Float64,
    grids_dict::AbstractDict{String,Dict{String,gridbackend}},
    keys_array::Vector{String},
    params::Dict{String,Any},
)
    function add_suffix_to_keys!(dict::Dict, suffix::String)
        new_dict = Dict{String,Any}()
        for (key, value) in dict
            new_key = string(key, suffix)
            new_dict[new_key] = value
        end
        empty!(dict)
        for (key, value) in new_dict
            dict[key] = value
        end
    end
    suffixes = collect(keys(grids_dict))
    for suffix in suffixes
        griddict = grids_dict[suffix]
        for i in eachindex(keys_array)
            key = keys_array[i]
            if !(haskey(griddict, key))
                if !(haskey(griddict, key * "m"))
                    error("ExtractionError: No match the column name $key !")
                else
                    keys_array[i] *= "m"
                end
            end
        end
    end
    column_names = create_column_names(keys_array, suffixes)
    column_dict = create_column_dict(column_names)
    sample_data = grids_dict[suffixes[1]][collect(keys(grids_dict[suffixes[1]]))[1]]
    axes = sample_data.axes

    prepare_dict = deepcopy(grids_dict)
    for suffix in suffixes
        griddict = prepare_dict[suffix]
        add_suffix_to_keys!(griddict, "_$suffix")
    end

    data_dict = Dict{Int,gridbackend}()
    for i in eachindex(column_names)
        column = column_names[i]
        if column == "phi"
            continue
        end
        for suffix in suffixes
            if haskey(prepare_dict[suffix], column)
                data_dict[i] = prepare_dict[suffix][column]
            end
        end
    end
    return Analysis_result_buffer(time, data_dict, axes, column_dict, params)
end

"""
    axes_self_check(input :: Analysis_result_buffer)
Checking whether the axes in the `data_dict` is equal to the `axes` field.

# Parameters
- `input :: Analysis_result_buffer`: The analysis buffer that would be used for checking.
"""
function axes_self_check(input::Analysis_result_buffer)
    input_data = input.data_dict
    iaxes = input.axes
    for key in keys(input.column_names)
        column = input.column_names[key]
        if occursin("phi", column)
            continue
        end
        if typeof(input_data[key]) <: gridbackend
            axes = input_data[key].axes
            for i in eachindex(iaxes)
                if !(iaxes[i] == axes[i])
                    error("AxesError: Mismatching of axes in $column")
                end
            end
        end
    end
    @info "Axes self-checking success! "
end

"""
    buffer2output(buffer_struct::Analysis_result_buffer)
Transfer the `Analysis_result_buffer` into `Analysis_result` for extracting the result to files.

# Parameters
- `buffer_struct :: Analysis_result_buffer`

# Returns
- `Analysis_result`
"""
function buffer2output(buffer_struct::Analysis_result_buffer)
    axes_self_check(buffer_struct)
    fields = fieldnames(Analysis_result_buffer)
    converted_values = [convert_field(getfield(buffer_struct, f)) for f in fields]
    return Analysis_result(converted_values...)
end

"""
    Write_HDF5(RawData_filename::String, data::Analysis_result, data_prefix :: String = "NaS")
Write the data into the HDF5 format. The output filename would be in `PPP_00XXX.h5`.

# Parameters
- `RawData_filename :: String`: The file name of the original dumpfile.
- `data :: Analysis_result`: The data that would be extracted.
- `data_prefix :: String = "NaS"`: The prefix of the output filename i.e. the `PPP` part. If this argument is not given manully, the value of `data.params["Analysis_type"]` would be selected.
"""
function Write_HDF5(
    RawData_filename::String,
    data::Analysis_result,
    data_prefix::String = "NaS",
)
    if data_prefix == "NaS"
        try
            data_prefix = data.params["Analysis_type"]
        catch
            data_prefix = ""
        end
    end
    numberfile = extract_number(RawData_filename)
    output_filename = data_prefix * "_" * numberfile * ".h5"
    h5open(output_filename, "w") do f
        for name in fieldnames(Analysis_result)
            val = getfield(data, name)
            if typeof(val) <: AbstractArray
                write(f, string(name), val)
            elseif typeof(val) <: Dict
                g = create_group(f, string(name))
                for (key, value) in val
                    write(g, string(key), value)
                end
            else
                write(f, string(name), val)
            end
        end
    end
end

"""
    Read_HDF5(filepath::String, close_warn::Bool=false)
Read the PhantomRevealer dumpfile into julia.

# Parameters
- `filepath :: String`: The filepath of dumpfile
- `close_warn :: Bool = false`: Disable the file version checking.

# Returns
- `Analysis_result`: The data in Analysis_result format
"""
function Read_HDF5(filepath::String, close_warn::Bool=false)
    function read_dict(f, name, K = String, V = Float64)
        dict = Dict{K,V}()
        g = f[name]
        for key in keys(g)
            val = read(g, key)
            if typeof(key) <: K
                dict[key] = val
            elseif (typeof(key) <: String) && (K <: Int)
                dict[parse(K, key)] = val
            else
                dict[convert(K, key)] = val
            end
        end
        return dict
    end
    function PR_version_match(file_identifier::String, close_warn::Bool=false)
        version = get_PhantomRevealer_version()
        pattern = r"PhantomRevealer:(\d+\.\d+\.\d+)"
        m = match(pattern, file_identifier)
        if m !== nothing
            file_version = String(m.captures[1])
            if version != file_version && !close_warn
                println("VersionMismatchWarning: The dump file was generated by PhantomRevealer version $file_version, but the loaded PhantomRevealer is version $version. Be cautious of potential unsupported features!")
            end
        else
            error("UnknownFileIdentifier: Is this a PhantomRevealer dumpfile?")
        end
    end
    data = nothing
    h5open(filepath, "r") do f
        time = read(f, "time")
        data_dict = read_dict(f, "data_dict", Int, Array)
        axes = read_dict(f, "axes", Int, Vector)
        column_names = read_dict(f, "column_names", Int, String)
        params = read_dict(f, "params", String, Any)
        PR_version_match(params["file_identifier"])
        data = Analysis_result(time, data_dict, axes, column_names, params)
    end
    return data
end

"""
    transfer_cgs!(data :: Analysis_result, year::Bool=true)
Transfer all the quantities into cgs unit, and also add another dictionary about the unit.

# Parameters
- `data :: Analysis_result`: The data.
- `year :: Bool = ture`: The flag about transfering `time` stamp into year or not.
- `au :: Bool = ture`: The flag about transfering distance into AU or not.

"""
function transfer_cgs!(data::Analysis_result, year::Bool = true)
    function replace_grident_exp!(latex_str)
        str = latex_str.s
        regex = r"cm\$\^\{(-?\d+)\}"
        m = match(regex, str)
        if m !== nothing
            exponent = parse(Int, m.captures[1]) - 1
            new_exponent_str = "cm\$^{$exponent}"
            new_str = replace(str, m.match => new_exponent_str)
            return LaTeXString(new_str)
        else
            return latex_str
        end
    end
    function extract_suffix(str::String)
        idx = findlast(isequal('_'), str)
        if idx === nothing
            return ""
        else
            result = str[idx+1:end]
            if endswith(result, "]")
                result = chop(result)
            end
            return result
        end
    end
    function extract_label(s::String)::String
        pattern = r"\[\d+\s+(.*?)\]"
        match_result = match(pattern, s)
        if match_result !== nothing
            return match_result.captures[1]
        else
            return ""
        end
    end
    
    umass = data.params["umass"]
    udist = data.params["udist"]
    utime = data.params["utime"]
    try
        _ = data.params["_cgs"]
    catch err
        usigma = umass / (udist^2)
        urho = umass / (udist^3)
        uv = udist / utime

        if year
            data.time *= (utime / 31536000)
            data.params["time"] *= (utime / 31536000)
        else
            data.time *= utime
            data.params["time"] *= utime
        end
        data.params["graindens"] *= urho
        data.params["grainsize"] *= udist


        column_unit = Dict{Int,LaTeXString}()
        for key in keys(data.column_names)
            column_name = data.column_names[key]
            suffix = extract_suffix(column_name)

            if occursin("sigma", column_name) || occursin("Sigma", column_name)
                data.data_dict[key] *= usigma
                header = LaTeXString(L"$\Sigma_")
                unit = LaTeXString(L"$ [g cm$^{-2}$]")
            elseif occursin("rho", column_name)
                data.data_dict[key] *= urho
                header = LaTeXString(L"$\rho_")
                unit = LaTeXString(L"$ [g cm$^{-3}$]")
            elseif occursin("vs", column_name)
                data.data_dict[key] *= uv
                header = LaTeXString(L"$v_{s,")
                suffix = LaTeXString(suffix* "}")
                unit = LaTeXString(L"$ [cm s$^{-1}$]")
            elseif occursin("vphi", column_name) || occursin("vϕ", column_name)
                data.data_dict[key] *= uv
                header = LaTeXString(L"$v_{\phi,")
                suffix = LaTeXString(suffix* "}")
                unit = LaTeXString(L"$ [cm s$^{-1}$]")
            elseif occursin("vz", column_name)
                data.data_dict[key] *= uv
                header = LaTeXString(L"$v_{z,")
                suffix = LaTeXString(suffix* "}")
                unit = LaTeXString(L"$ [cm s$^{-1}$]")
            else
                header = LaTeXString(extract_label(column_name))
                suffix = LaTeXString(L"")
                unit = LaTeXString(L"")
            end
            column_unit[key] = latexstring(header, suffix, unit)
            if occursin("∇", column_name)
                data.data_dict[key] /= udist
                column_unit[key] =
                    latexstring(L"$\nabla$", replace_grident_exp!(column_unit[key]))
            end
            data.params["column_units"] = deepcopy(column_unit)
        end
    end
    data.params["_cgs"] = true
end

"""
    add_more_label!(data::Analysis_result, column_index::Int, label::LaTeXString))
Add a extra label for a given column.

# Parameters
- `data :: Analysis_result`: The data.
- `column_index :: Int`: The index of column.
- `label :: LaTeXString`: The label.

"""
function add_more_label!(data::Analysis_result, column_index::Int, label::LaTeXString)
    if !haskey(data.params, "column_units")
        error("LookupError: You should calling the transfer_cgs!() function before calling this function!")
    else
        current_column_unit :: Dict = data.params["column_units"]
    end
    if haskey(current_column_unit, column_index)
        while true
            println(
                "Old column unit $(current_column_unit[column_index]) has found. Are you sure you want to replace it? [y/n]",
            )
            yn = String(readline())
            if yn == "y"
                break
            elseif yn == "n"
                return
            else
                nothing
            end
        end
    end
    data.params["column_units"][column_index] = label
end

"""
    add_dust2gas_ratio!(data::Analysis_result, column_index::Int64 = 61)
Add the column of dust-to-gas ratio

# Parameters
- `data :: Analysis_result`: The data.
- `column_index :: Int64 = 61`: The index of column.
"""
function add_dust2gas_ratio!(data::Analysis_result, column_index::Int64 = 61)
    transfer_cgs!(data)
    rhog :: Union{Nothing, Array{Float64}} = nothing
    rhod :: Union{Nothing, Array{Float64}} = nothing
    for key in keys(data.column_names)
        column_name = data.column_names[key]
        if occursin("rhom_g", column_name) || occursin("rho_g", column_name)
            rhog = data.data_dict[key]
        elseif occursin("rhom_d", column_name) || occursin("rho_d", column_name)
            rhod = data.data_dict[key]
        end
    end
    if isnothing(rhog) || isnothing(rhod)
        error("LookupError: The density column has not found!")
    end
    d2g = rhod./rhog
    data.column_names[column_index] = "[$column_index  dust-to-gas]"
    add_more_label!(data,column_index, L"\varepsilon")
    data.data_dict[column_index] = d2g
end

"""
    add_St!(data::Analysis_result, column_index::Int64 = 62)
Add the column of Stokes number.

# Parameters
- `data :: Analysis_result`: The data.
- `column_index :: Int64 = 61`: The index of column.
"""
function add_St!(data::Analysis_result, column_index::Int64 = 62)
    for column in values(data.column_names)
        if occursin("St",column)
            return
        end
    end
    transfer_cgs!(data)
    Sigmag :: Union{Nothing, Array{Float64}} = nothing
    for key in keys(data.column_names)
        column_name = data.column_names[key]
        if occursin("Sigmam_g", column_name) || occursin("Sigma_g", column_name)
            Sigmag = data.data_dict[key]
        end
    end
    if isnothing(Sigmag)
        error("LookupError: The surface density column has not found!")
    end
    grainsize = data.params["grainsize"]
    graindens = data.params["graindens"]
    St = (π/2)*(grainsize*graindens)./Sigmag
    data.column_names[column_index] = "[$column_index  St]"
    add_more_label!(data,column_index, L"St")
    data.data_dict[column_index] = St
end