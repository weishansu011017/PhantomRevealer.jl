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
- `axis :: Vector{LinRange}`: The axis of grid.
- `column_names :: Dict{Int,String}`: The column name of each data.
- `params :: Dict{String, Any}`: Other information of analysis/simulation.
"""
struct Analysis_result_buffer <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,gridbackend}
    axis::Vector{LinRange}
    column_names::Dict{Int,String}
    params::Dict{String,Any}
end

"""
    Analysis_result
The mutable structure that contains all of the analysis result for extraction.  

# Fields
- `time :: Float64`: The timestamp of simulation.
- `data_dict :: Dict{Int, Array{Float64}}`: The dictionary that contains all of the data.
- `axis :: Dict{Int, Vector{Float64}}`: The axis of grid.
- `column_names :: Dict{Int,String}`: The column name of each data.
- `params :: Dict{String, Any}`: Other information of analysis/simulation.
"""
mutable struct Analysis_result <: PhantomRevealerDataStructures
    time::Float64
    data_dict::Dict{Int,Array{Float64}}
    axis::Dict{Int,Vector{Float64}}
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
    axis = sample_data.axis

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
    return Analysis_result_buffer(time, data_dict, axis, column_dict, params)
end

"""
    axis_self_check(input :: Analysis_result_buffer)
Checking whether the axes in the `data_dict` is equal to the `axis` field.

# Parameters
- `input :: Analysis_result_buffer`: The analysis buffer that would be used for checking.
"""
function axis_self_check(input::Analysis_result_buffer)
    input_data = input.data_dict
    iaxis = input.axis
    for key in keys(input.column_names)
        column = input.column_names[key]
        if occursin("phi", column)
            continue
        end
        if typeof(input_data[key]) <: gridbackend
            axis = input_data[key].axis
            for i in eachindex(iaxis)
                if !(iaxis[i] == axis[i])
                    error("AxisError: Mismatching of axis in $column")
                end
            end
        end
    end
    @info "Axis self-checking success! "
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
    axis_self_check(buffer_struct)
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

function Read_HDF5(filepath::String)
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
    data = nothing
    h5open(filepath, "r") do f
        time = read(f, "time")
        data_dict = read_dict(f, "data_dict", Int, Array)
        axis = read_dict(f, "axis", Int, Vector)
        column_names = read_dict(f, "column_names", Int, String)
        params = read_dict(f, "params", String, Any)
        data = Analysis_result(time, data_dict, axis, column_names, params)
    end
    return data
end

"""
    transfer_cgs!(data :: Analysis_result, year::Bool=true,au::Bool=true)
Transfer all the quantities into cgs unit, and also add another dictionary about the unit.

# Parameters
- `data :: Analysis_result`: The data.
- `year :: Bool = ture`: The flag about transfering `time` stamp into year or not.
- `au :: Bool = ture`: The flag about transfering distance into AU or not.

"""
function transfer_cgs!(data::Analysis_result, year::Bool = true; kwargs...)
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


        column_unit = Dict{Int,LaTeXString}()
        for key in keys(data.column_names)
            column_name = data.column_names[key]
            suffix = extract_suffix(column_name)

            if occursin("sigma", column_name) || occursin("Sigma", column_name)
                data.data_dict[key] *= usigma
                header = latexstring(L"$\Sigma_")
                unit = latexstring(L"$ [g cm$^{-2}$]")
            elseif occursin("rho", column_name)
                data.data_dict[key] *= urho
                header = latexstring(L"$\rho_")
                unit = latexstring(L"$ [g cm$^{-3}$]")
            elseif occursin("vs", column_name)
                data.data_dict[key] *= uv
                header = latexstring(L"$v_{r,")
                suffix = latexstring(suffix, L"}")
                unit = latexstring(L"$ [cm s$^{-1}$]")
            elseif occursin("vphi", column_name) || occursin("vϕ", column_name)
                data.data_dict[key] *= uv
                header = latexstring(L"$v_{\phi,")
                suffix = latexstring(suffix, L"}")
                unit = latexstring(L"$ [cm s$^{-1}$]")
            else
                header = latexstring(column_name)
                suffix = latexstring(L"")
                unit = latexstring(L"")
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

function add_more_label(data::Analysis_result, column_index::Int, label::LaTeXString)
    try
        _ = data.params["column_units"][column_index]

        if (data.params["column_units"][column_index] != "")
            while true
                println(
                    "Old column unit $(column_unit[column_index]) has found. Are you sure you want to replace it? [y/n]",
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
    catch err
        nothing
    end
    data.params["column_units"][column_index] = label
end
