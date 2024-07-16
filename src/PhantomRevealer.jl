module PhantomRevealer
# Include the Julia Module
# With the order of level
_module_location = @__DIR__
#Level 1 (Package and information)
include("$_module_location/module_initialization.jl")
include("$_module_location/logging.jl")
#Level 2 (SPH Mathematics)
include("$_module_location/mathematical_tools.jl")
include("$_module_location/kernel_function.jl")
#Level 3 (Data Structure)
include("$_module_location/grid.jl")
include("$_module_location/PhantomRevealerDataFrame.jl")
#Level 4 (Sigal point analysis and File read)
include("$_module_location/physical_quantity.jl")
include("$_module_location/read_phantom.jl")
#Level 5 (Analysis)
include("$_module_location/Grid_interpolation.jl")
#Level 6 (Extract data)
include("$_module_location/Extract_data.jl")
# # include("$_module_location/result_toolkits.jl")

# Import Python plt backend
push!(pyimport("sys")."path", _module_location)
prplt = pyimport("pyplot_backend")

# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end
