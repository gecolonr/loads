cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
# Pkg.instantiate()
using PowerSystems
using PowerSimulationsDynamics
using ZIPE_loads
using TLmodels
using Combinatorics
using Sundials
using DifferentialEquations
using ArgCheck
using DataFrames
import Logging
using CSV
using Serialization

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

"""creates a copy of `injector` with the given name"""
function withName(injector::DynamicInjection, name::String)
    newinjector = deepcopy(injector)
    newinjector.name = name
    return newinjector
end

"""
makes all configurations of the given system and injectors, or specific given configurations.

## Args:
 - `sys::System` : the base system to work off of
 - `injectors` (array of `DynamicInjection`): injectors to use. if one dimensional, all configurations of the given injectors will be returned. If 2d, each row of `injectors` will represent output configuration.
 - `busgroups` (array of String or array of `Vector{String}`, optional): if just array of strings, represents list of buses to consider. If array of Vector{String}, each vector of bus names will be grouped together and always receive the same injector type. If not passed, just considers all busses with Generators attached.

## Returns:
 - `Dict{Vector{String}, System}`: dictionary of {generator names (ordered) => system}. contains all variations.
 """
function makeSystems(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)
    # first fix busgroups
    if isnothing(busgroups) # if nothing was passed
        # just get the name of every generator's bus in the system
        busgroups = [[get_bus(i).name] for i in get_components(Generator, sys)]
    elseif busgroups[1] isa String # if list of strings
        # make 2d. don't change data.
        busgroups = [[i] for i in busgroups]
    end # if 2d array, all work is done - that's our target form already.

    # now get groups of static generators to match busgroups
    gengroups = [[] for group in busgroups] # init empty array
    # iterate through preexisting generators
    for g in get_components(Generator, sys)
        # iterate through bus groups
        for (idx, group) in enumerate(busgroups)
            # if this generator's name is in this group,
            #       append this generator's name to this generator group
            if get_bus(g).name in group
                gengroups[idx] = [g.name, gengroups[idx]...]
                break
            end
        end
    end
    # now get injectors in correct form
    if length(size(injectors))>1 # this checks if it's 2d
        combos = (injectors[i, :] for i in 1:size(injectors)[1]) # if 2d, don't get permutations - user wants specific manually listed configurations
    else
        # otherwise, get all the permutations! length(busgroups) is the number of injectors in each configuration.
        combos = get_permutations(injectors, length(busgroups))
    end

    # now putting it all together

    # make a bunch of copies of `sys`
    systems = [deepcopy(sys) for _ in combos]
    # for each combination of injectors,
    for (comboIdx, combo) in enumerate(combos)
        # for each injector in this combination,
        #       make a new injector which shares a name with 
        #       each of the generators listed in the corresponding 
        #       generator group
        named_injectors = [
            withName(injector, gengroups[idx][jdx]) 
                for (idx, injector) in enumerate(combo) 
                    for jdx in 1:length(gengroups[idx]) 
                ]
        # now that we have all the generators, we can add them to the system at the static injector with the same name.
        for i in named_injectors
            add_component!(systems[comboIdx], i, first(get_components_by_name(Generator, systems[comboIdx], i.name)))
        end
    end

    # now compile all of the systems into a dictionary with the combos in a dictionary
    return Dict(zip([[i.name for i in combo] for combo in combos], systems))
end

"""
Struct to hold all the systems to iterate through.
Create a GridSearchSys with the constructor whose signature matches the [`makeSystems`](@ref) method. This will start you off with combinations of injectors. Then you can add a ZIPE load sweep with [`add_zipe_sweep`](@ref) or a line model sweep with [`add_lines_sweep`](@ref).
Use [`executeSims`](@ref) to run simulations for all the systems.

`length` and `size` are implemented to give the number of systems and (number of systems, number of columns in header).

## Attributes
`base` is the base system. It may be modified, but you should always be able to make any of the actual systems by simply adding to the base.
`header` is a list of the column titles for the dataframe. It should be human-readable.
`sysdict` stores the systems in the format config::Vector => System.
`results_header` is a list of column titles for results.
`results_getters` is a list of functions taking in a GridSearchSystem, Simulation, SmallSignalOutput, and String (error message) which returns results corresponding to `results_header`.
`df` is a DataFrame which will hold the results.
`chunksize` is the number of rows to hold in memory at a time (before saving to file).
"""
mutable struct GridSearchSys
    base::System
    header::Vector{String}
    sysdict::Dict{Vector{Any}, Function}
    results_header::Vector{String}
    results_getters::Vector{Function}
    df::DataFrame
    chunksize::Union{Int, Float64}
end

"""
WARNING: Generally, it's better to use Serialization-based storage, for speed and memory efficiency. This still works, but it's suboptimal.

Load data from Dataframe saved to TSV at `path`. stores this data in gss.df.
"""
function load_data!(gss::GridSearchSys, path::String)
    gss.df = load_data(path)
end

"""
helper function for `load_data`
tries to parse the string. if that fails, just return the string.
"""
function parse_result(item::String)
    try
        return eval(Meta.parse(item))
    catch
        return item
    end
end

"""
WARNING: Generally, it's better to use Serialization-based storage, for speed and memory efficiency. This still works, but it's suboptimal.

Load data from Dataframe saved to TSV at `path`.
"""
function load_data(path::String)
    println("Loading DataFrame...")
    df = CSV.read(path, DataFrame; delim='\t')
    println("DataFrame loaded!")
    for i in names(df)
        if i == "Eigenvectors"
            continue
        end
        println("Parsing column $i...")
        lk = ReentrantLock()
        Threads.@threads for j in 1:length(df[!, i])
            if !(df[j, i] isa String)
                continue
            end
            try
                x = eval(Meta.parse(df[j, i]))
                lock(lk) do 
                    df[j, i] = x
                end
            catch
                # don't do anything
            end
        end
        # df[!, i] .= map((x)->(x isa String ? parse_result(x) : x), df[:, i])
    end
    println("DataFrame loaded and parsed!")
    return df
end

"""
WARNING: This is very inefficient. Use Serialization with `save_serde_data` instead.
Save the data in `gss.df` to TSV.
"""
function save_data(gss::GridSearchSys, path::String)
    CSV.write(path, gss.df, delim='\t')
end
"""
Adds a column to the header. literally just `push!(gss.header, title)`. for cleanliness and encapsulation I guess.
"""
function add_column!(gss::GridSearchSys, title)
    push!(gss.header, title)
end

"""
Add a column to the output dataframe with a result to store.

`getter` must have the following signature:
(GridSearchSys, Simulation, SmallSignalOutput, String)->(Any)

Any of the inputs might be `missing`.
"""
function add_result!(gss::GridSearchSys, title::String, getter::Function)
    push!(gss.results_header, title)
    push!(gss.results_getters, getter)
end

"""
add multiple columns to the output dataframe with results to store. 

`getter` must have the following signature:
(GridSearchSys, Simulation, SmallSignalOutput, String)->(Vector{Any})

the output vector must be the same length as `titles` (the column titles), and any of the inputs might be `missing`.
"""
function add_result!(gss::GridSearchSys, titles::Vector{String}, getter::Function)
    push!(gss.results_header, titles...)
    push!(gss.results_getters, getter)
end
"""
set the number of rows save in each file (and thus how many to hold in memory before saving to file).

Set to `Inf` to hold all rows in memory (useful for small datasets and to allow use of the dataframe immediately after running the sims)
"""
function set_chunksize(gss::GridSearchSys, chunksize::Union{Int, Float64})
    @assert chunksize>=1.0 "invalid chunksize: can't save chunks of zero or negative size"
    gss.chunksize = chunksize
end

"""
constructor for GridSearchSys with the exact same behavior as [`makeSystems`](@ref).
"""
function GridSearchSys(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing, chunksize::Union{Int, Float64}=Inf)
    sysdict = makeSystems(sys, injectors, busgroups)
    newsysdict = Dict()
    for (key, val) in sysdict
        # instead of storing the `sys` object, we store a function which, when called, copies the system.
        # The external function allows us to save the pointer to the system (val) so that when we call this function
        # later, we'll get the value of val right now, not then.
        newsysdict[key] = ((value)->(()->deepcopy(value)))(val)
    end
    # first fix busgroups
    if isnothing(busgroups) # if nothing was passed
        # just get the name of every generator's bus in the system
        header = [get_bus(i).name for i in get_components(Generator, sys)]
    elseif busgroups[1] isa String # if list of strings
        header = busgroups
    else # if 2d array, all work is done - that's our target form already.
        header = [join(i, ", ") for i in busgroups]
    end
    header = (x->"injector at {$x}").(header)
    return GridSearchSys(sys, header, newsysdict, Vector(), Vector(), DataFrame(), chunksize)
end

"""
Adds a ZIPE load sweep to a GridSearchSys. Pass in a standard load to base the ZIPE load off of, and a vector of LoadParams structs to test.

The standard load should be a function which takes in a system and returns an appropriate load. if `missing`, won't add any standard loads.
"""
function add_zipe_sweep!(gss::GridSearchSys, standardLoadFunction::Union{Function, Missing}, zipe_params::Vector{LoadParams})
    sysdict = deepcopy(gss.sysdict)
    if !(standardLoadFunction isa Missing)
        for s in values(sysdict)
            # oh boy this is terrible code
            # something something monad side effect functional programming bind operator haskell lambda calculus
            
            # it's basically equivalent to this:
            #
            # s: &dyn Fn() -> System = make_loadadder_function(s: &dyn Fn() -> System)

            # where make_loadadder_function is defined like this:

            # function make_loadadder_function(system: &dyn Fn()->System) -> &dyn Fn(System)-> System
            #     function addload_noargs()
            #         sys = system() # instantiate system
            #         add_component!(sys, load_creator(sys))
            #         return sys
            #     end
            #     return addload_noargs
            # end

            s = (s->(()->(sys=s();add_component!(sys, standardLoadFunction(sys));sys)))(s)
        end
    end
    function withzipeload(sys, params)
        withzipe = deepcopy(sys)
        create_ZIPE_load(withzipe, params)
        return withzipe
    end
    gss.sysdict = Dict()
    for params in zipe_params
        for (key, val) in sysdict
            # withzipe = deepcopy(val)
            # create_ZIPE_load(withzipe, params)
            gss.sysdict[[key..., params]] = ((params, valarg)->(()->withzipeload(valarg(), params)))(params, val)
        end
    end
    add_column!(gss, "ZIPE Load Params")
end

"""
adds a sweep over line models and parameters.
`lineParams` should be a vector of LineModelParams structs, and `linemodelAdders` should be a dictionary of {model name (human readable) => func!(system, LineModelParams)} with all the model types to add.

This is based around the TLModels.jl package, so linemodelAdders could for example have the pair `"statpi" =>`[`create_statpi_system`](@ref)
"""
function add_lines_sweep!(gss::GridSearchSys, lineParams::Vector{LineModelParams}, linemodelAdders::Dict{String, Function}=Dict("statpi"=>create_statpi_system, "dynpi"=>create_dynpi_system, "mssb"=>create_MSSB_system))
    add_column!(gss, "Line Model")
    add_column!(gss, "Line Params")
    newsysdict = Dict()
    for (key, val) in gss.sysdict
        for (linename, creator) in linemodelAdders
            for params in lineParams
                newsysdict[[key..., linename, params]] = ((val, params, creator)->(()->creator(val(), params)))(val, params, creator)
            end
        end
    end
    gss.sysdict = newsysdict
end

import Base: length, size
"""
number of total systems in the gridsearch
"""
function length(gss::GridSearchSys)
    return length(gss.sysdict)
end

"""
(number of systems in gridsearch, number of columns in header)
"""
function size(gss::GridSearchSys)
    return (length(gss.sysdict), length(gss.header))
end

"""
run simulations on all of the systems in the gridsearch and store the results in a DataFrame.

## Args
 - `gss::GridSearchSys` : the systems
 - `change::Perturbation` : perturbation to apply to the system
 - `tspan::Tuple{Float64, Float64}` : time interval (in seconds) to simulate.
 - `dtmax::Float64` : max timestep for solver (make sure λh is in the feasible region for the solver)
 - `output_res::Float64` : resolution of saved timeseries results
 - `run_transient::Bool` : whether or not to run the transient simulations.
 - `log_path` : folder where outputs will be saved (when chunksize is reached). 

Whenever the number of rows in `gss.df` reaches `gss.chunksize`, results will be saved to file then deleted from gss.df in order to limit total memory usage.

If `gss.chunksize` is finite, the final dataframe will be saved. This way all results will be saved, even the last chunk or if `chunksize` was not reached.
In this case, the final chunk will also be deleted from gss.df.

To not save anything to file and keep the results in gss.df, make sure to `set_chunksize(gss, Inf)`.
"""
function executeSims!(gss::GridSearchSys, change::PSID.Perturbation, tspan::Tuple{Float64, Float64}=(0.48, 0.55), dtmax=0.0001, output_res=0.0001, run_transient::Bool=true, log_path::String="data/sims")
    # gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    # gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    # bus_voltages = (x->"Voltage at Bus $(get_number(x))").(get_components(Bus, gss.base))
    # inverter_currents = (x->"Inverter Current at Bus $(get_number(x))").(gen_busses)
    # generator_speeds = (x->"Generator Speed at Bus $(get_number(x))").(gen_busses)
    # cols = [gss.header..., "eigenvalues", "eigenvectors", bus_voltages..., inverter_currents..., generator_speeds..., "Simulation Status", "Simulation Time (ns)"]
    if !isdir(log_path)
        mkdir(log_path)
    end
    gss.df = DataFrame([i=>[] for i in [gss.header..., gss.results_header...]])
    counter = Threads.Atomic{Int}(0)
    total = length(gss)
    function inner(config::Vector{Any}, sys::System)
        (sim, sm, time, error) = runSim(sys, change, ResidualModel, tspan, IDA(; linear_solver=:Dense, max_convergence_failures=5), dtmax, output_res, run_transient)
        i = Threads.atomic_add!(counter, 1) + 1
        println("finished solve $i/$total in $(round(Int(time)/1e9, digits=2))s ($(round(100.0*i/total))%)")

        return vcat(config, reduce(vcat, (getter(gss, sim, sm, error) for getter in gss.results_getters)))
    end
    
    lk = ReentrantLock()
    
    chunk_counter = 0
    Threads.@threads for (key, val) in collect(gss.sysdict)
        results = inner(deepcopy(key), val())
        lock(lk) do 
            push!(gss.df, results)
            if size(gss.df, 1) >= gss.chunksize
                save_serde_data(gss, log_path*"/results$(chunk_counter).jls")
                deleteat!(gss.df, 1:size(gss.df, 1))
                chunk_counter += 1
            end
        end
    end
    if isfinite(gss.chunksize) && size(gss.df, 1)>0
        save_serde_data(gss, log_path*"/results$(chunk_counter).jls")
        deleteat!(gss.df, 1:size(gss.df, 1))
    end
end
"""
serializes the results dataframe and saves it to `path` to be read later using `Serialization.deserialize` (through `load_serde_data`)
"""
function save_serde_data(gss::GridSearchSys, path::String)
    if isfinite(gss.chunksize)
        @warn "chunksize is finite, so `gss.df` likely does not contain all results. Saving data currently contained in `gss.df`."
    end
    Serialization.serialize(path, gss.df)
end
"""
Loads serialized dataframe from file or folder.

If `path` is a folder, looks for all non-hidden .jls files, reads them, and concatenates them.
"""
function load_serde_data(path::String)
    if !isdir(path)
        return Serialization.deserialize(path)
    end

    files = [path*"/"*file for file in readdir(path) if (file[1]!='.')&&(file[end-3:end]==".jls")]
    dfs = Vector{DataFrame}(undef, length(files))
    counter = Threads.Atomic{Int}(0)
    progress_bar_width() = displaysize(stdout)[2]-28-Int(2*(floor(log10(length(files)))+1))-length(path)
    print("Reading files from $path: |"*(" "^progress_bar_width())*"| (0/$(length(files)))")
    Threads.@threads for i in 1:length(files)
        dfs[i] = load_serde_data(files[i])
        files_read = Threads.atomic_add!(counter, 1) + 1
        boxes = Int(round((files_read/length(files))*progress_bar_width()))
        print("\r"*(" "^(displaysize(stdout)[2])))
        print("\rReading files from $path: |"*("@"^boxes)*(" "^(progress_bar_width()-boxes))*"| ($files_read/$(length(files)))")
    end
    return vcat(dfs...)
end

"""
[This is a results getter function]

WARNING: mostly deprecated. just use `get_sim`.

Saves the `sim` object to file, then returns the filename to be added to the results dataframe.
"""
function get_serialized_sim_filename_function_builder(path::String)
    try
        mkdir(path)
    catch
        println("Failed to create directory $path (maybe it already exists)")
    end

    function get_serialized_sim_filename(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
        if sim isa Missing
            return missing
        end
        filename = "$(path)/SimObject$(time_ns()).jls"
        Serialization.serialize(filename, sim)
        return [filename]
    end
    return get_serialized_sim_filename
end

"""
[This is a results getter function]

gets system eigenvalues from small signal analysis.
"""
function get_eigenvalues(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    return sm isa Missing ? missing : [sm.eigenvalues]
end


"""
[This is a results getter function]

gets system eigenvectors from small signal analysis.
"""
function get_eigenvectors(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    return sm isa Missing ? missing : [sm.eigenvectors]
end

"""
[This is a results getter function]

**RETURNS A VECTOR!** use a vector of column titles if you want each bus to be a separate column.

gets voltage time series at every bus in the system.
"""
function get_bus_voltages(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Bus, gss.base)))
    else
        return [get_voltage_magnitude_series(sim.results, i)[2] for i in get_number.(get_components(Bus, gss.base))]
    end
end

"""
[This is a results getter function]

**RETURNS A VECTOR!** use a vector of column titles if you want each inverter to be a separate column.

Returns a vector with an entry for each `Generator` in the system. If there's an inverter there, the entry is a time series of the current magnitude. Otherwise, it's `missing`.
"""
function get_inverter_currents(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    throw("unimplemented: not working yet")
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Generator, gss.base)))
    end
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    inverters = [i for i in get_components(DynamicInverter, sys) if i.name in keys(gen_dict)]
    inverters = Dict(map(x->gen_dict[x], get_name.(inverters)) .=> inverters)
    return [(i in keys(inverters) ? current_magnitude(res, inverters[i].name) : missing) for i in get_number.(gen_busses)]
end


"""
[This is a results getter function]

**RETURNS A VECTOR!** use a vector of column titles if you want each generator to be a separate column.

gets speed (in rad/s) of all generators in the system.
"""
function get_generator_speeds(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    throw("unimplemented: not tested yet")
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Generator, gss.base)))
    end
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    generators = [i for i in get_components(DynamicGenerator, sys) if i.name in keys(gen_dict)]
    generators = Dict(map(x->gen_dict[x], get_name.(generators)) .=> generators)
    return [(i in keys(generators) ? generator_speed(res, generators[i].name) : missing) for i in get_number.(gen_busses)]
end

"""
[This is a results getter function]

**RETURNS A VECTOR!** use a vector of column titles if you want each load to be a separate column.

Gets the voltage magnitude time series at all ZIPE loads. 
"""
function get_zipe_load_voltages(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    if (sim isa Missing) || (sim.results isa Nothing)
        return Array{Missing}(missing, length(get_components(StandardLoad, gss.base)))
    end
    # the [2] here just gets the voltage instead of a tuple of (time, voltage) since we know the time divisions
    return [get_voltage_magnitude_series(sim.results, i)[2] for i in get_number.(get_bus.(get_components(StandardLoad, gss.base)))]
end

"""
[This is a results getter function]

**RETURNS A VECTOR!** use a vector of column titles if you want each inverter to be a separate column.

Gets the currrent magnitude time series at each ZIPE load.
"""
function get_zipe_load_current_magnitudes(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    throw("unimplemented: not tested yet.")
    if (sim isa Missing) || (sim.results isa Nothing)
        return Array{Missing}(missing, length(get_components(StandardLoad, gss.base)))
    end
    
    return ([current_magnitude(sim.results, "load_GFL_inverter"*string(get_number(get_bus(load)))).+ current_magnitude(sim.results, load.name) for load in get_components(StandardLoad, gss.base)])
end

"""
[This is a results getter function]

gets the value `sim.status` from the Simulation object (or missing).
"""
function get_sim_status(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    return (sim isa Missing) ? missing : sim.status
end

"""
[This is a results getter function]

gets the string representation of the error raised during simulation (or missing)
"""
function get_error(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    return error
end

"""
[This is a results getter function]

gets the whole Simulation object.
"""
function get_sim(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
    return sim
end

"""
extract array of current magnitude over time at the object with name `name` from the simulation results `res`.

Timestamps are thrown away. If you're using the GridSearchSystem, the returned array will correspond to an even time grid with spacing `dtmax` (the argument you passed to [`executeSims`](@ref), default 0.02s).
"""
function current_magnitude(res::SimulationResults, name::String)
    _, ir = get_state_series(res, (name, :ir_cnv))
    _, ii = get_state_series(res, (name, :ii_cnv))
    return sqrt.(ir.^2 .+ ii.^2)
end

"""
extracts generator speed time series from `res`. Time discretization is even `dtmax` (the argument you passed to [`executeSims`](@ref), default 0.02s)
"""
function generator_speed(res::SimulationResults, name::String)
    return get_state_series(res, (name, :ω))[2]
end

"""
expand line and load param columns into each individual parameter. This improves usability of the resulting dataframe, and makes saving to TSV cleaner.

Simply adds a column for each attribute of the line parameter struct and/or the load parameter struct.
"""
function expand_columns!(gss::GridSearchSys)
    columns_to_include = []
    if "Line Params" in names(gss.df)
        push!(columns_to_include, :"Line Params")
        for i in fieldnames(LineModelParams)
            gss.df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(gss.df.var"Line Params")
        end
    end

    if "ZIPE Load Params" in names(gss.df)
        push!(columns_to_include, :"ZIPE Load Params")
        for i in fieldnames(LoadParams)
            gss.df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(gss.df.var"ZIPE Load Params")
        end
    end

    select!(gss.df, Not(columns_to_include))
end
function expand_columns!(df::DataFrame)
    columns_to_include = []
    if "Line Params" in names(df)
        push!(columns_to_include, :"Line Params")
        for i in fieldnames(LineModelParams)
            df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"Line Params")
        end
    end

    if "ZIPE Load Params" in names(df)
        push!(columns_to_include, :"ZIPE Load Params")
        for i in fieldnames(LoadParams)
            df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"ZIPE Load Params")
        end
    end

    select!(df, Not(columns_to_include))
end

function unexpand_columns!(gss::GridSearchSys)
    # TODO: fill out this function
end

"""little wrapper to run simulations

## Args:
 - `system`: the system to simulate
 - `change`: perturbation to apply for the transient sim
 - `model`: model to pass into Simulation()
 - `tspan`: time span for transient simulation
 - `solver`: DE solver for transient simulation
 - `dtmax`: maximum timestep for DE solver
 - `output_res`: timestep of returned time series results
 - `run_transient`: whether or not to actually perform the transient simulation
 - `log_path`: path for simulation logs.

## Returns:
 - `(Simulation, SmallSignalOutput, UInt64, String)`: the Simulation object, the small signal analysis object, the time in nanoseconds the simulation took, and the error message. All but the time might be `missing` if things failed.
"""
function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(linear_solver=:LapackDense, max_convergence_failures=5), dtmax=0.001, output_res=0.02, run_transient=true, log_path::String=mktempdir())
    tic = Base.time_ns()
    local sim, sm

    sim = Simulation(
        model,
        system,
        log_path,
        tspan,
        change,
        disable_timer_outputs=true
        
    )
    try
        sm = small_signal_analysis(sim)
    catch err
        return (sim, missing, Base.time_ns()-tic, "Small Signal Analysis failed with error $err")
    end
    if run_transient
        try
            execute!(
                sim,
                solver,
                dtmax = dtmax,
                saveat = output_res,
                tstops = [0.5],
                # enable_progress_bar = false, # with multithreading it's meaningless anyways
            )
        catch err
            return (sim, sm, Base.time_ns()-tic, "Transient sim failed with error $err")
        end
    end

    return (sim, sm, Base.time_ns()-tic, missing)
end

"""
gets all permutations of `iterable` of length `k` with replacement.

## Args
- `iterable` : some iterable (with well-defined length and indexing)
- `k::Int` : desired permutation length

## Returns
- `Base.Generator`: Generator of vectors of length `k` of combinations of elements of `iterable`

## Examples
```jldoctest
julia> collect(get_permutations([1, 2], 2))
4-element Vector{Vector{Int64}}:
 [1, 1]
 [1, 2]
 [2, 1]
 [2, 2]
```
"""
function get_permutations(iterable, k)
    # yes this is really the best way I could think of.
    # computes all k-digit numbers in base-n, and then indexes `iterable` with the digits.
    n = length(iterable)
    return ([iterable[parse(Int, j, base=n)+1] for j in string(i, base=n, pad=k)] for i in 0:(n^k-1))
end

