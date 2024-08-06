# cd(@__DIR__)
# cd("..")
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()
# using PowerSystems
# using PowerSimulationsDynamics
# using PowerFlows
# using ZIPE_loads
# using TLmodels
# using Combinatorics
# using Sundials
# using DifferentialEquations
# using ArgCheck
# using DataFrames
# import Logging
# using CSV
# using Serialization
# using InfrastructureSystems

# const PSY = PowerSystems;
# const PSID = PowerSimulationsDynamics;

"""
    withName(injector::DynamicInjection, name::String)

creates a copy of `injector` with the given name
"""
function withName(injector::DynamicInjection, name::String)
    newinjector = deepcopy(injector)
    newinjector.name = name
    return newinjector
end



"""
    makeSystems(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)

makes all configurations of the given system and injectors, or specific given configurations.

## Args:
 - `sys::System` : the base system to work off of
 - `injectors` (array of `DynamicInjection`): injectors to use. if one dimensional, all configurations of the given injectors will be returned. If 2d, each row of `injectors` will represent output configuration.
 - `busgroups` (array of String or array of `Vector{String}`, optional): if just array of strings, represents list of buses to consider. If array of Vector{String}, each vector of bus names will be grouped together and always receive the same injector type. If not passed, just considers all busses with Generators attached.

## Returns:
 - `Dict{Vector{String}, System}`: dictionary of {generator names (ordered) => system}. contains all variations.
 """
function makeSystems(sys::System, injectors::Union{AbstractArray{T}, AbstractArray{T, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing) where T <: DynamicInjection
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
            # if this generator's bus is in this group,
            #       append this generator's name to this generator group
            if get_bus(g).name in group
                push!(gengroups[idx], g.name)
                # gengroups[idx] = [g.name, gengroups[idx]...]
                # break
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
```julia
mutable struct GridSearchSys
    base::System
    header::Vector{AbstractString}
    sysdict::Dict{Vector{Any}, Function}
    results_header::Vector{String}
    results_getters::Vector{Function}
    df::DataFrame
    chunksize::Union{Int, Float64}
    hfile::String
end
```

Struct to hold all the systems to iterate through.
Create a GridSearchSys with the constructor, whose signature matches the [`ExperimentUtils.makeSystems`](@ref) method. This will start you off with combinations of injectors. Then you can add sweeps with [`add_generic_sweep!`](@ref)
Use [`execute_sims!`](@ref) to run simulations for all the systems.

`length` and `size` are implemented to give the number of systems and (number of systems, number of columns in header).

## Attributes
 - `base` is the base system. It may be modified, but you should always be able to make any of the actual systems by simply adding to the base.  
 - `header` is a list of the column titles for the dataframe. It should be human-readable.  
 - `sysdict` stores the systems in the format config::Vector => System.  
 - `results_header` is a list of column titles for results.  
 - `results_getters` is a list of functions taking in a GridSearchSystem, Simulation, SmallSignalOutput, and String (error message) which returns results corresponding to `results_header`.  
 - `df` is a DataFrame which will hold the results.  
 - `chunksize` is the number of rows to hold in memory at a time (before saving to file).  
 - `hfile` is a string of function definitions. This allows us to load the serialized data even if user-defined functions aren't present in a new environment.
"""
mutable struct GridSearchSys
    base::System
    header::Vector{AbstractString}
    sysdict::Dict{Vector{Any}, Function}
    results_header::Vector{String}
    results_getters::Vector{Function}
    df::DataFrame
    chunksize::Union{Int, Float64}
    hfile::String
end

"""
    Base.show(io::IO, gss::GridSearchSys)

pretty printing for `GridSearchSys` objects.
"""
function Base.show(io::IO, gss::GridSearchSys)
    print(io, """
    GridSearchSys with $(length(gss)) systems
      base: System (Busses: $(PSID.get_n_buses(gss.base)))
      header: $(string(gss.header))
      sysdict: Dict{Vector{Any}, Function} with $(length(gss.sysdict)) entries
      results_header: $(string(gss.results_header))
      results_getters: $(string(gss.results_getters))
      df: $(nrow(gss.df))x$(ncol(gss.df)) DataFrame
      chunksize: $(gss.chunksize)
      hfile: $(count(';', gss.hfile)) function declarations""")
end

function _add_column!(gss::GridSearchSys, title) push!(gss.header, title) end

"""
    add_result!(gss::GridSearchSys, title::AbstractString, getter::Function)

Add a column to the output dataframe with a result to store. If `gss.df` isn't empty, computes the result and adds it as a column.

`getter` must have the following signature: `(GridSearchSys, Simulation, SmallSignalOutput, String)->(Any)`

Any of the inputs might be `missing`.
"""
function add_result!(gss::GridSearchSys, title::AbstractString, getter::Function)
    gss.hfile *= "function $(string(Symbol(getter))) end; "
    push!(gss.results_header, title)
    push!(gss.results_getters, getter)
    if !isempty(gss.df)
        gss.df[!, title] .= map(x->getter(gss, x...), eachrow(select(gss.df, "sim", "sm", "error", "dt")))
    end
end

"""
    add_result!(gss::GridSearchSys, titles::Vector{T}, getter::Function) where T <: AbstractString

add multiple columns to the output dataframe with results to store. If `gss.df` isn't empty, computes the result and adds it as a column.

`getter` must have the following signature: `(GridSearchSys, Simulation, SmallSignalOutput, String)->(Vector{Any})`

the output vector must be the same length as `titles` (the column titles), and any of the inputs might be `missing`.
"""
function add_result!(gss::GridSearchSys, titles::Vector{T}, getter::Function) where T <: AbstractString
    gss.hfile *= "function $(string(Symbol(getter))) end; "
    push!(gss.results_header, titles...)
    push!(gss.results_getters, getter)
    if !isempty(gss.df)
        data = map(x->getter(gss, x...), eachrow(select(gss.df, "sim", "sm", "error", "dt")))
        # println(data)
        for (idx, title) in enumerate(titles)
            gss.df[!, title] .= map(x->x[idx], data)
        end
    end
end

"""
    set_chunksize!(gss::GridSearchSys, chunksize::Union{Int, Float64})

set the number of rows save in each file (and thus how many to hold in memory before saving to file).

Set to `Inf` to hold all rows in memory (useful for small datasets and to allow use of the dataframe immediately after running the sims)
"""
function set_chunksize!(gss::GridSearchSys, chunksize::Union{Int, Float64})
    @assert chunksize>=1.0 "invalid chunksize: can't save chunks of zero or negative size"
    gss.chunksize = chunksize
end

"""
    GridSearchSys(
        sys::System, 
        injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, 
        busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing} = nothing
    )

constructor for GridSearchSys with the exact same behavior as [`makeSystems`](@ref).
automatically adds the columns `error`, `sim`, `sm`, and `dt` to allow results getter methods to be used after saving.

By default, injector configuration will be represented as a "injector at {bus name or bus names joined by ', '}" for each bus or busgroup.

## Args:
 - `sys::System` : the base system to work off of
 - `injectors` (array of `DynamicInjection`): injectors to use. if one dimensional, all configurations of the given injectors will be returned. If 2d, each row of `injectors` will represent output configuration.
 - `busgroups` (array of String or array of `Vector{String}`, optional): if just array of strings, represents list of buses to consider. If array of Vector{String}, each vector of bus names will be grouped together and always receive the same injector type. If not passed, just considers all busses with Generators attached.

## Returns:
 - `GridSearchSys`: properly initialized GridSearchSys with all the right injectors and the default columns `error`, `sim`, `sm`, and `dt`.
"""
function GridSearchSys(sys::System, injectors::Union{AbstractArray{T}, AbstractArray{T, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing) where T <: DynamicInjection
    if !(sys.data.time_series_storage isa InfrastructureSystems.InMemoryTimeSeriesStorage)
        @warn "Static time series storage detected. Saving to file may miss system time series storage."
    end
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
    gss = GridSearchSys(sys, header, newsysdict, Vector(), Vector(), DataFrame(), Inf, "")
    add_result!(gss, "error", get_error)
    add_result!(gss, "sim", get_sim)
    add_result!(gss, "sm", get_sm)
    add_result!(gss, "dt", get_dt)
    return gss
end


"""
    add_generic_sweep!(gss::GridSearchSys, title::String, adder::Function, params::Vector)

add arbitrary sweeps to a GridSearch. uses `adder` to set a parameter of a system to each of the values in `params`.
 - `title` is the name of the variable this sweep changes.
 - `adder`: `Fn(System, T) -> System`. Modifying the system in-place is OK.
 - `params`: `Vector<T>`
"""
function add_generic_sweep!(gss::GridSearchSys, title::String, adder::Function, params::Vector{T}) where T
    gss.hfile *= "function $(string(Symbol(adder))) end; "
    newsysdict = Dict()
    for (key, s) in gss.sysdict
        for p in params
            # this is called "currying" or something, idk im not a haskell programmer
            # im hungry, curry sounds good rn
            # https://xkcd.com/2210/
            
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

            # we're using higher order functions to insert saved information

            newsysdict[[key..., p]] = ((s, p, adder)->(()->adder(s(), p)))(s, p, adder)
        end
    end
    _add_column!(gss, title)
    gss.sysdict = newsysdict;
end

"""
    add_zipe_sweep!(gss::GridSearchSys, standardLoadFunction::Union{Function, Missing}, zipe_params::Vector{LoadParams})

!!! danger "Deprecation Warning"
    This function is deprecated; just use [`add_generic_sweep!`](@ref)

Adds a ZIPE load sweep to a GridSearchSys. Pass in a standard load to base the ZIPE load off of, and a vector of LoadParams structs to test.

The standard load should be a function which takes in a system and returns an appropriate load. if `missing`, won't add any standard loads.
"""
function add_zipe_sweep!(gss::GridSearchSys, standardLoadFunction::Union{Function, Missing}, zipe_params::Vector{LoadParams})
    sysdict = deepcopy(gss.sysdict)
    if !(standardLoadFunction isa Missing)
        gss.hfile *= "function $(string(Symbol(standardLoadFunction))) end; "
        for s in values(sysdict)
            # this is called "currying" or something, idk im not a haskell programmer
            # im hungry, curry sounds good rn
            # https://xkcd.com/2210/
            
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

            # we're using higher order functions to insert saved information

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
    _add_column!(gss, "ZIPE Load Params")
end

"""
    add_lines_sweep!(gss::GridSearchSys, lineParams::Vector{LineModelParams}, linemodelAdders::Dict{String, Function}=Dict("statpi"=>create_statpi_system, "dynpi"=>create_dynpi_system, "mssb"=>create_MSSB_system))

adds a sweep over line models and parameters.
`lineParams` should be a vector of LineModelParams structs, and `linemodelAdders` should be a dictionary of {model name (human readable) => func!(system, LineModelParams)} with all the model types to add.

This is based around the TLModels.jl package, so linemodelAdders could for example have the pair `"statpi" => create_statpi_system`.
"""
function add_lines_sweep!(gss::GridSearchSys, lineParams::Vector{LineModelParams}, linemodelAdders::Dict{String, Function}=Dict("statpi"=>create_statpi_system, "dynpi"=>create_dynpi_system, "mssb"=>create_MSSB_system))
    for f in values(linemodelAdders)
        gss.hfile *= "function $(string(Symbol(f))) end; "
    end
    _add_column!(gss, "Line Model")
    _add_column!(gss, "Line Params")
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
    length(gss::GridSearchSys)

number of total systems in the gridsearch
"""
function length(gss::GridSearchSys)
    return length(gss.sysdict)
end

"""
    size(gss::GridSearchSys)

returns tuple of (number of systems in gridsearch, number of columns in header)
"""
function size(gss::GridSearchSys)
    return (length(gss.sysdict), length(gss.header))
end


"""
```julia
execute_sims!(
    gss::GridSearchSys, 
    change::PSID.Perturbation; 
    
    tspan::Tuple{Float64, Float64}=(0.48, 0.55), 
    tstops::Vector{Float64}=[0.5]
    dtmax=0.0005, 
    run_transient::Bool=true, 
    log_path::String="sims",
    ida_opts::Dict{Symbol, Any} = Dict(:linear_solver=>:Dense, :max_convergence_failures=>5),
)
```

run simulations on all of the systems in the gridsearch and store the results in a DataFrame.
Fully parallelized.

## Args
 - `gss::GridSearchSys` : the systems
 - `change::Perturbation` : perturbation to apply to the system
 - `tspan::Tuple{Float64, Float64}` : time interval (in seconds) to simulate.
 - `tstops::Vector{Float64}` : `tstops` argument for DifferentialEquations.jl solver - places to force a solver step.
 - `dtmax::Float64` : max timestep for solver (make sure λh is in the feasible region for the solver)
 - `run_transient::Bool` : whether or not to run the transient simulations.
 - `log_path` : folder where outputs will be saved (when chunksize is reached). 
 - `ida_opts` : options for IDA integrator. Defaults are usually OK.

Whenever the number of rows in `gss.df` reaches `gss.chunksize`, results will be saved to file then deleted from gss.df in order to limit total memory usage.

If `gss.chunksize` is finite, the final dataframe will be saved. This way all results will be saved, even the last chunk or if `chunksize` was not reached.

To not save anything to file and keep the results in gss.df, make sure to `set_chunksize(gss, Inf)`.
"""
function execute_sims!(
    gss::GridSearchSys, 
    change::PSID.Perturbation; 
    
    tspan::Tuple{Float64, Float64}=(0.48, 0.55), 
    tstops::Vector{Float64} = [0.5],
    dtmax=0.0005, 
    run_transient::Bool=true, 
    log_path::String="sims",
    ida_opts::Dict{Symbol, Any} = Dict(:linear_solver=>:Dense, :max_convergence_failures=>5),
)
    start = time()

    if !isdir(log_path)
        mkdir(log_path)
    end
    
    gss.df = DataFrame([i=>[] for i in [gss.header..., gss.results_header...]])

    counter = Threads.Atomic{Int}(0)
    total = length(gss)
    function inner(config::Vector{Any}, sys::System)
        (sim, sm, dt, error) = runSim(sys, change, ResidualModel, tspan, tstops, IDA(; ida_opts...), dtmax, run_transient)
        i = Threads.atomic_add!(counter, 1) + 1
        results = vcat(config, reduce(vcat, (getter(gss, sim, sm, error, dt) for getter in gss.results_getters)))
        println("finished solve $i/$total in $(round(Int(dt)/1e9, digits=2))s ($(round(100.0*i/total))%) (runtime: $(round(time()-start))s)")
        return results
    end

    lk = Threads.ReentrantLock()
    chunk_counter = 0
    Threads.@threads for (key, val) in collect(gss.sysdict)
        # key = configuration
        # val = System builder function
        
        # actually run the solve
        results = inner(deepcopy(key), val())

        lock(lk) do
            push!(gss.df, results)
            if nrow(gss.df) >= gss.chunksize
                _save_serde_data(gss.df, log_path*"/results$(chunk_counter).jls")
                deleteat!(gss.df, 1:nrow(gss.df))
                chunk_counter += 1
            end
        end
    end

    # save final chunk if necessary
    if isfinite(gss.chunksize) && nrow(gss.df)>0
        _save_serde_data(gss.df, log_path*"/results$(chunk_counter).jls")
        gss.df = last(gss.df, 0)
        Serialization.serialize(joinpath(log_path, ".gss"), gss)
        open(joinpath(log_path, ".hfile"), "w") do file
            write(file, gss.hfile)
        end
    end
end

function _save_serde_data(df::DataFrame, path::String)
    println("dirname: ", dirname(path))
    println("path: ", path)
    if !isdir(dirname(path))
        mkdir(dirname(path))
    end
    Serialization.serialize(path, df)
end


"""
    save_serde_data(gss::GridSearchSys, path::String)

serializes the GridSearchSys and saves it to `path` to be read later using `Serialization.deserialize` (through `load_serde_data`)

deletes all .jls files and the .gss file in `path` if `path` is a directory containing those files.
"""
function save_serde_data(gss::GridSearchSys, path::String)
    if !(gss.base.data.time_series_storage isa InfrastructureSystems.InMemoryTimeSeriesStorage)
        @warn "Static time series storage detected. System objects may not retain time series data."
    end
    if !isdir(path)
        mkdir(path)
    else
        for file in filter(x->occursin(r"\.jls", x) || (x ∈ [".gss", ".hfile"]), readdir(path))
            rm(joinpath(path, file))
        end
    end
    counter = Threads.Atomic{Int}(0)
    numfiles = Int(ceil(nrow(gss.df)/gss.chunksize))+1
    progress_bar_width() = displaysize(stdout)[2]-28-Int(2*(floor(log10(numfiles))+1))-length(path)
    print("Writing files to $path: |"*(" "^progress_bar_width())*"| (0/$(numfiles))")
    Threads.@threads for (i, df) in collect(enumerate(Iterators.partition(gss.df, isfinite(gss.chunksize) ? gss.chunksize : nrow(gss.df))))
        Serialization.serialize(joinpath(path, "results$(i).jls"), df)
        files_written = Threads.atomic_add!(counter, 1) + 1
        boxes = Int(round((files_written/numfiles)*progress_bar_width()))
        print("\r"*(" "^(displaysize(stdout)[2])))
        print("\rWriting files to $path: |"*("@"^boxes)*(" "^(progress_bar_width()-boxes))*"| ($files_written/$numfiles)")
    end

    Serialization.serialize(
        joinpath(path, ".gss"),
        GridSearchSys(
            gss.base,
            gss.header,
            gss.sysdict,
            gss.results_header,
            gss.results_getters,
            last(gss.df, 0),
            gss.chunksize,
            gss.hfile,
        )
    )
    print("\rWriting files to $path: |"*("@"^progress_bar_width())*"| ($numfiles/$numfiles)\n")

    open(joinpath(path, ".hfile"), "w") do file
        write(file, gss.hfile);
    end
    println("Done!")
end

"""
    load_serde_data(path::String)

Loads serialized dataframe from file or folder.

If `path` is a folder, looks for all non-hidden .jls files, reads them, and concatenates them.
"""
function load_serde_data(path::String)
    if !isdir(path)
        if !isfile(path) AssertionError("Could not load data from path: `$path` does not exist.") end
        return Serialization.deserialize(path)
    end
    files = collect(
    map(   file -> path*"/"*file, 
    filter(file -> (file[1]!='.')&&(file[end-3:end]==".jls"), 
        readdir(path)
    )))
    # files = [joinpath(path, file) for file in readdir(path) if (file[1]!='.')&&(file[end-3:end]==".jls")]
    if ".hfile" ∈ readdir(path)
        # println(abspath(path))
        # println(joinpath(path, ".hfile"))
        include(joinpath(abspath(path), ".hfile"))
    end
    # println(files)
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
    println("\nDone!")
    if ".gss" ∈ readdir(path)
        gss = load_serde_data(joinpath(path, ".gss"))
        Meta.parse(gss.hfile)
        gss.df = vcat(dfs...)
        return gss
    else
        return vcat(dfs...)
    end
end


# """
# extract array of current magnitude over time at the object with name `name` from the simulation results `res`.

# Timestamps are thrown away. If you're using the GridSearchSystem, the returned array will correspond to an even time grid with spacing `dtmax` (the argument you passed to [`execute_sims`](@ref), default 0.02s).
# """
# function current_magnitude(res::SimulationResults, name::String)
#     _, ir = PSID.post_proc_real_current_series(res, name, nothing)
#     _, ii = PSID.post_proc_imaginary_current_series(res, name, nothing)
#     # _, ir = get_state_series(res, (name, :ir_cnv))
#     # _, ii = get_state_series(res, (name, :ii_cnv))
#     return sqrt.(ir.^2 .+ ii.^2)
# end


# """
# extracts generator speed time series from `res`. Time discretization is even `dtmax` (the argument you passed to [`execute_sims`](@ref), default 0.02s)
# """
# function generator_speed(res::SimulationResults, name::String)
#     return get_state_series(res, (name, :ω))[2]
# end

# """
# expand line and load param columns into each individual parameter. This improves usability of the resulting dataframe, and makes saving to TSV cleaner.

# Simply adds a column for each attribute of the line parameter struct and/or the load parameter struct.
# """
# function expand_columns!(gss::GridSearchSys, remove_expanded=false)
#     expand_columns!(gss.df, remove_expanded)
# end
# function expand_columns!(df::DataFrame, remove_expanded=false)
#     columns_to_include = []
#     if "Line Params" in names(df)
#         push!(columns_to_include, :"Line Params")
#         for i in fieldnames(LineModelParams)
#             df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"Line Params")
#         end
#     end

#     if "ZIPE Load Params" in names(df)
#         push!(columns_to_include, :"ZIPE Load Params")
#         for i in fieldnames(LoadParams)
#             df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"ZIPE Load Params")
#         end
#     end
#     if remove_expanded
#         select!(df, Not(columns_to_include))
#     end
# end

"""
    runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(linear_solver=:LapackDense, max_convergence_failures=5), dtmax=0.001, run_transient=true, log_path::String=mktempdir())

little wrapper to run simulations

## Args:
 - `system`: the system to simulate
 - `change`: perturbation to apply for the transient sim
 - `model`: model to pass into Simulation()
 - `tspan`: time span for transient simulation
 - `tstops`: places to force solver timesteps
 - `solver`: DE solver for transient simulation
 - `dtmax`: maximum timestep for DE solver
 - `run_transient`: whether or not to actually perform the transient simulation
 - `log_path`: path for simulation logs.

## Returns:
 - `(Simulation, SmallSignalOutput, UInt64, String)`: the Simulation object, the small signal analysis object, the time in nanoseconds the simulation took, and the error message. All but the time might be `missing` if things failed.
"""
function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), tstops=[0.5], solver=IDA(linear_solver=:LapackDense, max_convergence_failures=5), dtmax=0.001, run_transient=true, log_path::String=mktempdir())
    tic = Base.time_ns()
    local sim, sm
    solve_ac_powerflow!(system)
    # println("HERE1")
    sim = Simulation(
        model,
        system,
        log_path,
        tspan,
        change,
        disable_timer_outputs=true, # needed for multiprocessing
        # console_level=Logging.Error,
        # file_level=Logging.Error
        # initialize_simulation=false,
    )
    # println("HERE2")
    try
        sm = small_signal_analysis(sim)
        # println(sm)
    catch err
        return (sim, missing, Base.time_ns()-tic, "Small Signal Analysis failed with error $err")
    end
    if run_transient
        try
            execute!(
                sim,
                solver,
                dtmax = dtmax,
                # saveat = output_res,
                tstops = tstops,
                enable_progress_bar = false, # with multithreading it's meaningless anyways
            )
        catch err
            return (sim, sm, Base.time_ns()-tic, "Transient sim failed with error $err")
        end
    end

    return (sim, sm, Base.time_ns()-tic, missing)
end

"""
    get_permutations(iterable::Iterable, k::Int)

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