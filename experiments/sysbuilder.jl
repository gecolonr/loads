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
using ArgCheck
using DataFrames
import Logging

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
 - `busgroups` (array of String or array of `Vector{String}`, optional): if just array of strings, represents list of buses to consider. If array of Vector{String}, each vector of bus names will be grouped together and always receive the same injector type. 

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
"""
mutable struct GridSearchSys
    base::System
    header::Vector{String}
    sysdict::Dict{Vector{Any}, System}
    results_header::Vector{String}
    results_getters::Vector{Function}
end

"""
Adds a column to the header. literally just `push!(gss.header, title)`. for cleanliness and encapsulation I guess.
"""
function add_column!(gss::GridSearchSys, title)
    push!(gss.header, title)
end

function add_result!(gss::GridSearchSys, title::String, getter::Function)
    push!(gss.results_header, title)
    push!(gss.results_getters, getter)
end


"""
constructor for GridSearchSys with the exact same behavior as [`makeSystems`](@ref).
"""
function GridSearchSys(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)
    sysdict = makeSystems(sys, injectors, busgroups)
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
    return GridSearchSys(sys, header, sysdict, Vector(), Vector())
end

"""
Adds a ZIPE load sweep to a GridSearchSys. Pass in a standard load to base the ZIPE load off of, and a vector of LoadParams structs to test.

The standard load should be a function which takes in a system and returns an appropriate load.
"""
function add_zipe_sweep!(gss::GridSearchSys, standardLoadFunction::Function, zipe_params::Vector{LoadParams})
    sysdict = deepcopy(gss.sysdict)
    for s in values(sysdict)
        add_component!(s, standardLoadFunction(s))
    end
    gss.sysdict = Dict()
    for params in zipe_params
        for (key, val) in sysdict
            withzipe = deepcopy(val)
            create_ZIPE_load(withzipe, params)
            gss.sysdict[[key..., params]] = withzipe
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
                newsysdict[[key..., linename, params]] = creator(val, params)
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

fully parallelized. Be aware that the memory usage is quite high.

## Args
 - gss::GridSearchSys : the systems
 - change : perturbation to apply to the system
 - tspan::Tuple{Float64, Float64} : time interval (in seconds) to simulate.
 - dtmax::Float64 : max timestep for solver, and resolution of saved timeseries results
 - run_transient::Bool : whether or not to run the transient simulations.

## Returns
dataframe of results, with the following columns:
| column name                 | content                                                                                                                                |
| --------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| gss.header                  | Everything in the GridSearchSys's header (all each individual columns). This is the system configuration.                              |
| eigenvalues                 | λ's from the small signal analysis.                                                                                                    |
| eigenvectors                | Eigenvectors from the small signal analysis.                                                                                           |
| Voltage at Bus \$n          | Voltage time series, with timesteps `dtmax`, at each bus. One column per bus.                                                          |
| Inverter Current at Bus \$n | Inverter current magnitude at bus \$n, or `missing` if there is no inverter there. One column per bus.                                 |
| Generator Speed at Bus \$n  | Angular velocity of generator at bus \$n, or `missing` if there is no generator there. One column per bus.                             |
| Simulation Status           | String representation of either the simulation status or the error caught during the call to [`runSim`](@ref).                         |
| Simulation Time (ns)        | Time in nanoseconds to build the simulation, perform the small signal analysis, and (if `run_transient`) run the transient simulation. | 

the columns `eigenvalues` and `eigenvectors` may be `missing` if the small signal analysis failed to converge.
The columns with transient simulation results and the Simulation Time column may be `missing` if the small signal analysis or the transient simulation failed.
"""
function executeSims(gss::GridSearchSys, change, tspan::Tuple{Float64, Float64}=(0.0, 3.0), dtmax=0.02, run_transient::Bool=true)
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    bus_voltages = (x->"Voltage at Bus $(get_number(x))").(get_components(Bus, gss.base))
    inverter_currents = (x->"Inverter Current at Bus $(get_number(x))").(gen_busses)
    generator_speeds = (x->"Generator Speed at Bus $(get_number(x))").(gen_busses)
    cols = [gss.header..., "eigenvalues", "eigenvectors", bus_voltages..., inverter_currents..., generator_speeds..., "Simulation Status", "Simulation Time (ns)"]
    transient_results_length = length(bus_voltages) + length(inverter_currents) + length(generator_speeds)
    df = DataFrame([i=>[] for i in cols])
    counter = Threads.Atomic{Int}(0)
    total = length(gss)
    function inner(config::Vector{Any}, sys::System)
        local sim, sm, time
        try
            (sim, sm, time) = runSim(sys, change, ResidualModel, tspan, IDA(linear_solver=:LapackBand, max_convergence_failures=10), dtmax, run_transient)
        catch error
            i = Threads.atomic_add!(counter, 1) + 1
            println("finished solve $i/$total in ----s ($(round(100.0*i/total))%)")
            return (config..., missing, missing, Array{Missing}(missing, transient_results_length)..., string(error), missing)
        else
            i = Threads.atomic_add!(counter, 1) + 1
            println("finished solve $i/$total in $(round(Int(time)/1e9, digits=2))s ($(round(100.0*i/total))%)")
            if string(sim.status) != "SIMULATION_FINALIZED"
                return (config..., sm.eigenvalues, sm.eigenvectors, Array{Missing}(missing, transient_results_length)..., string(sim.status), time)
            else
                res = read_results(sim)
                inverters = [i for i in get_components(DynamicInverter, sys) if i.name in keys(gen_dict)]
                inverters = Dict(map(x->gen_dict[x], get_name.(inverters)) .=> inverters)
                generators = [i for i in get_components(DynamicGenerator, sys) if i.name in keys(gen_dict)]
                generators = Dict(map(x->gen_dict[x], get_name.(generators)) .=> generators)
                
                return (
                    config..., 
                    sm.eigenvalues, 
                    sm.eigenvectors, 
                    [get_voltage_magnitude_series(res, i)[2] for i in get_number.(get_components(Bus, gss.base))]...,
                    [(i in keys(inverters) ? current_magnitude(res, inverters[i].name) : missing) for i in get_number.(gen_busses)]...,
                    [(i in keys(generators) ? generator_speed(res, generators[i].name) : missing) for i in get_number.(gen_busses)]...,
                    string(sim.status),
                    time
                    )
            end
        end
    end
        
    lk = ReentrantLock()
    Threads.@threads for (key, val) in collect(gss.sysdict)
        results = inner(key, val)
        lock(lk) do 
            push!(df, results)
        end
    end
    return df
end

function get_eigenvalues(_gss::GridSearchSys, _sim::Simulation, sm::PSID.SmallSignalOutput)
    return sm isa Missing ? missing : [sm.eigenvalues]
end

function get_eigenvectors(_gss::GridSearchSys, _sim::Simulation, sm::PSID.SmallSignalOutput)
    return sm isa Missing ? missing : [sm.eigenvectors]
end

function get_bus_voltages(gss::GridSearchSys, sim::Simulation, _sm::PSID.SmallSignalOutput)
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Bus, gss.base)))
    else
        return [get_voltage_magnitude_series(sim.results, i)[2] for i in get_number.(get_components(Bus, gss.base))]
    end
end

function get_inverter_currents(gss::GridSearchSys, sim::Simulation, _sm::PSID.SmallSignalOutput)
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Generator, gss.base)))
    end
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    inverters = [i for i in get_components(DynamicInverter, sys) if i.name in keys(gen_dict)]
    inverters = Dict(map(x->gen_dict[x], get_name.(inverters)) .=> inverters)
    return [(i in keys(inverters) ? current_magnitude(res, inverters[i].name) : missing) for i in get_number.(gen_busses)]
end

function get_generator_speeds(gss::GridSearchSys, sim::Simulation, _sm::PSID.SmallSignalOutput)
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Generator, gss.base)))
    end
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    generators = [i for i in get_components(DynamicGenerator, sys) if i.name in keys(gen_dict)]
    generators = Dict(map(x->gen_dict[x], get_name.(generators)) .=> generators)
    return [(i in keys(generators) ? generator_speed(res, generators[i].name) : missing) for i in get_number.(gen_busses)]
end

function get_sim_status(_gss::GridSearchSys, sim::Simulation, _sm::PSID.SmallSignalOutput)
    return (sim isa Missing) ? missing : sim.status
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
function expand!(df)
    if "Line Params" in names(df)
        for i in fieldnames(LineModelParams)
            df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"Line Params")
        end
    end

    if "ZIPE Load Params" in names(df)
        for i in fieldnames(LoadParams)
            df[!, i] = (x->x isa Missing ? missing : getfield(x, i)).(df.var"ZIPE Load Params")
        end
    end

    select!(df, Not([:"Line Params", :"ZIPE Load Params"]))
end

"""little wrapper to run simulations

## Args:
 - `system`: the system to simulate
 - `change`: perturbation to apply for the transient sim
 - `model`: model to pass into Simulation()
 - `tspan`: time span for transient simulation
 - `solver`: DE solver for transient simulation
 - `dtmax`: maximum timestep for DE solver
 - `run_transient`: whether or not to actually perform the transient simulation

## Returns:
 - `(Simulation, SmallSignalOutput, UInt64)`: the Simulation object, the small signal analysis object, and the time in nanoseconds the simulation took
"""
function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(linear_solver=:LapackDense, max_convergence_failures=5), dtmax=0.02, run_transient=true)
    tic = Base.time_ns()
    sim = Simulation(
        model,
        system,
        mktempdir(),
        tspan,
        change,
        file_level = Logging.Warn,
    )
    sm = small_signal_analysis(sim)
    if run_transient
        execute!(
            sim,
            solver,
            dtmax = dtmax,
            saveat = dtmax,
            enable_progress_bar = false, # with multithreading it's meaningless anyways
        )
    end
    toc = Base.time_ns()
    # println("FINISHED A SOLVE")
    return (sim, sm, (toc-tic)::UInt64)
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

