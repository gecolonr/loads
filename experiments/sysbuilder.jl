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

function withName(injector::DynamicInjection, name::String)
    """creates a copy of `injector` with the given name"""
    newinjector = deepcopy(injector)
    newinjector.name = name
    return newinjector
end

function makeSystems(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)
    """makes all configurations of the given system and injectors
    
    Args:
        sys (System): the base system to work off of
        injectors (array of DynamicInjection): injectors to use. if one dimensional, all configurations of the given injectors will be returned. If 2d, each row of `injectors` will represent output configuration.
        busgroups (array of string or array of Vector{String}, optional): if just array of strings, represents list of buses to consider. If array of Vector{String}, each vector of bus names will be grouped together and always receive the same injector type. 
    
    Returns:
        Dict{Vector{String}=>System}: dictionary of {generator names, ordered => system}. contains all variations.
    """
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
        combos = injectors # if 2d, don't get permutations - user wants specific manually listed configurations
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

mutable struct GridSearchSys
    base::System
    header::Vector{Any}
    sysdict::Dict{Vector{Any}, System}
end

function GridSearchSys(base::System)
    return GridSearchSys(base, Vector(), Dict())
end

function add_column!(gss::GridSearchSys, title)
    push!(gss.header, title)
end

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
    return GridSearchSys(sys, header, sysdict)
end


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
function length(gss::GridSearchSys)
    return length(gss.sysdict)
end
function size(gss::GridSearchSys)
    return (length(gss.sysdict), length(gss.header))
end
function executeSims(gss::GridSearchSys, change, tspan::Tuple{Float64, Float64}=(0.0, 3.0), dtmax=0.02, run_transient::Bool=true)
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    bus_voltages = (x->"Voltage at Bus $(get_number(x))").(get_components(Bus, gss.base))
    inverter_currents = (x->"Inverter Current at Bus $(get_number(x))").(gen_busses)
    generator_speeds = (x->"Generator Speed at Bus $(get_number(x))").(gen_busses)
    cols = [gss.header..., "eigenvalues", "eigenvectors", bus_voltages..., inverter_currents..., generator_speeds..., "Simulation Status"]
    transient_results_length = length(bus_voltages) + length(inverter_currents) + length(generator_speeds)
    df = DataFrame([i=>[] for i in cols])
    lk = ReentrantLock()
    Threads.@threads :static for (config, sys) in collect(gss.sysdict)
        local sim, sm
        try
            (sim, sm) = runSim(sys, change, ResidualModel, tspan, IDA(), dtmax, run_transient)
        catch error
            lock(lk) do
                push!(df, (config..., missing, missing, Array{Missing}(missing, transient_results_length)..., string(error)))
            end
        else
            if string(sim.status) != "SIMULATION_FINALIZED"
                lock(lk) do
                    push!(df, (config..., sm.eigenvalues, sm.eigenvectors, Array{Missing}(missing, transient_results_length)..., string(sim.status)))
                end
            else
                res = read_results(sim)
                inverters = [i for i in get_components(DynamicInverter, sys) if i.name in keys(gen_dict)]
                inverters = Dict(map(x->gen_dict[x], get_name.(inverters)) .=> inverters)
                generators = [i for i in get_components(DynamicGenerator, sys) if i.name in keys(gen_dict)]
                generators = Dict(map(x->gen_dict[x], get_name.(generators)) .=> generators)

                lock(lk) do
                    push!(df, (
                        config..., 
                        sm.eigenvalues, 
                        sm.eigenvectors, 
                        [get_voltage_magnitude_series(res, i)[2] for i in get_number.(get_components(Bus, gss.base))]...,
                        [(i in keys(inverters) ? current_magnitude(res, inverters[i].name) : missing) for i in get_number.(gen_busses)]...,
                        [(i in keys(generators) ? generator_speed(res, generators[i].name) : missing) for i in get_number.(gen_busses)]...,
                        string(sim.status)
                    ))
                end
            end
        end
    end
    return df
end

function current_magnitude(res, name)
    _, ir = get_state_series(res, (name, :ir_cnv))
    _, ii = get_state_series(res, (name, :ii_cnv))
    return sqrt.(ir.^2 .+ ii.^2)
end

function generator_speed(res, name)
    return get_state_series(res, (name, :ω))[2]
end

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

function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(), dtmax=0.02, run_transient=true)
    """little wrapper to run simulations
    
    Args:
        system: the system to simulate
        change: perturbation to apply for the transient sim
        model: model to pass into Simulation()
        tspan: time span for transient simulation
        solver: DE solver for transient simulation
        dtmax: maximum timestep for DE solver
        run_transient: whether or not to actually perform the transient simulation
    
    Returns:
        (sim, sm): the Simulation object (sim) and the small signal analysis object (sm)
    """
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
    return (sim, sm)
end

function get_permutations(iterable, k)
    """gets all permutations of `iterable` of length `k` with replacement.

    Args:
        iterable (Iterable): some iterable (with well-defined length and indexing)
        k (Int): desired permutation length
    
    Returns:
        Base.Generator: Generator of vectors of length `k` of combinations of elements of `iterable`
    """
    # yes this is really the best way I could think of.
    # computes all k-digit numbers in base-n, and then indexes `iterable` with the digits.
    n = length(iterable)
    return ([iterable[parse(Int, j, base=n)+1] for j in string(i, base=n, pad=k)] for i in 0:(n^k-1))
end