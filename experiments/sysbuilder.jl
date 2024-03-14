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
    )
    sm = small_signal_analysis(sim)
    if run_transient
        execute!(
            sim,
            solver,
            dtmax = dtmax,
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