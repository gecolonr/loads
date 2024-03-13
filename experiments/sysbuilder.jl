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

const PSY = PowerSystems
const PSID = PowerSimulationsDynamics

function withName(injector::DynamicInjection, name::String)
    newinjector = deepcopy(injector)
    newinjector.name = name
    return newinjector
end

function makeSystems(sys::System, injectors::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)
    if isnothing(busgroups)
        busgroups = [[get_bus(i).name] for i in get_components(Generator, sys)]
    elseif busgroups[1] isa String
        busgroups = [[i] for i in busgroups]
    end
    # println(busgroups)
    gengroups = [[] for _ in busgroups]
    for g in get_components(Generator, sys)
        for (idx, group) in enumerate(busgroups)
            if get_bus(g).name in group
                gengroups[idx] = [g.name, gengroups[idx]...]
                break
            end
        end
    end
    if length(size(injectors))>1
        combos = injectors
        print("HERE1")
    else
        combos = with_replacement_combinations(injectors, length(busgroups))
        combos = reduce(vcat, [collect(permutations(i)) for i in combos])
        print("HERE2")
    end
    # println(combos)
    systems = Array{System}(undef, length(combos))
    for (comboIdx, combo) in enumerate(combos)
        systems[comboIdx] = deepcopy(sys)
        injectors = [
            withName(injector, gengroups[idx][jdx]) 
                for (idx, injector) in enumerate(combo) 
                    for jdx in 1:length(gengroups[idx]) 
                ]
        for i in injectors
            add_component!(systems[comboIdx], i, first(get_components_by_name(Generator, systems[comboIdx], i.name)))
        end
        
        
        println("here, %d", length(systems))
    end
    return systems, [[i.name for i in combo] for combo in combos]
end

function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(), dtmax=0.02, run_transient=true)
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
        return (sim, sm)
    else
    # results = read_results(sim)
    # plot(sm.eigenvalues, seriestype=:scatter, label=L"eigenvalues $\lambda$")
        return (0, sm)
    end
end

function get_permutations(iterable, k)
    if k == 1
        return iterable
    end
    return Iterators.flatten(
        (Iterators.flatten((i, j)) for j in get_permutations(iterable, k-1))
            for i in iterable
    )
end