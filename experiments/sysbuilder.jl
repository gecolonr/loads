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

function makeSystems(sys::System, injectors::AbstractArray{DynamicInjection}, busgroups::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing}=nothing)
    if isnothing(busgroups)
        busgroups = [[get_bus(i).name] for i in get_components(Generator, sys)]
    elseif busgroups[1] isa String
        busgroups = [[i] for i in busgroups]
    end
    println(busgroups)
    gengroups = [[] for _ in busgroups]
    for g in get_components(Generator, sys)
        for (idx, group) in enumerate(busgroups)
            if get_bus(g).name in group
                gengroups[idx] = [g.name, gengroups[idx]...]
                break
            end
        end
    end
    combos = with_replacement_combinations(injectors, length(busgroups))
    println(combos)
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

function runSim(system, change=BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), model=ResidualModel, tspan=(0., 5.), solver=IDA(), dtmax=0.02)
    sim = Simulation(
        model,
        system,
        mktempdir(),
        tspan,
        change,
    )
    execute!(
        sim,
        solver,
        dtmax = dtmax,
    )
    # results = read_results(sim)
    sm = small_signal_analysis(sim)
    # plot(sm.eigenvalues, seriestype=:scatter, label=L"eigenvalues $\lambda$")
    return (sim, sm)
end