cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
# using Wandb, Logging
using InfrastructureSystems
using Plots
using PlotlyJS, DataFrames

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

# system builder
sys() = System(joinpath(pwd(), "data/raw_data/OMIB.raw"))
# instantiate system
s = sys()

# add a second generator (initially it's just one generator and an infinite bus)
# taken from some json in the ZIPE_loads repo
add_component!(s, ThermalStandard(
    "generator2", 
    true, 
    true, 
    first(get_components(ACBus, s)), 
    0.5, 
    0.0, 
    1.4142135623730951, 
    (min = 0.0, max = 1.0), 
    (min = -1.0, max = 1.0), 
    (up = 1.0, down = 1.0), 
    ThreePartCost(VariableCost((0.0, 1.0)), 0.0, 0.0, 0.0), 
    100.0, 
    nothing, 
    false, 
    19, 
    14, 
    Service[],
    10000.0,
    nothing,
    Dict("z_source" => (r = 0.0, x = 1.0)))
)

# Adding a second line to the OMIB system so we can trip one without it all breaking
line_cp = deepcopy(first(get_components(Line, s)))
line_cp.name = "otherline" # change name so they're not identical
# add new time series container so they're not both trying to use the same one
line_cp.time_series_container = InfrastructureSystems.TimeSeriesContainer();
# add the line to the system!
add_component!(s, line_cp)

# functions to create our machines: one inverter and one generator
case_inv() = DynamicInverter(
    "I", # stands for "Inverter"
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)

case_gen() = DynamicGenerator(
    "G", # stands for "Generator"
    1.0, # ω_ref,
    AF_machine(), #machine
    shaft_no_damping(), #shaft
    avr_type1(), #avr
    tg_none(), #tg
    pss_none(), #pss
)


# get all combinations of generators on this system
sysdict = makeSystems(s, [case_inv(), case_gen()]);

# function to create a standard load for a particular system
# this is what the ZIPE load will attach itself to (?)
load(s) = StandardLoad(
    name="load",
    available=true,
    bus=first(get_components(ACBus, s)),
    base_power=100.0,
    constant_active_power=2.0,
    constant_reactive_power=0.1,
)

# add standard load to all systems
for s in values(sysdict)
    add_component!(s, load(s))
end

function gridsearch(dx=0.1, Zmax=1.0, Imax=1.0, Pmax=1.0)
    """returns generator of all Z, I, P, and E combinations"""
    return ([i j k 1.0-i-j-k] for i in 0.0:dx:Zmax for j in 0.0:dx:min(1.0-i, Imax) for k in 0.0:dx:min(1.0-i-j, Pmax))
end

# TODO: parallelize this
function zipe_gridsearch(systems)
    """performs a sweep over all ZIPE parameters and all systems in `systems`.

    Args:
        systems (Dict{Vector{String}=>System}): dictionary of all the system variants to test
    
    Returns:
        DataFrame: df of results with columns [combo (`systems` key type), Z, I, P, E, number of positive eigenvalues, max eigenvalue, simulation status]
    """
    # initialize df to hold all data
    df = DataFrame(combo=Vector{String}[], Z=Float64[], I=Float64[], P=Float64[], E=Float64[], n_pos_eigs=Int[], max_eig=Float64[], sim_status=String[]);

    # loop through all combinations of ZIPE parameters
    for params in gridsearch()
        # loop through all system configurations
        for (combo, s) in systems
            # make a new system for this simulation
            sys = deepcopy(s)
            # add ZIPE load to the system with `params`
            create_ZIPE_load(sys, LoadParams(params...))

            # try/catch block to catch simulation convergence errors
            try
                # runSim is from sysbuilder.jl
                (sim, sm) = runSim(
                    sys, # system to simulate
                    BranchTrip(0.5, ACBranch, first(get_components(ACBranch, sys)).name), # perturbation
                    ResidualModel, # model to use (we're using IDA so this is the right one)
                    (0.0, 5.0), # time span
                    IDA(), # DE solver
                    0.02, # max dt
                    false, # true=run transient simulation, false=don't
                )
                # add results as a new row of the df
                push!(df, (combo, params[1], params[2], params[3], params[4], count(x->(x.re>0.0), sm.eigenvalues), findmax(map(x->x.re, sm.eigenvalues))[1], string(sim.status)))
                # status print (removed for cleanliness)
                # println(combo, " had ", count(x->(x.re>=0.0), sm.eigenvalues), " positive eigenvalue(s), with max real part ", findmax(map(x->x.re, sm.eigenvalues))[1])
            catch e
                # this catches cases where no equilibrium point could be found
                # we exclude these cases from the output data because they are infeasible.
                continue
            end
            println(params) # so users know approximately how far along we are in the sim
        end
    end
    return df
end

# run gridsearch
df = zipe_gridsearch(sysdict)

# create trace for plotlyjs
mytrace = parcoords(
    ;line = attr(color=df.n_pos_eigs),
    dimensions = [
        # system configuration
        attr(range = [1, length(combos)], label = "generators", values = map(x->Dict(zip(keys(sysdict), 1:(length(sysdict))))[x], df.combo),tickvals = 1:length(sysdict),ticktext=[i[1]*"-"*i[2] for i in keys(sysdict)]),
        # ZIPE parameters
        attr(range = [findmin(df.Z)[1], findmax(df.Z)[1]], label = "Z", values = df.Z),
        attr(range = [findmin(df.I)[1], findmax(df.I)[1]], label = "I", values = df.I),
        attr(range = [findmin(df.P)[1], findmax(df.P)[1]], label = "P", values = df.P),
        attr(range = [findmin(df.E)[1], findmax(df.E)[1]], label = "E", values = df.E),
        # number of positive eigenvalues
        attr(range = [0,findmax(df.n_pos_eigs)[1]], label = "n pos eigs" , values = df.n_pos_eigs),
        # maximum real component
        attr(range = [log(findmin(df.max_eig)[1]), log(findmax(df.max_eig)[1])], label="ln(max eig)", values=log.(df.max_eig))
    ]);

# define layout
layout = Layout(
    title_text="ZIPE Load Sweep",
    title_x=0.5,
    title_y=0.99,
)

# make the plot!
myplot = PlotlyJS.plot(mytrace,layout)