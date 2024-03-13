cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
using Wandb, Dates, Logging

lg = WandbLogger(; project = "Wandb.jl", name = nothing)


include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

sys() = System(joinpath(pwd(), "data/raw_data/OMIB.raw"))
s = sys()
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
case_inv() = DynamicInverter(
    "DynamicInverter",
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)

case_gen() = DynamicGenerator(
    "DynamicGenerator",
    1.0, # ω_ref,
    AF_machine(), #machine
    shaft_no_damping(), #shaft
    avr_type1(), #avr
    tg_none(), #tg
    pss_none(), #pss
)

s, combos = makeSystems(s, [case_inv(), case_gen()]);
s = [s[1], s[3], s[5]]
combos = [combos[1], combos[3], combos[5]]
load(s) = StandardLoad(
    name="load",
    available=true,
    bus=first(get_components(ACBus, s)),
    base_power=100.0,
    constant_active_power=2.0,
    constant_reactive_power=0.1,
)

load_params() = LoadParams(
    z_percent = 0.8,
    i_percent = 0.05,
    p_percent = 0.05,
    e_percent = 0.1,
)
for i in s
    add_component!(i, load(i))
    # create_ZIPE_load(i, load_params())
end

function gridsearch(dx=0.1)
    return ([i j k 1-i-j-k] for i in 0.0:dx:1.0 for j in 0.0:dx:(1.0-i) for k in 0.0:dx:(1.0-i-j))
end

function zipe_gridsearch(systems)
    for params in gridsearch()
        for (idx, s) in enumerate(systems)
            sys = deepcopy(s)
            create_ZIPE_load(sys, LoadParams(params...))
            (sim, sm) = runSim(
                s, 
                BranchTrip(0.5, ACBranch, first(get_components(ACBranch, sys)).name),
                ResidualModel,
                (0.0, 5.0),
                IDA(),
                0.02,
                false,
            )
            # Wandb.log(lg, Dict(
            #     "z"=>params[1],
            #     "i"=>params[2],
            #     "p"=>params[3],
            #     "e"=>params[4],
            #     "power in"=>["i-i", "i-g", "g-g"][idx],
            #     "n pos eigs"=>count(x->(x>0), sm.eigenvalues),
            # ))
        end
    end
end

