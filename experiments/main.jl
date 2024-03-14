using Distributed
addprocs(14)
@everywhere cd(@__DIR__)
@everywhere cd("..")
@everywhere using Pkg
@everywhere Pkg.activate(".")

@everywhere include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

sys() = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"))


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

machine_oneDoneQ() = OneDOneQMachine(
    0.1, #R
    1.3125, #Xd
    1.2578, #Xq
    0.1813, #Xd_p
    0.25, #Xq_p
    5.89, #Td0_p
    0.6, #Tq0_p
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
sysdict = makeSystems(sys(), [case_gen(), case_inv()])
combos, syses = keys(sysdict), values(sysdict)
results = pmap(runSim, syses)

results = Dict(zip(combos, results))
Dict(zip(keys(results), [i[2] for i in values(results)]))
