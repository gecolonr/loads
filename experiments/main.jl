cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")

using PowerSystems
using PowerSimulationsDynamics
using ZIPE_loads
using TLmodels

const PSY = PowerSystems
const PSID = PowerSimulationsDynamics

include("../data/build_scripts/device_models.jl")

sys = System(joinpath(pwd(), "../raw_data/WSCC_9bus.raw"))


case_inv = DynamicInverter(
    get_name(g),
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)

case_gen = DynamicGenerator(
    get_name(g),
    1.0, # ω_ref,
    AF_machine(), #machine
    shaft_no_damping(), #shaft
    avr_type1(), #avr
    tg_none(), #tg
    pss_none(), #pss
)