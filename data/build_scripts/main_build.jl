cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("device_models.jl")

sys = System(joinpath(pwd(), "../raw_data/OMIB.raw"))

slack_bus = [b for b in get_components(Bus, sys) if get_bustype(b) == BusTypes.REF][1]

inf_source = Source(
           name = "InfBus", #name
           available = true, #availability
           active_power = 0.0,
           reactive_power = 0.0,
           bus = slack_bus, #bus
           R_th = 0.0, #Rth
           X_th = 5e-6, #Xth
       )

add_component!(sys, inf_source)



for g in get_components(Generator, sys)
    case_inv = DynamicInverter(
        "generator-102-1",
        1.0, # ω_ref,
        converter_high_power(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source_lv(), #dc source
        pll(), #pll
        filt(), #filter
    )

#Attach the dynamic inverter to the system
    add_component!(sys, case_inv, g)
end

to_json(sys, joinpath(pwd(), "../json_data/SIIB.json"), force = true)