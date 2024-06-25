cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sys = System(joinpath(pwd(), "../raw_data/WSCC_9bus.raw"))

slack_bus =[b for b in get_components(Bus, sys) if get_bustype(b) == BusTypes.REF][1]

include("device_models.jl")

for g in get_components(Generator, sys)
    if get_number(get_bus(g)) == 1
        #Create the dynamic inverter
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
        #Attach the dynamic inverter to the system
        add_component!(sys, case_inv, g)
    #Find the generator at bus 102
    elseif get_number(get_bus(g)) == 2
        #Create the dynamic generator
        case_gen = DynamicGenerator(
            get_name(g),
            1.0, # ω_ref,
            AF_machine(), #machine
            shaft_no_damping(), #shaft
            avr_type1(), #avr
            tg_none(), #tg
            pss_none(), #pss
        )
        #Attach the dynamic generator to the system by
        # specifying the dynamic and static components
        add_component!(sys, case_gen, g)
    #Find the generator at bus 103
    elseif get_number(get_bus(g)) == 3
        #Create the dynamic inverter
        case_inv = DynamicInverter(
            get_name(g),
            1.0, # ω_ref,
            converter_high_power(), #converter
            GFL_outer_control(), #outer control
            GFL_inner_control(), #inner control voltage source
            dc_source_lv(), #dc source
            pll(), #pll
            lcl_filt(), #filter
        )
        #Attach the dynamic inverter to the system
        add_component!(sys, case_inv, g)
    end
end

to_json(sys, joinpath(pwd(), "../json_data/9bus_VSM_SM_GFL_.json"), force = true)