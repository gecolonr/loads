cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
using InfrastructureSystems
using Plots
using PlotlyJS, DataFrames
using TLmodels
using CSV
using DataFramesMeta
using LaTeXStrings
using Logging
Logging.disable_logging(Logging.Error)

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

s = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"))

gfm_inj() = DynamicInverter(
    "I", # stands for "Inverter"
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)

gfl_inj() = DynamicInverter(
    "I", # stands for "Inverter"
    1.0, # ω_ref,
    converter_high_power(), #converter
    GFL_outer_control(), #outer control
    GFL_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)

sm_inj() = DynamicGenerator(
    "G", # stands for "Generator"
    1.0, # ω_ref,
    AF_machine(), #machine
    shaft_no_damping(), #shaft
    avr_type1(), #avr
    tg_none(), #tg
    pss_none(), #pss
)

# taken from TLModels.jl `TLmodels_tutorial.ipynb`
impedance_csv = "data/cable_data/dommel_data.csv"
capacitance_csv = "data/cable_data/dommel_data_C.csv"

M = 3
z_km, y_km, z_km_ω, Z_c = get_line_parameters_from_data(impedance_csv, capacitance_csv, M)

line_length_dict = Dict(
    "Bus 5-Bus 4-i_1" => 90,
    "Bus 7-Bus 8-i_1" => 80,
    "Bus 6-Bus 4-i_1" => 100,
    "Bus 7-Bus 5-i_1" => 170,
    "Bus 8-Bus 9-i_1" => 110,
    "Bus 9-Bus 6-i_1" => 180,
)

line_params = LineModelParams(
    z_km, 
    y_km, 
    z_km_ω, 
    Z_c,
    M,
    line_length_dict,    
    get_name(first(get_components(Line, s))),
    10.0,
    1.0,
    1.0
)

function no_change(sys::System, params::LineModelParams)
    return sys
end

function create_simple_dynpi(sys::System, params::LineModelParams)
    for ll in get_components(Line, sys)
        if (ll.name != "Bus 5-Bus 4-i_1")
            dyn_branch = DynamicBranch(ll)
            add_component!(sys, dyn_branch)
        end
    end
    return sys
end

line_adders = Dict{String, Function}([
    "statpi (dommel)"=>create_statpi_system,
    "dynpi (dommel)"=>create_dynpi_system,
    # "MSSB"=>create_MSSB_system,
])

η_combos = Dict(
    "Constant Impedance" => [1.0, 0.0, 0.0, 0.0],
)

function set_power_setpt!(sys::System, scale::Real)
    for load in get_components(StandardLoad, sys)
        set_impedance_active_power!(load, get_impedance_active_power(load)*scale)
        set_current_active_power!(load, get_current_active_power(load)*scale)
        set_constant_active_power!(load, get_constant_active_power(load)*scale)
        
        set_impedance_reactive_power!(load, get_impedance_reactive_power(load)*scale)
        set_current_reactive_power!(load, get_current_reactive_power(load)*scale)
        set_constant_reactive_power!(load, get_constant_reactive_power(load)*scale)
    end
    # if scale <= 1.0 return sys end
    for gen in get_components(Generator, sys)
        if gen.bus.bustype == ACBusTypes.PV
            # set_base_power!(g, g.base_power * p.load_scale)
            set_active_power!(gen, get_active_power(gen) * scale)
            set_reactive_power!(gen, get_reactive_power(gen) * scale)
        end
    end
    return sys
end

### Defining a function to scale the line impedance
function scale_line_impedance!(sys::System, scale::Real)
    for line in get_components(Line, sys)
        if line.name == "Bus 5-Bus 4-i_1" || line.name == "Bus 6-Bus 4-i_1"
            line.r = line.r * scale
            line.x = line.x * scale
        end
    end
    return sys
end

# Defining a function to change generation portfolio, maybe this needs to be done in two
function change_gen_portfolio!(sys::System, gen_scale::Real)
    p_sum = 0
    q_sum = 0
    gen_num = 0
    for gen in get_components(Generator, sys)
        p_sum += get_active_power(gen)
        q_sum += get_reactive_power(gen)
        gen_num += 1
    end
    for gen in get_components(Generator, sys)
        # if (typeof(gen.dynamic_injector) == DynamicGenerator{})
        if (gen.bus. number == 1)  
            set_active_power!(gen, gen_scale/gen_num*p_sum)
            set_reactive_power!(gen, gen_scale/gen_num*q_sum)
        end
    end
    return sys
end

function change_inv_portfolio!(sys::System, gfm_scale::Real)
    p_inv = 0
    q_inv = 0
    for gen in get_components(Generator, sys)
        # if (typeof(gen.dynamic_injector) == DynamicInverter{})
        if (gen.bus.number == 2 || gen.bus.number == 3)
            p_inv += get_active_power(gen)
            q_inv += get_reactive_power(gen)
        end
    end
    for gen in get_components(Generator, sys)
        # if (typeof(gen.dynamic_inverter.inner_control) == VoltageModeControl)
        if (gen.bus.number == 2)    
            set_active_power!(gen, gfm_scale*p_inv)
            set_reactive_power!(gen, gfm_scale*q_inv)      
        # elseif (typeof(gen.dynamic_inverter.inner_control) == CurrentModeControl)
        elseif (gen.bus. number == 3)
            set_active_power!(gen, (1 - gfm_scale)*p_inv)
            set_reactive_power!(gen, (1 - gfm_scale)*q_inv)
        else
            continue
        end
    end
    return sys
end

# function small_signal_tripped(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
#     if isnothing(sim) return missing end
#     sys = deepcopy(sim.sys)
#     remove_component!(Line, sys, sim.perturbations[1].branch_name)
#     newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
#     sm = small_signal_analysis(newsim)
#     return sm
# end

gss = GridSearchSys(s, [sm_inj() gfl_inj() gfm_inj();],
                        ["Bus1", "Bus 2", "Bus 3"]) # just make sure the busses are in the right order
set_chunksize!(gss, 200)

power_stpt_range = collect(1.0:0.1:1.0)
scale_impedance_range = collect(1.0:0.1:1.0)
gen_portfolio_range = collect(0.1:0.1:0.3)
gfm_portfolio_range = collect(0.2:0.1:0.6)

add_generic_sweep!(gss, "Power Setpoint", set_power_setpt!, power_stpt_range)
add_generic_sweep!(gss, "Line impedance increase", scale_line_impedance!, scale_impedance_range)
add_generic_sweep!(gss, "SM percent", change_gen_portfolio!, gen_portfolio_range)
add_generic_sweep!(gss, "GFM percent", change_inv_portfolio!, gen_portfolio_range)
add_lines_sweep!(gss, [line_params], line_adders)
# add_zipe_sweep!(gss, missing, (x->LoadParams(x...)).(values(η_combos))) # no standard load adder. already in the system.
# add_result!(gss, "initial_sm", get_sm)
# add_result!(gss, "final_sm", small_signal_tripped)

# add_result!(gss, "time", get_time)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
add_result!(gss, "Eigs", get_eigenvalues)
execute_sims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), tspan=(0.48, 5.5), dtmax=0.05, run_transient=true, log_path="data/gab_tests")
