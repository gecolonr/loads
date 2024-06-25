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
# using PyPlot
# const plt = PyPlot
# pygui(true)

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")


##################################################################
####################### MAKE BASE SYSTEM #########################
##################################################################

# system from 9 bus base system
s = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"), time_series_in_memory=true)
set_runchecks!(s, false)
##################################################################
############ MACHINES FOR INV/GEN SWEEP ##########################
##################################################################

# functions to create our machines: one inverter and one generator
gfm_inj() = DynamicInverter(
    "GFM", # Grid forming control
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)
gfl_inj() = DynamicInverter(
    "GFL", # Grid following inverter
    1.0, # ω_ref,
    converter_high_power(), #converter
    GFL_outer_control(), #outer control
    GFL_inner_control(), # this one's GFL
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)
sm_inj() = DynamicGenerator(
    "SM", # Synchronous Machine
    1.0, # ω_ref,
    AF_machine(), #machine
    shaft_no_damping(), #shaft
    avr_type1(), #avr
    tg_none(), #tg
    pss_none(), #pss
)

##################################################################
############ LINE PARAMS FOR TLMODELS SWEEP ######################
##################################################################

# taken from TLModels.jl `TLmodels_tutorial.ipynb`
impedance_csv = "../TLModels.jl/data/cable_data/dommel_data.csv"
capacitance_csv = "../TLModels.jl/data/cable_data/dommel_data_C.csv"

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
    "statpi (default)"=>no_change,
    "dynpi (default)"=>create_simple_dynpi,
    # "MSSB"=>create_MSSB_system,
])

##################################################################
############ ZIPE STUFF ##########################################
##################################################################

"""
returns generator of all Z, I, P, and E combinations in steps of `dx`.

NOTE: previously this code did i∈[0, 1], j∈[0, 1-i], k∈[0, 1-i-j], 
which makes the if statement at the end implicit, makes the code cleaner, 
and makes it align more with the typical mathematical definition. However,
it turns out that sometimes float math would give you 1:dx:(5*dx-ε), which
would return one less element than expected. This is BAD since we'd lose 
quite a few elements.
The current code works much better. The `dx/10` is just very safe tolerance.
"""
function gridsearch(dx=0.1)
    return ([i j k 1.0-i-j-k] for i in 0:dx:1 for j in 0:dx:1 for k in 0:dx:1 if ((i+j+k-1.0) < (dx/10)) )
end


zipe_combos = [
#     Z    I    P    E
    # [1.0, 0.0, 0.0, 0.0],
    # [0.5, 0.1, 0.2, 0.2],
    # [0.2, 0.1, 0.5, 0.2],
    # [0.2, 0.1, 0.2, 0.5],
    # [0.2, 0.1, 0.7, 0.0],
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.3, 0.3, 0.3, 0.1],
    [0.2, 0.2, 0.2, 0.4],
    [0.1, 0.1, 0.1, 0.7],
    [0.0, 0.0, 0.0, 1.0],
]
η_combos = Dict(
    "Constant Impedance" => [1.0, 0.0, 0.0, 0.0],
    "Constant Power"     => [0.0, 0.0, 1.0, 0.0],
    "Full E Load"        => [0.0, 0.0, 0.0, 1.0],
    "High E Load"        => [0.1, 0.1, 0.1, 0.7],
    "Medium E Load"      => [0.2, 0.2, 0.2, 0.4],
    "Low E Load"         => [0.3, 0.3, 0.3, 0.1],
    "Low P High E Load"  => [0.15, 0.15, 0.15, 0.55],
    "Low P Low E Load"   => [0.15, 0.15, 0.55, 0.15],

    "ZI-0PE" => [0.5, 0.5, 0.0, 0.0],
    "ZI-0.1P" => [0.45, 0.45, 0.1, 0.0],
    "ZI-0.1E" => [0.45, 0.45, 0.0, 0.1],
    "ZI-0.2P" => [0.4, 0.4, 0.2, 0.0],
    "ZI-0.2E" => [0.4, 0.4, 0.0, 0.2],
    "ZI-0.3P" => [0.35, 0.35, 0.3, 0.0],
    "ZI-0.3E" => [0.35, 0.35, 0.0, 0.3],
    "ZI-0.4P" => [0.3, 0.3, 0.4, 0.0],
    "ZI-0.4E" => [0.3, 0.3, 0.0, 0.4],
    "ZI-0.5P" => [0.25, 0.25, 0.5, 0.0],
    "ZI-0.5E" => [0.25, 0.25, 0.0, 0.5],
    "ZI-0.6P" => [0.2, 0.2, 0.6, 0.0],
    "ZI-0.6E" => [0.2, 0.2, 0.0, 0.6],
    "ZI-0.7P" => [0.15, 0.15, 0.7, 0.0],
    "ZI-0.7E" => [0.15, 0.15, 0.0, 0.7],
    "ZI-0.8P" => [0.1, 0.1, 0.8, 0.0],
    "ZI-0.8E" => [0.1, 0.1, 0.0, 0.8],
    "ZI-0.9P" => [0.05, 0.05, 0.9, 0.0],
    "ZI-0.9E" => [0.05, 0.05, 0.0, 0.9],
    # "ZI-1.0E" => [0.0, 0.0, 1.0, 0.0],
    # "ZI-1.0E" => [0.0, 0.0, 0.0, 1.0]
)

# η_combos = vcat((x->[round.([(1-x)/2, (1-x)/2, x, 0], digits=3), 
#                      round.([(1-x)/2, (1-x)/2, 0, x], digits=3)]).(0:0.1:1)...)

# η_combos = [round.([0.1, 0.1, x, 0.8-x], digits=3) for x in (0:0.1:0.8)]
# η_combos = [round.([0.2, 0.2, x, 0.6-x], digits=3) for x in (0:0.1:0.8)]

##################################################################
###################### POWER SETPOINT ############################
##################################################################

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

##################################################################
################### RUNNING ALL OF IT ############################
##################################################################

# get all combinations of generators on this system
# gss = GridSearchSys(s, [gfl_inj(), gfm_inj(), sm_inj()])

"""
function to get small signal analysis 
"""
function small_signal_tripped(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
    sm = small_signal_analysis(newsim)
    return sm
end


gss = GridSearchSys(s, [sm_inj() gfl_inj() gfm_inj(); 
                        sm_inj() gfm_inj() gfl_inj(); 
                        gfm_inj() gfl_inj() sm_inj();
                        gfm_inj() gfm_inj() sm_inj();
                        gfm_inj() sm_inj() sm_inj();
                        sm_inj() sm_inj() sm_inj();
                        gfm_inj() gfm_inj() gfl_inj();
                        gfm_inj() gfl_inj() gfm_inj();],
                        ["Bus1", "Bus 2", "Bus 3"]) # just make sure the busses are in the right order
set_chunksize!(gss, 200)

add_generic_sweep!(gss, "Power Setpoint", set_power_setpt!, collect(0.2:0.1:1.4))
add_lines_sweep!(gss, [line_params], line_adders)
add_zipe_sweep!(gss, missing, (x->LoadParams(x...)).(values(η_combos))) # no standard load adder. already in the system.
add_result!(gss, "initial_sm", get_sm)
add_result!(gss, "final_sm", small_signal_tripped)
# add_result!(gss, "Eigenvalues", get_eigenvalues)
# add_result!(gss, "Eigenvectors", get_eigenvectors)
add_result!(gss, "time", get_time)
add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)


execute_sims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), tspan=(0.48, 5.5), dtmax=0.05, run_transient=true, log_path="/data/reiddye/loads/fivesecondsbetter")
