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

# currently unused
function no_change(sys::System, params::LineModelParams) return sys end
# currently unused
function create_simple_dynpi(sys::System, params::LineModelParams)
    for ll in get_components(Line, sys)
        if (ll.name != "Bus 5-Bus 4-i_1")
            dyn_branch = DynamicBranch(ll)
            add_component!(sys, dyn_branch)
        end
    end
    return sys
end

# add more line models if you want. adds more runtime.
line_adders = Dict{String, Function}([
    "statpi"=>create_statpi_system,
    "dynpi"=>create_dynpi_system,
    # "statpi (default)"=>no_change,
    # "dynpi (default)"=>create_simple_dynpi,
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

# add whatever you want here. 
# currenlt set is for plots in paper and will work with
# `make_plots_from_paper`. For `make_plots_more_data`,
# use everything except the first four.
# 
# If you plan to make your own plots, you can do whatever
# you like
η_combos = [
    [1.0, 0.0, 0.0, 0.0],

    [0.5, 0.5, 0.0, 0.0],
    # [0.45, 0.45, 0.1, 0.0],
    # [0.45, 0.45, 0.0, 0.1],
    [0.4, 0.4, 0.2, 0.0],
    [0.4, 0.4, 0.0, 0.2],
    # [0.35, 0.35, 0.3, 0.0],
    # [0.35, 0.35, 0.0, 0.3],
    [0.3, 0.3, 0.4, 0.0],
    [0.3, 0.3, 0.0, 0.4],
    # [0.25, 0.25, 0.5, 0.0],
    # [0.25, 0.25, 0.0, 0.5],
    [0.2, 0.2, 0.6, 0.0],
    [0.2, 0.2, 0.0, 0.6],
    # [0.15, 0.15, 0.7, 0.0],
    # [0.15, 0.15, 0.0, 0.7],
    [0.1, 0.1, 0.8, 0.0],
    [0.1, 0.1, 0.0, 0.8],
    # [0.05, 0.05, 0.9, 0.0],
    # [0.05, 0.05, 0.0, 0.9],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0]
]

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

#* For more information on how the `GridSearchSys` struct works, See
#* the [guide](https://reid.xz.ax/loads/gsstutorial). The docstrings
#* are helpful too.

"""
function to get small signal analysis after the branch trip.

Not currently used.
"""
function small_signal_tripped(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
    sm = small_signal_analysis(newsim)
    return sm
end

# enable or add more cases here if you want. The second line is the case from the paper.
gss = GridSearchSys(s, [#sm_inj() gfl_inj() gfm_inj(); 
                        sm_inj() gfm_inj() gfl_inj(); 
                        sm_inj() sm_inj() sm_inj();
                        # gfm_inj() gfl_inj() sm_inj();
                        # gfm_inj() gfm_inj() sm_inj();
                        # gfm_inj() sm_inj() sm_inj();
                        # sm_inj() sm_inj() sm_inj();
                        # gfm_inj() gfm_inj() gfl_inj();
                        # gfm_inj() gfl_inj() gfm_inj();
                        ],
                        ["Bus1", "Bus 2", "Bus 3"]) # just make sure the busses are in the right order
set_chunksize!(gss, 200)

# we swept power setpoint with steps of 0.1 and 0.2. Your choice. 0.1
# will double the data size but it makes the slider on the plots smoother.
add_generic_sweep!(gss, "Power Setpoint", set_power_setpt!, collect(0.05:0.05:1.0))
add_lines_sweep!(gss, [line_params], line_adders)
add_zipe_sweep!(gss, missing, (x->LoadParams(x...)).(η_combos))
# add_generic_sweep!(gss, "ZIPE Load Params", create_ZIPE_load, (x->LoadParams(x...)).(η_combos))

add_result!(gss, "Eigenvalues", get_eigenvalues)
add_result!(gss, ["Bus 3 Injector Current", "Bus 1 Injector Current", "Bus 2 Injector Current"], get_injector_currents)
add_result!(gss, "time", get_time)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)


execute_sims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), tspan=(0.48, 1.0), dtmax=0.05, run_transient=true, log_path="data/paper_results")
