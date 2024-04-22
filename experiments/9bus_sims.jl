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
using PyPlot
const plt = PyPlot
# pygui(true)

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")


##################################################################
################## MAKE 2-MACHINE BASE SYSTEM ####################
##################################################################

# system from 9 bus base system
s = System(joinpath(pwd(), "data/raw_data/WSCC_9bus.raw"))

##################################################################
############ MACHINES FOR INV/GEN SWEEP ##########################
##################################################################

# functions to create our machines: one inverter and one generator
gfm_inj() = DynamicInverter(
    "GFM",
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFM_inner_control(), #inner control voltage source
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)
gfl_inj() = DynamicInverter(
    "GFL",
    1.0, # ω_ref,
    converter_high_power(), #converter
    VSM_outer_control(), #outer control
    GFL_inner_control(), # this one's GFL
    dc_source_lv(), #dc source
    pll(), #pll
    filt(), #filter
)
sm_inj() = DynamicGenerator(
    "SM", # stands for "Generator"
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

line_adders = Dict{String, Function}([
    "statpi"=>create_statpi_system,
    "dynpi"=>create_dynpi_system,
    # "MSSB"=>create_MSSB_system,
])

##################################################################
############ ZIPE STUFF ##########################################
##################################################################

"""
returns generator of all Z, I, P, and E combinations in steps of `dx`.
"""
function gridsearch(dx=0.1, Zmax=1.0, Imax=1.0, Pmax=1.0)
    return ([i j k 1.0-i-j-k] for i in 0.0:dx:Zmax for j in 0.0:dx:min(1.0-i, Imax) for k in 0.0:dx:min(1.0-i-j, Pmax))
end

zipe_combos = [
#     Z    I    P    E
    [1.0, 0.0, 0.0, 0.0],
    [0.5, 0.1, 0.2, 0.2],
    [0.2, 0.1, 0.5, 0.2],
    [0.2, 0.1, 0.2, 0.5]
]

##################################################################
################### RUNNING ALL OF IT ############################
##################################################################

# get all combinations of generators on this system
gss = GridSearchSys(s, [gfl_inj(), gfm_inj(), sm_inj()])
set_chunksize(gss, 20)

add_lines_sweep!(gss, [line_params], line_adders)
add_zipe_sweep!(gss, missing, (x->LoadParams(x...)).(zipe_combos)) # no standard load adder. already in the system.

add_result!(gss, "Eigenvalues", get_eigenvalues)
add_result!(gss, "Eigenvectors", get_eigenvectors)
add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
add_result!(gss, "Simulation Status", get_sim_status)
add_result!(gss, "Error", get_error)
add_result!(gss, "sim", get_sim)

executeSims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), (0.48, 0.55), 0.001, 0.0001, true)
# expand_columns!(gss)
# save_serde_data(gss, "data/results.jls")

