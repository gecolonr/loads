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

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")


##################################################################
################## MAKE 2-MACHINE BASE SYSTEM ####################
##################################################################

# system from OMIB base system
s = System(joinpath(pwd(), "data/raw_data/OMIB.raw"))

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

##################################################################
############ MACHINES FOR INV/GEN SWEEP ##########################
##################################################################

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

##################################################################
############ LINE PARAMS FOR TLMODELS SWEEP ######################
##################################################################

# taken from TLModels.jl `TLmodels_tutorial.ipynb`
impedance_csv = "../TLModels.jl/data/cable_data/dommel_data.csv"
capacitance_csv = "../TLModels.jl/data/cable_data/dommel_data_C.csv"

M = 3
z_km, y_km, z_km_ω, Z_c = get_line_parameters_from_data(impedance_csv, capacitance_csv, M)

line_length_dict = Dict(map(x->x.name, get_components(Line, s)) .=> 100)

line_params = LineModelParams(
    z_km, 
    y_km, 
    z_km_ω, 
    Z_c,
    M,
    line_length_dict,    
    "otherline", # this is never used i think
    10.0,
    1.0,
    1.0
)

line_adders = Dict{String, Function}([
    "statpi"=>create_statpi_system,
    "dynpi"=>create_dynpi_system,
    "MSSB"=>create_MSSB_system,
])

##################################################################
############ ZIPE STUFF ##########################################
##################################################################

load(s) = StandardLoad(
    name="load",
    available=true,
    bus=first(get_components(ACBus, s)),
    base_power=100.0,
    constant_active_power=2.0,
    constant_reactive_power=0.1,
)

"""
returns generator of all Z, I, P, and E combinations in steps of `dx`.
"""
function gridsearch(dx=0.1, Zmax=1.0, Imax=1.0, Pmax=1.0)
    return ([i j k 1.0-i-j-k] for i in 0.0:dx:Zmax for j in 0.0:dx:min(1.0-i, Imax) for k in 0.0:dx:min(1.0-i-j, Pmax))
end

##################################################################
################### RUNNING ALL OF IT ############################
##################################################################

# get all combinations of generators on this system
gss = GridSearchSys(s, [case_inv(), case_gen()])
add_lines_sweep!(gss, [line_params], line_adders)
add_zipe_sweep!(gss, load, map(x->LoadParams(x...), gridsearch()))

results_df = executeSims(gss, BranchTrip(0.5, ACBranch, "otherline"), (0.0, 2.0))

expand!(results_df)

CSV.write("data/results.tsv", results_df, delim='\t')