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
pygui(true)

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

add_lines_sweep!(gss, [line_params], line_adders)
add_zipe_sweep!(gss, missing, (x->LoadParams(x...)).(gridsearch())) # no standard load adder. already in the system.

add_result!(gss, "Eigenvalues", get_eigenvalues)
add_result!(gss, "Eigenvectors", get_eigenvectors)
add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
add_result!(gss, "Simulation Status", get_sim_status)
add_result!(gss, "Error", get_error)
add_result!(gss, "sim_filename", get_serialized_sim_filename_function_builder("data/saved_objects"))
# add_result!(gss, "sim.solution", get_sim)

executeSims!(gss, BranchTrip(0.5, ACBranch, line_params.alg_line_name), (0.48, 0.55), 0.001, 0.001, true, mktempdir())
expand_columns!(gss)
save_data(gss, "data/results.tsv")


df = load_data("data/results.tsv")
names(df)
# df = gss.df

# Eigenvalue Plot
function eigplot()
    eigs = [
        @subset(df, :z_percent.==1.0).Eigenvalues,
        @subset(df, :z_percent.==0.5).Eigenvalues,
        @subset(df, :z_percent.==0.2, :p_percent.==0.5).Eigenvalues,
        @subset(df, :z_percent.==0.2, :p_percent.==0.5).Eigenvalues
    ]

    for i in 1:length(eigs)
        plt.scatter(
            real.(reduce(vcat, eigs[i])), 
            imag.(reduce(vcat, eigs[i])),
            label=string(zipe_combos[i]),
            s=5,
            alpha=0.25
        )
    end
    plt.xlabel(L"\mathrm{Re}(\lambda)")
    plt.ylabel(L"\mathrm{Im}(\lambda)")
    plt.legend()
    plt.show()
end

# Max Eig box
function maxeigbox()
    sta = map(x->maximum(real.(x)), @subset(df, :"Line Model" .== "statpi").Eigenvalues)
    dyn  = map(x->maximum(real.(x)), @subset(df, :"Line Model" .== "dynpi").Eigenvalues)
    print(dyn)
    (fig, axs) = plt.subplots(2, 1, sharex=true)
    axs[1].boxplot(dyn, vert=false)
    axs[2].boxplot(sta, vert=false)
    axs[1].set_title("dynpi")
    axs[2].set_title("statpi")
    fig.suptitle(L"\max_{i}\: \mathrm{Re}(\lambda_i)")
    plt.show()
end

function transient()
    data = @subset df begin
        # :"Load Voltage at Bus 5"
        :"Line Model" .== "statpi"
        :"injector at {Bus 3}" .== "GFM"
        :"injector at {Bus 2}" .== "SM"
        :"injector at {Bus1}" .== "GFM"
    end
    plt.plot(0.48:0.0001:0.55, data.var"Load Voltage at Bus 5"[1])
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage")
    plt.show()
end