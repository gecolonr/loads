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
using Sundials
using ForwardDiff

pygui(true)


include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

# one or the other
# you have to run the sims first to get these data files
# since they can be big files we shouldn't push them to git
df = load_serde_data("data/results.jls")
df = load_serde_data("data/sims")

expand_columns!(df)

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
    axs[1].set_xscale("log")
    fig.suptitle(L"\max_{i}\: \mathrm{Re}(\lambda_i)")
    plt.show()
end

function transient()
    dynpidata = @subset df begin
        # :"Load Voltage at Bus 5"
        :"Line Model" .== "dynpi"
        :"injector at {Bus 3}" .== "GFM"
        :"injector at {Bus 2}" .== "SM"
        :"injector at {Bus1}" .== "GFM"
    end
    statpidata = @subset df begin
        # :"Load Voltage at Bus 5"
        :"Line Model" .== "statpi"
        :"injector at {Bus 3}" .== "GFM"
        :"injector at {Bus 2}" .== "SM"
        :"injector at {Bus1}" .== "GFM"
    end

    # busses = ["Load Voltage at $i" for i in get_name.(get_components(Bus, df.sim[1].sys))]
    busses = ["Load Voltage at $i" for i in get_name.(get_bus.(get_components(StandardLoad, first(df.sim).sys)))]
    fig, axs = plt.subplots(2, length(busses), sharex=true)
    fig.suptitle("Load Voltage at Bus 5 for BranchTrip on line Bus 5-Bus 4-i₁")
    axs[1].set_ylabel("dynpi")
    axs[2].set_ylabel("statpi")
    for (j, bus) in enumerate(busses)
        for row in eachrow(dynpidata)
            if row[bus] isa Missing
                println("dynpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                continue
            end
            axs[1, j].plot(0.48:0.0001:0.55, row[bus], label="Z=$(round(row.z_percent, digits=3)), I=$(round(row.i_percent, digits=3)), P=$(round(row.p_percent, digits=3)), E=$(round(row.e_percent, digits=3))")
        end
        for row in eachrow(statpidata)
            if row[bus] isa Missing
                println("statpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                continue
            end
            axs[2, j].plot(0.48:0.0001:0.55, row[bus], label="Z=$(round(row.z_percent, digits=3)), I=$(round(row.i_percent, digits=3)), P=$(round(row.p_percent, digits=3)), E=$(round(row.e_percent, digits=3))")
        end
        axs[1, j].set_title(bus)
        axs[2, j].set_xlabel("Time (s)")
        
        axs[1, j].legend()
        axs[2, j].legend()
    end
    plt.xlabel("Time (s)")
    # axs[1].set_ylabel("Voltage")
    # axs[2].set_ylabel("Voltage")
    plt.show()
end
