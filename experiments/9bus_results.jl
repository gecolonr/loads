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
using SciMLBase
using ModelingToolkit

using ForwardDiff
using SparseArrays
using SymbolicIndexingInterface

pygui(true)

include("sysbuilder.jl")
include("../data/build_scripts/device_models.jl")

df = load_serde_data("data/results.jls")


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
