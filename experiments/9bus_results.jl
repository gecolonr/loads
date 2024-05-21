cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
using InfrastructureSystems
using PlotlyJS
using DataFrames
using TLmodels
using DataFramesMeta
using LaTeXStrings
using PyPlot
const plt = PyPlot
# using PyCall
# mplc = pyimport("mplcursors")
# pygui(:qt5)
pygui(true)


include("sysbuilder.jl")

# one or the other
# you have to run the sims first to get these data files
# since they can be big files we shouldn't push them to git
# df = load_serde_data("data/results.jls")
# df = load_serde_data("data/sims")

# expand_columns!(df)

# Eigenvalue Plot
function eigplot(df)
    df2 = @subset df begin
        :"injector at {Bus 3}" .== "GFM"
        :"injector at {Bus1}" .== "SM"
        :"injector at {Bus 2}" .== "GFL"
    end

    # select loads:
    #  1. 0. 0. 0
    #  0. 0. 1. 0
    #  0. 0. 0. 1
    #  .3 .3 .3 .1
    #  .2 .2 .2 .4
    #  .1 .1 .1 .7
    combos = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.3, 0.3, 0.3, 0.1],
        [0.2, 0.2, 0.2, 0.4],
        [0.1, 0.1, 0.1, 0.7],
        [0.0, 0.0, 0.0, 1.0],
    ]
    geteigs(df) = [
        @subset(df, :z_percent.≈1.0                                    ).Eigenvalues,
        @subset(df,                                   :p_percent.≈1.0).Eigenvalues,
        @subset(df, :z_percent.≈0.3, :i_percent.≈0.3, :p_percent.≈0.3).Eigenvalues,
        @subset(df, :z_percent.≈0.2, :i_percent.≈0.2, :p_percent.≈0.2).Eigenvalues,
        @subset(df, :z_percent.≈0.1, :i_percent.≈0.1, :p_percent.≈0.1).Eigenvalues,
        @subset(df, :z_percent.≈0.0, :i_percent.≈0.0, :p_percent.≈0.0).Eigenvalues,
    ]
    stateigs = geteigs(@subset(df2, :"Line Model" .== "statpi"))
    dyneigs = geteigs(@subset(df2, :"Line Model" .== "dynpi"))
    fig, axs = plt.subplots(2, 1, sharex=true)
    axs[1].set_title("Static Pi Model")
    axs[2].set_title("Dynamic Pi Model")
    for i in 1:length(combos)
        axs[1].scatter(
            real.(reduce(vcat, stateigs[i])), 
            imag.(reduce(vcat, stateigs[i])),
            label=string(combos[i]),
            s=5,
            # alpha=0.25
        )
        axs[2].scatter(
            real.(reduce(vcat, dyneigs[i])), 
            imag.(reduce(vcat, dyneigs[i])),
            label=string(combos[i]),
            s=5,
            # alpha=0.25
        )
    end
    axs[2].set_xlabel(L"\mathrm{Re}(\lambda)")
    axs[1].set_ylabel(L"\mathrm{Im}(\lambda)")
    axs[2].set_ylabel(L"\mathrm{Im}(\lambda)")
    axs[1].legend(prop=Dict("size"=>5))
    axs[2].legend(prop=Dict("size"=>5))
    axs[1].axvline(x=0, color="black", ls="--", lw=1)
    axs[2].axvline(x=0, color="black", ls="--", lw=1)
    plt.show()
end

# Max Eig box
function maxeigbox(df)
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

function transient(df)
    df2 = @subset df begin
        :"injector at {Bus 3}" .== "GFM"
        :"injector at {Bus1}" .== "SM"
        :"injector at {Bus 2}" .== "GFL"
    end
    combos = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.3, 0.3, 0.3, 0.1],
        [0.2, 0.2, 0.2, 0.4],
        [0.1, 0.1, 0.1, 0.7],
        [0.0, 0.0, 0.0, 1.0],
    ]
    selectZIPE(df) = vcat([
        @subset(df, :z_percent.≈1.0                                  ),
        @subset(df,                                   :p_percent.≈1.0),
        @subset(df, :z_percent.≈0.3, :i_percent.≈0.3, :p_percent.≈0.3),
        @subset(df, :z_percent.≈0.2, :i_percent.≈0.2, :p_percent.≈0.2),
        @subset(df, :z_percent.≈0.1, :i_percent.≈0.1, :p_percent.≈0.1),
        @subset(df, :z_percent.≈0.0, :i_percent.≈0.0, :p_percent.≈0.0),
    ]...)
    statpidata = (@subset(df2, :"Line Model" .== "statpi"))
    dynpidata = (@subset(df2, :"Line Model" .== "dynpi"))

    # return statpidata

    # busses = ["Load Voltage at $i" for i in get_name.(get_components(Bus, df.sim[1].sys))]
    busses = ["Load Voltage at $i" for i in get_name.(get_bus.(get_components(StandardLoad, first(df.sim).sys)))]
    fig, axs = plt.subplots(2, length(busses), sharex=true)
    fig.suptitle("Load Voltages for BranchTrip on line Bus 5-Bus 4-i₁")
    axs[1].set_ylabel("dynpi")
    axs[2].set_ylabel("statpi")
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "b", "g", "r", "c", "m", "y", "k", "w"]
    colordict = Dict()
    # axs = [Plot() Plot() Plot(); Plot() Plot() Plot()]
    for (j, bus) in enumerate(busses)
        for row in eachrow(dynpidata)
            if row[bus] isa Missing
                println("dynpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                continue
            end

            params = (row.z_percent, row.i_percent, row.p_percent, row.e_percent)
            if !(params in keys(colordict))
                colordict[params] = colors[length(colordict)+1]
            end
            # addtraces!(axs[1, j], PlotlyJS.scatter(; x=0.48:0.0001:0.55, y=row[bus]))
            axs[1, j].plot(0.48:0.00005:0.55, row[bus], color=colordict[params], label="Z=$(round(row.z_percent, digits=3)), I=$(round(row.i_percent, digits=3)), P=$(round(row.p_percent, digits=3)), E=$(round(row.e_percent, digits=3))")
        end
        for row in eachrow(statpidata)
            if row[bus] isa Missing
                println("statpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                continue
            end

            params = (row.z_percent, row.i_percent, row.p_percent, row.e_percent)
            if !(params in keys(colordict))
                colordict[params] = colors[length(colordict)+1]
            end
            # addtraces!(axs[2, j], PlotlyJS.scatter(; x=0.48:0.0001:0.55, y=row[bus]))
            axs[2, j].plot(0.48:0.00005:0.55, row[bus], color=colordict[params], label="Z=$(round(row.z_percent, digits=3)), I=$(round(row.i_percent, digits=3)), P=$(round(row.p_percent, digits=3)), E=$(round(row.e_percent, digits=3))")
        end
        axs[1, j].set_title(bus)
        axs[2, j].set_xlabel("Time (s)")
        
        axs[1, j].legend()#prop=Dict("size"=>5))
        axs[2, j].legend()#prop=Dict("size"=>5))
    end
    plt.xlabel("Time (s)")
    # axs[1].set_ylabel("Voltage")
    # axs[2].set_ylabel("Voltage")
    plt.show()
    plt.tight_layout()
end

function sillies(df)
    normal = @subset df begin
        :"Line Model" .== "statpi"
        ((x)->(x isa Missing ? false : (x[end]>0.5))).(:"Load Voltage at Bus 5")
    end
    bad = @subset df begin
        :"Line Model" .== "statpi"
        ((x)->(x isa Missing ? false : (x[end]<0.5))).(:"Load Voltage at Bus 5")
    end
    missingdata = @subset df begin
        :"Line Model" .== "statpi"
        (x->(x isa Missing)).(:"Load Voltage at Bus 5")
    end
    data = [bad, normal, missingdata]
    data = [((x->x.z_percent).(data)),
            ((x->x.i_percent).(data)),
            ((x->x.p_percent).(data)),
            ((x->x.e_percent).(data))]
    fig, axs = plt.subplots(4, 1, sharex=true)
    weights = x->((n->(2/(10*(1.1-n)*(10*(1.1-n)+1)))).(x))
    axs[1, 1].hist(data[1], histtype="bar", stacked=true, bins=0:0.1:1.1 .- 0.01, weights=weights.(data[1]), label=["anomalous", "normal", "convergence failure"])
    axs[2, 1].hist(data[2], histtype="bar", stacked=true, bins=0:0.1:1.1 .- 0.01, weights=weights.(data[2]), label=["anomalous", "normal", "convergence failure"])
    axs[3, 1].hist(data[3], histtype="bar", stacked=true, bins=0:0.1:1.1 .- 0.01, weights=weights.(data[3]), label=["anomalous", "normal", "convergence failure"])
    axs[4, 1].hist(data[4], histtype="bar", stacked=true, bins=0:0.1:1.1 .- 0.01, weights=weights.(data[4]), label=["anomalous", "normal", "convergence failure"])
    axs[1, 1].hist
    fig.legend()
    map(x->x.set_yscale("log"), axs)
    plt.show()
end

function endvoltage(df)
    statpi = @subset df begin
        :"Line Model" .== "statpi"
        ((x)->!(x isa Missing)).(:"Load Voltage at Bus 5")
    end
    dynpi = @subset df begin
        :"Line Model" .== "dynpi"
        ((x)->!(x isa Missing)).(:"Load Voltage at Bus 5")
    end
    fig, axs = plt.subplots(3, sharex=false, sharey=true)
    axs[1].hist([(x->(x[end]-x[1])).(statpi.var"Load Voltage at Bus 5"),
                 (x->(x[end]-x[1])).(dynpi.var"Load Voltage at Bus 5")], label=[raw"Static $\pi$ model", raw"Dynamic $\pi$ model"], bins=100)
    axs[2].hist([(x->(x[end]-x[2])).(statpi.var"Load Voltage at Bus 6"),
                 (x->(x[end]-x[2])).(dynpi.var"Load Voltage at Bus 6")], label=[raw"Static $\pi$ model", raw"Dynamic $\pi$ model"], bins=100)
    axs[3].hist([(x->(x[end]-x[3])).(statpi.var"Load Voltage at Bus 8"),
                 (x->(x[end]-x[4])).(dynpi.var"Load Voltage at Bus 8")], label=[raw"Static $\pi$ model", raw"Dynamic $\pi$ model"], bins=100)
    for ax in axs
        ax.set_yscale("log")
        ax.legend(prop=Dict("size"=>6))
    end
    axs[3].set_xlabel(raw"$V(t_f)-V(t_0)$ (pu)")
    axs[1].set_ylabel("Bus 5\n\n")
    axs[2].set_ylabel("Bus 6\n\n")
    axs[3].set_ylabel("Bus 8\n\n")
    fig.suptitle("Net Voltage Change During Transient Simulation (0.07s)")
    plt.show()


end
# transient()
