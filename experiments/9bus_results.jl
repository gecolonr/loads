cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")
using PowerSystems
const PSY = PowerSystems
using ZIPE_loads
using InfrastructureSystems
# using PlotlyJS
using DataFrames
using TLmodels
using DataFramesMeta
using LaTeXStrings
using PyPlot
const plt = PyPlot

# using PyCall
# @pyimport mpld3 # if this doesn't work do `pip install mpld3` with the right pip for your pycall
using Plots
plotlyjs()
pygui(true)
include("sysbuilder.jl")

# set_units_base_system!

# LOAD DATA
df = load_serde_data("data/fineresults")
df2 = load_serde_data("data/results")
expand_columns!(df)
expand_columns!(df2)
df = df3
############### CASES ###################
case1 = @subset df begin
    :"injector at {Bus1}" .== "SM"
    :"injector at {Bus 2}" .== "GFL"
    :"injector at {Bus 3}" .== "GFM"
end

case2 = @subset df begin
    :"injector at {Bus1}" .== "SM"
    :"injector at {Bus 2}" .== "GFM"
    :"injector at {Bus 3}" .== "GFL"
end

case3 = @subset df begin
    :"injector at {Bus1}" .== "GFM"
    :"injector at {Bus 2}" .== "GFL"
    :"injector at {Bus 3}" .== "SM"
end

cases = Dict(
    "SM, GFL, GFM" => case1,
    "SM, GFM, GFL" => case2, 
    "GFM, GFL, SM" => case3,
)

η_combos = Dict(
    "Constant Impedance" => [1.0, 0.0, 0.0, 0.0],
    "Constant Power"     => [0.0, 0.0, 1.0, 0.0],
    "Full E Load"        => [0.0, 0.0, 0.0, 1.0],
    "High E Load"        => [0.1, 0.1, 0.1, 0.7],
    "Medium E Load"      => [0.2, 0.2, 0.2, 0.4],
    "Low E Load"         => [0.3, 0.3, 0.3, 0.1],
    "Low P High E Load"  => [0.15, 0.15, 0.15, 0.55],
    "High P Low E Load"   => [0.15, 0.15, 0.55, 0.15],
)


################## PLOTTING FUNCTIONS #################
# Eigenvalue Plot
function eigplot(df)
    geteigs(df) = [(name=name, eigs=@subset(df, :z_percent.≈η[1], :i_percent.≈η[2], :p_percent.≈η[3]).Eigenvalues) for (name, η) in (η_combos)]

    stateigs = [(name=casename, data=geteigs(@subset(casedf, :"Line Model" .== "statpi"))) for (casename, casedf) in cases]
    dyneigs = [(name=casename, data=geteigs(@subset(casedf, :"Line Model" .== "dynpi"))) for (casename, casedf) in cases]
    
    fig, axs = plt.subplots(2, 3, sharex=true, sharey=true)
    # axs[1].set_title("Static Pi Model")
    # axs[2].set_title("Dynamic Pi Model")
    fig.tight_layout()
    # fig, ax = plt.subplots(4, 4, sharex="col", sharey="row", figsize=(8, 8))
    # fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
    #                 hspace=0.1, wspace=0.1)
    reals = []
    imags = []
    labels = []
    for η_combo_idx in 1:length(η_combos)
        for inj_case_idx in 1:length(cases)
            println(stateigs[inj_case_idx].name, stateigs[inj_case_idx].data[η_combo_idx].name, length(stateigs[inj_case_idx].data[η_combo_idx].eigs))
            if stateigs[inj_case_idx].data[η_combo_idx].eigs isa Missing
                println("No data for statpi, case $inj_case_idx ($(stateigs[inj_case_idx].name)), $(stateigs[inj_case_idx].data[η_combo_idx].name)")
            else
                # scatter!(
                #     real.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)), 
                #     imag.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)),
                #     lab=string(stateigs[inj_case_idx].data[η_combo_idx].name),
                # )
                push!(reals, real.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)))
                push!(imags, imag.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)))
                push!(labels, string(stateigs[inj_case_idx].data[η_combo_idx].name))

                axs[2*inj_case_idx-1].scatter(
                    real.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)), 
                    imag.(reduce(vcat, stateigs[inj_case_idx].data[η_combo_idx].eigs)),
                    label=string(stateigs[inj_case_idx].data[η_combo_idx].name),
                    s=50,
                    alpha=0.5
                )
            end

            if dyneigs[inj_case_idx].data[η_combo_idx].eigs isa Missing
                println("No data for dynpi, case $inj_case_idx ($(dyneigs[inj_case_idx].name)), $(dyneigs[inj_case_idx].data[η_combo_idx].name)")
            else
                push!(reals, real.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)))
                push!(imags, imag.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)))
                push!(labels, string(dyneigs[inj_case_idx].data[η_combo_idx].name))
                axs[2*inj_case_idx].scatter(
                    real.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)), 
                    imag.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)),
                    label=string(dyneigs[inj_case_idx].data[η_combo_idx].name),
                    s=50,
                    alpha=0.5
                )
                # scatter!(
                #     real.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)), 
                #     imag.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs)),
                #     lab=string(dyneigs[inj_case_idx].data[η_combo_idx].name),
                # )
            end
        end
    end
    println(labels)
    scat = Plots.scatter(reals, imags, layout=(2, 3), label=vec(labels), markeralpha=0.75, size=(1500, 1200))
    axs[1].set_title("Case 1 ($(collect(keys(cases))[1]))")
    axs[3].set_title("Case 2 ($(collect(keys(cases))[2]))")
    axs[5].set_title("Case 3 ($(collect(keys(cases))[3]))")
    axs[2].set_xlabel(L"\mathrm{Re}(\lambda)")
    axs[4].set_xlabel(L"\mathrm{Re}(\lambda)")
    axs[6].set_xlabel(L"\mathrm{Re}(\lambda)")
    axs[1].set_ylabel("Statpi Model\n\n\n\n"*L"\mathrm{Im}(\lambda)")
    axs[2].set_ylabel("Dynpi Model\n\n\n\n"*L"\mathrm{Im}(\lambda)")
    fig.suptitle("System Eigenvalues")
    map(x->x.legend(), axs)
    # map(x->x.set_xscale("log"))
    map(x->x.axvline(x=0, color="black", ls="--", lw=1), axs)
    plt.show()
    return scat
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
    selectZIPE(df) = vcat([@subset(df, :z_percent.≈η[1], :i_percent.≈η[2], :p_percent.≈η[3]) for (name, η) in (η_combos)]...)

    statpidata = (@subset(selectZIPE(df), :"Line Model" .== "statpi"))
    dynpidata = (@subset(selectZIPE(df), :"Line Model" .== "dynpi"))

    # busses = ["Load Voltage at $i" for i in get_name.(get_components(Bus, df.sim[1].sys))]
    busses = ["Load Voltage at $i" for i in get_name.(get_bus.(get_components(StandardLoad, first(df.sim).sys)))]
    fig, axs = plt.subplots(2, length(busses), sharex=true, sharey=true)#, figsize=(25, 15))
    fig.tight_layout()
    fig.suptitle("$(first(busses)) for BranchTrip on line Bus 5-Bus 4-i₁")
    axs[1].set_ylabel("dynpi")
    axs[2].set_ylabel("statpi")
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "b", "g", "r", "c", "m", "y", "k", "w"]
    colordict = Dict()
    # axs = [Plot() Plot() Plot(); Plot() Plot() Plot()]
    
    # for (j, bus) in enumerate(busses)
    bus = first(busses)
    for (j, case) in enumerate(cases)
        statpidata = (@subset(selectZIPE(last(case)), :"Line Model" .== "statpi"))
        dynpidata = (@subset(selectZIPE(last(case)), :"Line Model" .== "dynpi"))
        for row in eachrow(dynpidata)
            η = round.((row.z_percent, row.i_percent, row.p_percent, row.e_percent), digits=3)
            if !(η in keys(colordict))
                colordict[η] = colors[length(colordict)+1]
            end

            if row[bus] isa Missing
                println("dynpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                axs[1, j].plot([0.48], [0.97], color=colordict[η], label="Z=$(η[1]), I=$(η[2]), P=$(η[3]), E=$(η[4]) [FAIL]")
                continue
            end

            # addtraces!(axs[1, j], PlotlyJS.scatter(; x=0.48:0.0001:0.55, y=row[bus]))
            η = round.(η, digits=3)
            axs[1, j].plot((0.48:0.00005:1.0)[begin : length(row[bus])], row[bus], color=colordict[η], alpha=0.75, label="Z=$(η[1]), I=$(η[2]), P=$(η[3]), E=$(η[4])")
        end
        for row in eachrow(statpidata)
            η = round.((row.z_percent, row.i_percent, row.p_percent, row.e_percent), digits=3)
            if !(η in keys(colordict))
                colordict[η] = colors[length(colordict)+1]
            end

            if row[bus] isa Missing
                println("statpi failure: Z=$(row.z_percent), I=$(row.i_percent), P=$(row.p_percent), E=$(row.e_percent)")
                axs[2, j].plot([0.48], [0.97], color=colordict[η], label="Z=$(η[1]), I=$(η[2]), P=$(η[3]), E=$(η[4]) [FAIL]")
                continue
            end

            # addtraces!(axs[2, j], PlotlyJS.scatter(; x=0.48:0.0001:0.55, y=row[bus]))
            axs[2, j].plot(collect(0.48:0.00005:1.0)[begin : length(row[bus])], row[bus], color=colordict[η], alpha=0.75, label="Z=$(η[1]), I=$(η[2]), P=$(η[3]), E=$(η[4])")
        end
        axs[1, j].set_ylim([0, 10])
        axs[2, j].set_ylim([0, 10])
        axs[1, j].set_title("Case $j ($(first(case)))")
        axs[2, j].set_xlabel("Time (s)")
        
        axs[1, j].legend()#prop=Dict("size"=>5))
        axs[2, j].legend()#prop=Dict("size"=>5))
    end
    plt.xlabel("Time (s)")
    # axs[1].set_ylabel("Voltage")
    # axs[2].set_ylabel("Voltage")
    plt.show()
    # plt.tight_layout()
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
