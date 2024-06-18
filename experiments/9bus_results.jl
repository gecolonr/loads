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
using OrderedCollections

# pygui(true)
include("sysbuilder.jl")
include("plotting.jl")
include("9bus_sims.jl") # to make sure we have all the methods for deserialization

############## LOAD DATA ###############
gss = load_serde_data("data/fineresults_powersetpt")


######################################################
################### SCRATCH SPACE ####################
######################################################

add_result!(gss, ["Bus 3 Inverter Current", "Bus 1 Inverter Current", "Bus 2 Inverter Current"], get_injector_currents)


η_names = OrderedDict(
    [0.1, 0.1, 0.1, 0.7]     => "High E Load",
    [0.0, 0.0, 0.0, 1.0]     => "Full E Load",
    [0.3, 0.3, 0.3, 0.1]     => "Low E Load",
    [0.2, 0.2, 0.2, 0.4]     => "Medium E Load",
    [0.15, 0.15, 0.15, 0.55] => "Low P High E Load",
    [0.0, 0.0, 1.0, 0.0]     => "Constant Power",
    [0.15, 0.15, 0.55, 0.15] => "High P Low E Load",
    [1.0, 0.0, 0.0, 0.0]     => "Constant Impedance",
)
inj_case_names = OrderedDict(
    "SM, GFL, GFM" => "Case 1",
    "SM, GFM, GFL" => "Case 2", 
    "GFM, GFL, SM" => "Case 3",
    "SM, SM, SM" => "Case 4"
)
df = gss.df
# first we add columns for all the data we want to include in the plot
getvec(x::LoadParams) = [x.z_percent, x.i_percent, x.p_percent, x.e_percent]
df[!, "ZIPE Parameters"] = (x->"$(η_names[round.(x, digits=3)]): $(round.(x, digits=3))").(getvec.(df."ZIPE Load Params"))
df[!, "hovertext"] = (x->"η=$x").(getvec.(df."ZIPE Load Params"))
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ")
df[!, "tracename"] = map(x->inj_case_names[x[1]]*", "*x[2], zip(df[!, "Injector Setup"], df[!, "Line Model"]))
df[!, "Injector Setup"] = map(x->inj_case_names[x]*" ($x)", df[!, "Injector Setup"])
df[!, "real_eigs"] = real.(df[!, "Eigenvalues"])
df[!, "imag_eigs"] = imag.(df[!, "Eigenvalues"])
function pfactors(sm)
    pf = summary_participation_factors(sm)
    return Dict(pf.Name .=> eachrow(select(pf, Not(:Name))))
end
df[!, "participation_factors"] = map(sm->get(pfactors(sm), "generator-3-1 ir_cnv", missing), df.sm)
df[!, "participation_factors"] = map(x->x isa Missing ? x : "PF: ".*string.(collect(x)), df.participation_factors)
# then we plot it!
p = makeplots(
    df;

    rows="Line Model",
    cols="Injector Setup",
    color="ZIPE Parameters",
    trace_names="tracename",
    hovertext="hovertext",

    slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x=(x0=0.48, dx=0.00005),
    y=["Bus 3 Inverter Current", "Load Voltage at Bus 5"][1],
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Current}\:\:[\mathrm{p.u.}]",

    col_titles="Injector Setup",
    row_titles="Line Model",
    supertitle="Bus 3 Injector Current",

    yaxis_home_range = (min=0, max=10),
    xaxis_home_range = nothing,

    image_export_filename = "transient_current_plot",
)


savehtmlplot(p, "media/transient_currents_bus3")


p = makeplots(
    df;

    rows="Line Model",
    cols="Injector Setup",
    color=:"ZIPE Parameters",
    trace_names="tracename",
    hovertext="participation_factors",

    slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x="real_eigs",
    y="imag_eigs",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Volage}\:\:[\mathrm{p.u.}]",

    col_titles="Injector Setup",
    row_titles="Line Model",
    supertitle="System Eigenvalues",

    # yaxis_home_range = (min=0, max=10),
    xaxis_home_range = nothing,

    image_export_filename = "eigenvalue_plot",
    scatterplot_args = Dict(:mode=>"markers"),
)

savehtmlplot(p, "media/eigplot_4case")