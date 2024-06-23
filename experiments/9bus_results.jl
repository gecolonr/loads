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
using Colors

# pygui(true)
include("sysbuilder.jl")
include("plotting.jl")
# include("9bus_sims.jl") # to make sure we have all the methods for deserialization

############## LOAD DATA ###############
gss = load_serde_data("data/fineresults_powersetpt")


######################################################
################### SCRATCH SPACE ####################
######################################################
function get_zipe_load_voltages_sm(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
    sm = small_signal_analysis(newsim)
    out = zeros(length(get_components(StandardLoad, gss.base)))
    for (idx, i) in enumerate(get_number.(get_bus.(get_components(StandardLoad, gss.base))))
        bus_lookup = PSID.get_bus_lookup(sim.results)
        bus_ix = bus_lookup[i]
        V_R = sm.operating_point[bus_ix]
        V_I = sm.operating_point[bus_ix+PSID.get_bus_count(sim.results)]
        out[idx] = sqrt(V_R^2 + V_I^2)
    end

    return out
end
add_result!(gss, ["Bus 3 Inverter Current", "Bus 1 Inverter Current", "Bus 2 Inverter Current"], get_injector_currents)
add_result!(gss, "posttripsm", small_signal_tripped)
add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
add_result!(gss, ["Final Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages_sm)

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
η_names = Dict(
    unique(getvec.(gss.df."ZIPE Load Params")) .=> ["η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η","η"]
)
inj_case_names = OrderedDict(
    "SM, GFL, GFM" => "Case 1",
    "SM, GFM, GFL" => "Case 2", 
    "GFM, GFL, SM" => "Case 3",
    "SM, SM, SM" => "Case 4",
    "GFM, GFM, SM" => "Case 5",
    "GFM, SM, SM" => "Case 6",
)

# first we add columns for all the data we want to include in the plot
getvec(x::LoadParams) = round.([x.z_percent, x.i_percent, x.p_percent, x.e_percent], digits=3)
df = @subset gss.df begin
    [i ∈ keys(η_names) for i in getvec.(:"ZIPE Load Params")]
    [i ∈ keys(inj_case_names) for i in join.(eachrow(hcat(:"injector at {Bus1}", 
                                                          :"injector at {Bus 2}",
                                                          :"injector at {Bus 3}")), ", ")]
end
df[!, "ZIPE Parameters"] = (x->"$(η_names[round.(x, digits=3)]): $(round.(x, digits=3))").(getvec.(df."ZIPE Load Params"))
df[!, "hovertext"] = (x->"η=$x").(getvec.(df."ZIPE Load Params"))
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ")
df[!, "tracename"] = map(x->inj_case_names[x[1]]*", "*x[2], zip(df[!, "Injector Setup"], df[!, "Line Model"]))
df[!, "Injector Setup"] = map(x->inj_case_names[x]*" ($x)", df[!, "Injector Setup"])
df[!, "real_eigs"] = real.(df[!, "Eigenvalues"])
df[!, "imag_eigs"] = imag.(df[!, "Eigenvalues"])
df[!, "op_pt"] = 
function pfactors(sm)
    pf = summary_participation_factors(sm)
    return Dict(pf.Name .=> eachrow(select(pf, Not(:Name))))
end
relu = x->((x+6.0)>0.0) ? (x+6.0)*5 : 0.0
df[!, "participation_factors"] = map(sm->get(pfactors(sm), "generator-3-1 ir_cnv", missing), df.sm)
df[!, "markersize"] = map(x->x isa Missing ? x : relu.(log10.(collect(x))), df.participation_factors)
df[!, "participation_factors"] = map(x->x isa Missing ? x : "PF: ".*string.(round.(relu.(log10.(collect(x))), sigdigits=3)), df.participation_factors)
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

    col_title_func=identity,
    row_title_func=identity,
    supertitle="Bus 3 Injector Current",

    yaxis_home_range = (min=0, max=10),
    xaxis_home_range = nothing,

    image_export_filename = "transient_current_plot",
    colorlist = Colors.distinguishable_colors(21),
)


savehtmlplot(p, "media/transient_currents_bus3_manycases")


p = makeplots(
    df;

    rows="Line Model",
    cols="Injector Setup",
    color="ZIPE Parameters",
    trace_names="Injector Setup",
    hovertext="participation_factors",

    slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x="real_eigs",
    y="imag_eigs",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Volage}\:\:[\mathrm{p.u.}]",

    row_title_func=identity, #"Line Model",
    col_title_func=identity, #"Injector Setup",
    supertitle="System Eigenvalues",

    xaxis_home_range = nothing,

    image_export_filename = "eigenvalue_plot",
    scatterplot_args = Dict(:mode=>"markers"),
)

savehtmlplot(p, "media/eigplot_manycase")


p = makeplots(
    df;
	
    rows="Line Model",
    cols="Injector Setup",
    color="ZIPE Parameters",
    slider="Power Setpoint",
    slider_current_value_prefix="Power Setpoint: ",
	
    x=(x0=0.48, dx=0.00005),
    y="Bus 3 Inverter Current",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Current}\:\:[\mathrm{p.u.}]",
	
    supertitle="Bus 3 Injector Current",
    yaxis_home_range = (min=0, max=10),
)
savehtmlplot(p, "media/transient_demo")




gss = load_serde_data("data/fineresults_powersetpt_nolines")
expand_columns!(gss)

df = @subset gss.df begin
    :z_percent .< 0.5
    [i ∈ keys(inj_case_names) for i in join.(eachrow(hcat(:"injector at {Bus1}", 
                                                          :"injector at {Bus 2}",
                                                          :"injector at {Bus 3}")), ", ")]
end
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ")
df[!, "Injector Setup"] = map(x->inj_case_names[x]*" ($x)", df[!, "Injector Setup"])
df[!, "Max Real Eigenvalue"] = [[i] for i in maximum.(real.(df[!, "Eigenvalues"]))]
df[!, "Power Setpoint 2"] = [[i] for i in df[!, "Power Setpoint"]]
df[!, "cpload"] = (x->x≈0.0 ? "No E Load" : "No P Load").(df[!, "e_percent"])
df[!, "zplusi"] = df[!, "z_percent"] .+ df[!, "i_percent"]
df[!, "ZIPE Parameters"] = "η=".*string.(getvec.(df[!, "ZIPE Load Params"]))
p = makeplots(
    df,
    rows="Line Model",
    cols="Injector Setup",
    color="ZIPE Parameters",
    x="Power Setpoint 2",
    y="Max Real Eigenvalue",
    # slider="cpload",
    # slider_label_func=identity,
    # slider_current_value_prefix="P or E load:",
    color_sort_func=x->parse(Float64, first(last(x, 9), 3))*0.9 + parse(Float64, first(last(x, 4), 3)),
    opacity=0.5,
    markersize=12,
    trace_names="Power Setpoint",
    hovertext="ZIPE Parameters",

    x_title="Power Setpoint",
    y_title="Maximum Real Eigenvalue Part",
    supertitle="Maximum Real Eigenvalue Part vs. Power Setpoint (Reds: P load, Blues: E load)",
    colorlist=collect(Iterators.flatten(zip(colormap("reds", 20; mid=0.5, logscale=false)[10:end], colormap("blues", 20; mid=0.5, logscale=false)[10:end])))

)
savehtmlplot(p, "media/maxeigvspwrsetpt_noslider_2")
df[!, "e_percent"] = [[i] for i in df[!, "e_percent"]]
df[!, "z_percent_nice"]="η_z=η_i=" .* string.(df[!, "z_percent"])
df[!, "Power Setpoint 3"] = "Load Scale: " .* string.(df[!, "Power Setpoint"])
p = makeplots(
    df,
    rows="Line Model",
    cols="Injector Setup",
    color="Power Setpoint 3",
    x="e_percent",
    y="Max Real Eigenvalue",
    # slider="Power Setpoint",
    slider_current_value_prefix="Power Setpoint: ",
    # hovertext="Power Setpoint",
    

    opacity=0.5,
    markersize=7,
    trace_names="ZIPE Parameters",
    hovertext="ZIPE Parameters",

    x_title="E Load Percent",
    y_title="Maximum Real Eigenvalue Part",
    supertitle="Maximum Real Eigenvalue Part vs. E/P Load Tradeoff",

    colorlist=collect(range(colorant"red", stop=colorant"blue", length=9))

)

savehtmlplot(p, "maxeigvsepercent")
colors = range(colorant"red", stop=colorant"blue", length=5)
colordict = OrderedDict(sort(unique(df[!, "z_percent"])) .=> colors)
df[!, Symbol("markershape")] = (x->x≈0.0 ? "square" : "circle").(df[!, "e_percent"])
df[!, "markercolor"] = [colordict[i] for i in df[!, "z_percent"]]
df[!, "Max Real Eigenvalue"]=[first(i) for i in df[!, "Max Real Eigenvalue"]]
p = PlotlyJS.plot(

    df,

    x=Symbol("Power Setpoint"),

    y=Symbol("Max Real Eigenvalue"),

    color=Symbol("ZIPE Parameters"),

    facet_col=Symbol("Injector Setup"),

    # facet_row=Symbol("Line Model"),

    # symbol=Symbol("cpload"),

    kind="scatter",

    # mode="markers",

    # marker=attr(

    #     size=7,

    #     symbol=Symbol("markershape"),
    #     # color=Symbol("markercolor")

    # ),
    # category_orders=attr(

    #     day=["Thur", "Fri", "Sat", "Sun"],

    #     smoker=["Yes", "No"],

    #     sex=["Male", "Female"]

    # )

)
