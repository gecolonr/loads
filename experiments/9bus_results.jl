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
using PyCall
using Logging
Logging.disable_logging(Logging.Error)
# using Pandas
using CategoricalArrays
# px = pyimport("plotly.express")
# pd = pyimport("pandas")
# pygui(true)
include("sysbuilder.jl")
include("plotting.jl")
# include("9bus_sims.jl") # to make sure we have all the methods for deserialization

############## LOAD DATA ###############
gss = load_serde_data("data/forplot")
gss = load_serde_data("data/fineresults_powersetpt")
# gss = load_serde_data("/data/reiddye/loads/fivesecondsbetter")
map!(getvec, gss.df[!, "ZIPE Load Params"], gss.df[!, "ZIPE Load Params"])
# gss.df[!, "injector at {Bus1}"] = categorical(gss.df[!, "injector at {Bus1}"])
# gss.df[!, "injector at {Bus 2}"] = categorical(gss.df[!, "injector at {Bus 2}"])
# gss.df[!, "injector at {Bus 3}"] = categorical(gss.df[!, "injector at {Bus 3}"])
# gss.df[!, "Power Setpoint"] = categorical(gss.df[!, "Power Setpoint"])
# gss.df[!, "Line Model"] = categorical(gss.df[!, "Line Model"])

# df = @subset gss.df begin
    # unwrap.(:"Power Setpoint") .< 1.01
# end

######################################################
################### SCRATCH SPACE ####################
######################################################
function get_final_injector_currents(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    # Base.show(sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    # Base.show(sys)
    
    (newsim, newsm, dt, newerror) = runSim(sys, PerturbState(0.5, 1, 0.0), ResidualModel, (0.4, 0.6), IDA(), 0.01)
    # println("<here>")
    # Base.show(read_results(newsim))
    # Base.show(newsm)
    # Base.show(newerror)
    # println("</here>")
    if isnothing(read_results(newsim))
        println("ERROR: $(newerror)")
        return [missing, missing, missing]
    end
    return map(x->ismissing(x) ? x : first(x), get_injector_currents(gss, newsim, newsm, newerror))
end
add_result!(gss, ["Bus 3 Injector Current", "Bus 1 Injector Current", "Bus 2 Injector Current"], get_injector_currents)
# add_result!(gss, ["Bus 3 Inverter Current Final", "Bus 1 Inverter Current Final", "Bus 2 Inverter Current Final"], get_final_injector_currents)
# add_result!(gss, "posttripsm", small_signal_tripped)
# add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
# add_result!(gss, ["Final Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages_sm)

η_names = OrderedDict(
    [0.1, 0.1, 0.1, 0.7]     => "High E Load",
    [0.0, 0.0, 0.0, 1.0]     => "Full E Load",
    [0.3, 0.3, 0.3, 0.1]     => "Low E Load",
    [0.2, 0.2, 0.2, 0.4]     => "Medium E Load",
    # [0.15, 0.15, 0.15, 0.55] => "Low P High E Load",
    # [0.0, 0.0, 1.0, 0.0]     => "Constant Power",
    # [0.15, 0.15, 0.55, 0.15] => "High P Low E Load",
    # [1.0, 0.0, 0.0, 0.0]     => "Constant Impedance",
)
getvec(x::LoadParams) = round.([x.z_percent, x.i_percent, x.p_percent, x.e_percent], digits=3)
η_names = Dict(
    [[0.2, 0.2, 0.0, 0.6],
    #  [0.35, 0.35, 0.0, 0.3],
     [0.3, 0.3, 0.0, 0.4],
    #  [0.15, 0.15, 0.7, 0.0],
     [0.3, 0.3, 0.4, 0.0],
    #  [0.05, 0.05, 0.9, 0.0],
    #  [0.25, 0.25, 0.0, 0.5],
     [0.1, 0.1, 0.0, 0.8],
    #  [0.45, 0.45, 0.0, 0.1],
    #  [0.25, 0.25, 0.5, 0.0],
     [0.4, 0.4, 0.2, 0.0],
     [0.1, 0.1, 0.8, 0.0],
     [0.5, 0.5, 0.0, 0.0],
    #  [0.45, 0.45, 0.1, 0.0],
     [0.2, 0.2, 0.6, 0.0],
     [0.4, 0.4, 0.0, 0.2],
    #  [0.05, 0.05, 0.0, 0.9],
    #  [0.35, 0.35, 0.3, 0.0],
    #  [0.15, 0.15, 0.0, 0.7],
     [0.0, 0.0, 0.0, 1.0],
     [0.0, 0.0, 1.0, 0.0],
     [1.0, 0.0, 0.0, 0.0]] .=> ["η","η","η","η","η","η","η","η","η","η","η","η",]#"η","η","η","η","η","η","η","η","η","η"]
)
inj_case_names = OrderedDict(
    "SM, GFL, GFM" => "Case 1",
    "SM, GFM, GFL" => "Case 2", 
    "SM, GFM, GFL" => "Case 2", 
    "GFM, GFL, SM" => "Case 3",
    "SM, SM, SM" => "Case 4",
    "GFM, GFM, SM" => "Case 5",
    "GFM, SM, SM" => "Case 6",
)

# first we add columns for all the data we want to include in the plot
df = @subset gss.df begin
#     :"Power Setpoint" .≈ 0.2

#     # @byrow :"Line Model" ∈ ["statpi (dommel)", "dynpi (dommel)"]
    @byrow getvec(:"ZIPE Load Params") ∈ keys(η_names)
    [i ∈ keys(inj_case_names) for i in join.(eachrow(hcat(:"injector at {Bus1}", 
                                                          :"injector at {Bus 2}",
                                                          :"injector at {Bus 3}")), ", ")]
end

df = @subset gss.df begin
    # :"Power Setpoint" .≈ 0.2
end
df[!, "ZIPE Parameters"] = (x->"$(η_names[round.(x, digits=3)]): $(round.(x, digits=3))").(getvec.(df."ZIPE Load Params"))
df[!, "hovertext"] = (x->"η=$x").(getvec.(df."ZIPE Load Params"))
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ")
# df[!, "tracename"] = map(x->inj_case_names[x[1]]*", "*x[2], zip(df[!, "Injector Setup"], unwrap.(df[!, "Line Model"])))
df[!, "Injector Setup"] = map(x->inj_case_names[x]*" ($x)", df[!, "Injector Setup"])
df[!, "real_eigs"] = real.(df[!, "Eigenvalues"])
df[!, "imag_eigs"] = imag.(df[!, "Eigenvalues"])


# df_plotable=to_plottable!(select(df, 
# ["Power Setpoint", "ZIPE Parameters", "Line Model", "real_eigs", "imag_eigs"]), 
# ["real_eigs", "imag_eigs", "Load Voltage at Bus 5"])

# # select(df, ["Power Setpoint", "ZIPE Parameters", "Injector Setup", "Line Model", "real_eigs", "imag_eigs"])

# p = PlotlyJS.plot(
#     subset(select(df_plotable, ["Power Setpoint", "ZIPE Parameters", "Injector Setup", "Line Model", "real_eigs", "imag_eigs"]), ["real_eigs"=>ByRow(x->!isnothing(x))]), 
#     x=:real_eigs,
#     y=:imag_eigs,
#     color=Symbol("ZIPE Parameters"),
#     facet_col=Symbol("Injector Setup"),
#     facet_row=Symbol("Line Model"),
#     animation_frame=Symbol("Power Setpoint"),
#     mode="markers",
#     hover_name=Symbol("hovertext"),
# )
# df_plotable[!, "Injector Setup"] .= categorical(df_plotable[!, "Injector Setup"])
# df_plotable[!, "Line Model"] .= categorical(df_plotable[!, "Line Model"])
# df_plotable[!, "Power Setpoint"] .= categorical(df_plotable[!, "Power Setpoint"])
# df_plotable[!, "ZIPE Parameters"] .= categorical(df_plotable[!, "ZIPE Parameters"])
# df_plotable[!, "real_eigs"] .= Vector{Float64}(map(x->isnothing(x) ? missing : x, df_plotable[!, "real_eigs"]))
# df_jlpandas = Pandas.DataFrame(select(df_plotable, Not(["real_eigs", "imag_eigs"])))
# df_pd = pd.DataFrame(df_jlpandas)
# fig = px.scatter(df_plotable, x="real_eigs", y="imag_eigs", animation_frame="Power Setpoint", #animation_group="country",
#            color="ZIPE Parameters", hover_name="hovertext")
# fig["layout"].pop("updatemenus") # optional, drop animation buttons
# fig.show()


function pfactors(sm)
    pf = summary_participation_factors(sm)
    return Dict(pf.Name .=> eachrow(select(pf, Not(:Name))))
end
# relu = x->((x+6.0)>0.0) ? (x+6.0)*5 : 0.0
df[!, "participation_factors"] = map(sm->get(pfactors(sm), "generator-3-1 ir_cnv", missing), df.sm)
# df[!, "markersize"] = map(x->x isa Missing ? x : relu.(log10.(collect(x))), df.participation_factors)
df[!, "participation_factors"] = map(x->x isa Missing ? "" : "\nPF: ".*string.(round.((collect(x)), sigdigits=5)), df.participation_factors)
# then we plot it!
df[!, "hovertext_eig"] = [i.*j for (i, j) in zip(df[!, "hovertext"], df[!, "participation_factors"])]
p = makeplots(
    # df;
    @subset df :"Power Setpoint" .≈ 0.2;

    rows="Line Model",
    # cols="Injector Setup",
    color="ZIPE Parameters",
    trace_names="Line Model",
    hovertext="hovertext",

    # slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x=[(x0=0.48, dx=0.00005), "time"][1],
    y=["Bus 3 Injector Current", "Load Voltage at Bus 5"][1],
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"||i_\text{filter}||\:\:[\mathrm{p.u.}]",

    col_title_func=identity,
    row_title_func=identity,
    # color_sort_func=x->-Dict(
    #     "η: [1.0, 0.0, 0.0, 0.0]" => 1,
    #     "η: [0.0, 0.0, 1.0, 0.0]" => 2,
    #     "η: [0.1, 0.1, 0.8, 0.0]" => 3,
    #     "η: [0.2, 0.2, 0.6, 0.0]" => 4,
    #     "η: [0.3, 0.3, 0.4, 0.0]" => 5,
    #     "η: [0.4, 0.4, 0.2, 0.0]" => 6,
    #     "η: [0.5, 0.5, 0.0, 0.0]" => 7,
    #     "η: [0.4, 0.4, 0.0, 0.2]" => 8,
    #     "η: [0.3, 0.3, 0.0, 0.4]" => 9,
    #     "η: [0.2, 0.2, 0.0, 0.6]" => 10,
    #     "η: [0.1, 0.1, 0.0, 0.8]" => 11,
    #     "η: [0.0, 0.0, 0.0, 1.0]" => 12,
    #  )[x],
    color_sort_func=x->last(x, 4),
    supertitle=L"\text{Bus 3 (GFL) Current Magnitude}",

    yaxis_home_range = (min=0.317, max=0.365),
    xaxis_home_range = (min=0.498, max=0.53),

    image_export_filename = "transient_current_plot",
    # colorlist = reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))))
    colorlist = Colors.distinguishable_colors(21),
)


savehtmlplot(p, "media/transient_currents_bus3_final_fourcase_slider")


p = makeplots(
    df;
    # @subset df :"Power Setpoint" .≈ 0.2;

    rows="Line Model",
    # cols="Injector Setup",
    color="ZIPE Parameters",
    trace_names="Line Model",
    hovertext="hovertext_eig",

    slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x="real_eigs",
    y="imag_eigs",
    x_title=L"\mathrm{Re}(\lambda)",
    y_title = L"\mathrm{Im}(\lambda)",

    col_title_func=identity,
    row_title_func=identity,
    # color_sort_func=x->-Dict(
    #     "η: [1.0, 0.0, 0.0, 0.0]" => 1,
    #     "η: [0.0, 0.0, 1.0, 0.0]" => 2,
    #     "η: [0.1, 0.1, 0.8, 0.0]" => 3,
    #     "η: [0.2, 0.2, 0.6, 0.0]" => 4,
    #     "η: [0.3, 0.3, 0.4, 0.0]" => 5,
    #     "η: [0.4, 0.4, 0.2, 0.0]" => 6,
    #     "η: [0.5, 0.5, 0.0, 0.0]" => 7,
    #     "η: [0.4, 0.4, 0.0, 0.2]" => 8,
    #     "η: [0.3, 0.3, 0.0, 0.4]" => 9,
    #     "η: [0.2, 0.2, 0.0, 0.6]" => 10,
    #     "η: [0.1, 0.1, 0.0, 0.8]" => 11,
    #     "η: [0.0, 0.0, 0.0, 1.0]" => 12,
    #  )[x],
    color_sort_func=x->last(x, 4),
    supertitle=L"\text{System Eigenvalues}",

    # yaxis_home_range = (min=0.317, max=0.358),
    # xaxis_home_range = (min=0.496, max=0.53),

    image_export_filename = "eigenvalue_plot",
    # colorlist = reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11)))),
    scatterplot_args = Dict(:mode=>"markers"),

    colorlist = Colors.distinguishable_colors(21),
)


savehtmlplot(p, "media/eigplot_final_fourcase_slider")

p = makeplots(
    df;

    rows="Line Model",
    # cols="Injector Setup",
    color="ZIPE Parameters",
    trace_names="Injector Setup",
    # hovertext="participation_factors",

    # slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x="real_eigs",
    y="imag_eigs",
    x_title=L"\mathrm{Re}(\lambda)",
    y_title = L"\mathrm{Im}(\lambda)",

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