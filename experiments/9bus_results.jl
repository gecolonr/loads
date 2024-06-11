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

"""
all things to make functions here backwards compatible with earlier sims we've run, like adding power setpoint column if it's not present
"""
function set_power_setpt! end
function preprocess!(df::Union{DataFrame, GridSearchSys})

    if df isa GridSearchSys
        df = df.df
    end

    # expand ZIPE params and LineModelParams into individual columns
    expand_columns!(df)

    if !("Power Setpoint" ∈ names(df))
        df[!, Symbol("Power Setpoint")] .= 1.0
    end
    if "Error" ∈ names(df)
        rename!(df, "Error"=>"error")
    end



    return df
end

############## LOAD DATA ###############
## all cases from paper, with power variation
gss = load_serde_data("data/fineresults_powersetpt")


"""
utility to save `plot` to `filename`.html
"""
function savehtmlplot(plot, filename::String=nothing)
    # if we pass the output of Plot() it's the wrong type
    if plot isa PlotlyJS.SyncPlot
        plot = plot.plot
    end
    # default filename
    if isnothing(filename)
        try # use try/catch here because it could be nothing or it could be a dict with no :filename key or it could be ok
            filename = plot.config.toImageButtonOptions[:filename]
        catch
            filename = "plot.html"
        end
    end
    # ensure proper ending
    if !occursin(r"\.html$", filename)
         filename *= ".html"
    end
    # save to file
    open(filename, "w") do io
        PlotlyBase.to_html(io, plot)
    end
end


######################################################
################### SCRATCH SPACE ####################
######################################################

# function get_participation_factors(sim, bus_idx)
#     sm = small_signal_analysis(sim)
#     jacwrapper = get_jacobian(ResidualModel, sim.inputs, sm.operating_point, 0)
#     eigen_vals, eigen_vecs = PSID._get_eigenvalues(jacwrapper.Jv, true)
#     statemap = deepcopy(PSID.make_global_state_map(sim.inputs))
#     statemap["busses"] = Dict([Symbol(string(i))=>i for i in 1:19]...)
#     participation_factors = PSID._get_participation_factors(eigen_vecs, PSID._make_reduced_jacobian_index(statemap, [true for _ in 0:length(sim.results.solution.u[1])]))

#     # covariance matrix of participation factors shows high correlation between v_i and v_r at each bus.
#     # so we'll just take the average. (surely this is the participation factor towards the voltage magnitude)

#     ## code to get cov matrix
#     # p = participation_factors
#     # P = hcat([p[Symbol(string(i))].-(sum(p[Symbol(string(i))])/length(p[Symbol(string(i))])) for i in 1:19]...)
#     # for i in 1:19
#     #     P[:, i] /= sqrt(sum(P[:, i].^2))
#     # end
#     # cov_p = P'*P
#     return eigen_vals, (participation_factors["busses"][Symbol(string(bus_idx))] .+ participation_factors["busses"][Symbol(string(bus_idx+PSID.get_bus_count(sim.inputs)))]).*0.5
# end

# function add_eigs!(df)
#     x = hcat([[a, b] for (a, b) in get_participation_factors.(df.sim, 5)]...)
#     df.Eigenvalues .= x[1, :]
#     df[!, :"Participation Factors"] .= real.(x[2, :]) # type is complex but im component is always zero. just from hcat making types funny
# end


"""
Plot data from a dataframe. allows high dimensional data through a grid of subplots, trace colors, and a slider.

// TODO: Finish these docs

## Arguments
 - `df`: the DataFrame to get data from.
 - `rows`: the column of the dataframe specifying which row of subplots this data should be on
 - `cols`: the column of the dataframe specifying which column of subplots this data should be on
 - `color`: the column of the dataframe specifying which color this trace should be
 - `trace_names`: the column of the dataframe specifying the text name of this trace
 - `hovertext`: the column of the dataframe specifying any additional text that appears for this trace on hover
 - `slider`: the column of the dataframe specifying which slider value this trace corresponds to
"""
function makeplots(
    df::DataFrame;

    rows::Union{AbstractString, Symbol, Nothing}=nothing,
    cols::Union{AbstractString, Symbol, Nothing}=nothing,
    color::Union{AbstractString, Symbol, Nothing}=nothing,
    trace_names::Union{AbstractString, Symbol, Nothing}=nothing,
    hovertext::Union{AbstractString, Symbol, Nothing}=nothing,

    slider::Union{AbstractString, Symbol, Nothing}=nothing,
    slider_sort_func::Function = identity,
    slider_label_func::Function = x->round(x, sigdigits=2),
    slider_current_value_prefix::Union{AbstractString, Nothing} = nothing,

    x::Union{AbstractString, Symbol, NamedTuple{(:x0, :dx), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    y::Union{AbstractString, Symbol, NamedTuple{(:y0, :dy), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    x_title::Union{AbstractString, Nothing} = nothing,
    y_title::Union{AbstractString, Nothing} = nothing,

    col_titles::Union{AbstractString, Symbol, Vector{Union{AbstractString}}, Nothing}=nothing,
    row_titles::Union{AbstractString, Symbol, Vector{Union{AbstractString}}, Nothing}=nothing,
    supertitle::Union{AbstractString, Nothing} = nothing,

    yaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,
    xaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,

    image_export_filename::String = "saved_plot",
    colorlist::Vector{String} = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
    scatterplot_args::Dict = Dict(),
)

    if !isnothing(color) && length(unique(df[!, color]))>length(colorlist)
        AssertionError("Not enough colors to give each series a unique color. Need at least $(length(unique(df[!, color]))).")
    end

    spkwargs = Dict()
    if !isnothing(col_titles) spkwargs[:column_titles] = col_titles isa Vector ? col_titles : unique(df[!, col_titles]) end
    if !isnothing(row_titles) spkwargs[:row_titles] = row_titles isa Vector ? row_titles : unique(df[!, row_titles]) end
    if !isnothing(x_title) spkwargs[:x_title] = x_title end
    if !isnothing(y_title) spkwargs[:y_title] = y_title end

    sp = Subplots(
        rows=isnothing(rows) ? 1 : length(unique(df[!, rows])), 
        cols=isnothing(cols) ? 1 : length(unique(df[!, cols])),
        shared_xaxes="all", shared_yaxes="all",
        start_cell = "top-left",
        horizontal_spacing=0.03,
        vertical_spacing=0.03;
        spkwargs...
    )

    layout_kwargs = Dict()
    if !isnothing(supertitle) layout_kwargs[:title_text] = supertitle end
    if !isnothing(xaxis_home_range) layout_kwargs[:xaxis_range] = [xaxis_home_range.min, xaxis_home_range.max] end
    if !isnothing(yaxis_home_range) layout_kwargs[:yaxis_range] = [yaxis_home_range.min, yaxis_home_range.max] end

    slider_values = isnothing(slider) ? [0.0] : sort(unique(df[!, slider]), by=slider_sort_func)
    layout_kwargs[:sliders] = [attr(
        steps=[
            attr(
                label=string(slider_label_func(slider_values[i])),
                method="restyle",
                # args=["visible", [j%2==i%2 for j in 1:48]]
                # args = [Dict(:visible=>(i%2==0)) for i in 1:48],
                # args=[attr(z=(zz[i, :, :],), text=(text[i, :, :],)), 0]
            )
            for i in 1:length(slider_values)
        ],
        active=1,
        currentvalue_prefix=slider_current_value_prefix,
        pad_t=40,
    )]

    l = Layout(sp; layout_kwargs...)

    fig = PlotlyJS.plot(l, config = PlotConfig(
        scrollZoom=true,
        modeBarButtonsToAdd=[
            "drawline",
            "drawcircle",
            "drawrect",
            "eraseshape"
        ],
        toImageButtonOptions=attr(
            format="svg", # one of png, svg, jpeg, webp
            filename=image_export_filename,
            height=900,
            width=1600,
            scale=1
        ).fields,
        displayModeBar=true,
    ))

    layoutkeys = keys(fig.plot.layout)
    layoutkeys = collect(layoutkeys)[(x->(occursin("axis", string(x)))).(layoutkeys)]
    for ax in layoutkeys
        fig.plot.layout[ax][:showticklabels] = true
    end
    fig.plot.layout.margin[:l] *= 2
    fig.plot.layout.margin[:b] *= 2
    fig.plot.layout.margin[:r] *= 2
    fig.plot.layout.margin[:t] *= 2
    
    # [(name=colorname, data=dict(sliderval=>data for each data of color)) for (colorname, color) in colors]
    sliderdata = []
    for (rowidx, rowval) in enumerate(isnothing(rows) ? [1] : unique(df[!, rows]))
        for (colidx, colval) in enumerate(isnothing(cols) ? [1] : unique(df[!, cols]))
            for (coloridx, colorval) in enumerate(isnothing(color) ? [1] : unique(df[!, color]))
                # data = @subset(df, Symbol(rows) .== rowval, Symbol(cols) .== colval, Symbol(color) .== colorval)
                data = subset(df, [i=>ByRow(x->x==ival) for (i, ival) in [(rows, rowval), (cols, colval), (color, colorval)] if !isnothing(i)])
                # println([i for i in [x, y, hovertext] if i isa Union{AbstractString, Symbol}])
                push!(sliderdata, Dict((isnothing(slider) ? [i for i in 1:length(eachrow(data))] : data[!, slider]) .=> eachrow(select(data, [i for i in [x, y, hovertext, trace_names] if i isa Union{AbstractString, Symbol}]))))

                scatargs = Dict()
                if (x isa Vector) 
                    scatargs[:x] = x 
                elseif (x isa NamedTuple) 
                    scatargs[:x0] = x.x0
                    scatargs[:dx] = x.dx
                # elseif isnothing(slider)
                #     scatargs[:x] = data[!, x]
                end # if it's a column, we'll update it, so don't put it there now unless we aren't doing slider stuff

                if (y isa Vector) 
                    scatargs[:y] = y
                elseif (y isa NamedTuple)
                    scatargs[:y0] = y.y0
                    scatargs[:dy] = y.dy
                # elseif isnothing(slider)
                #     scatargs[:y] = data[!, y]
                end
                
                if !isnothing(trace_names) scatargs[:name]=data[!, trace_names][1] end
                if !isnothing(color)
                    scatargs[:legendgroup]=coloridx
                    scatargs[:legendgrouptitle_text]=color
                end
                if !isnothing(hovertext) scatargs[:hovertext]=data[!, hovertext][1] end

                
                scatargs[:opacity]=0.7
                scatargs[:marker]=isnothing(color) ? attr(size=7) : attr(size=7, color=colorlist[coloridx])

                if isnothing(slider)
                    for i in values(sliderdata[end])
                        if x isa Union{AbstractString, Symbol} scatargs[:x]=i[x] end
                        if y isa Union{AbstractString, Symbol} scatargs[:y]=i[y] end
                        # println(typeof(scatargs[:y]))
                        if !isnothing(trace_names) scatargs[:name]=i[trace_names] end
                        if !isnothing(hovertext) scatargs[:hovertext]=i[hovertext] end
                        PlotlyJS.add_trace!(fig, PlotlyJS.scatter(; merge(scatargs, scatterplot_args)...), row=rowidx, col=colidx)
                    end
                else
                    PlotlyJS.add_trace!(fig, PlotlyJS.scatter(; merge(scatargs, scatterplot_args)...), row=rowidx, col=colidx)
                end
            end
        end
    end
    if isnothing(slider)
        return fig
    end
    for (step, sliderval) in zip(fig.plot.layout[:sliders][1].steps, sort(unique(df[!, slider]), by=slider_sort_func))
        step.args = [Dict()]
        if x isa Union{AbstractString, Symbol}
            step.args[1]["x"]=[data[sliderval][x] for data in sliderdata]
        end
        if y isa Union{AbstractString, Symbol}
            step.args[1]["y"]=[data[sliderval][y] for data in sliderdata]
        end
        if !isnothing(trace_names) step.args[1]["name"]      = [data[sliderval][trace_names] for data in sliderdata] end
        if !isnothing(hovertext)   step.args[1]["hovertext"] = [data[sliderval][hovertext]   for data in sliderdata] end
    end
    return fig
end

add_result!(gss, ["Bus 3 Inverter Current", "Bus 1 Inverter Current", "Bus 2 Inverter Current"], get_inverter_currents)


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
# then we plot it!
p = makeplots(
    df;

    rows="Line Model",
    cols="Injector Setup",
    color=:"ZIPE Parameters",
    trace_names="tracename",
    hovertext="hovertext",

    slider="Power Setpoint",
    slider_sort_func = identity,
    slider_label_func = x->round(x, digits=3),
    slider_current_value_prefix="Power Setpoint: ",

    x=(x0=0.48, dx=0.00005),
    y=["Bus 3 Inverter Current", "Load Voltage at Bus 5"][2],
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Volage}\:\:[\mathrm{p.u.}]",

    col_titles="Injector Setup",
    row_titles="Line Model",
    supertitle="Load Voltage at Bus 5",

    yaxis_home_range = (min=0, max=10),
    xaxis_home_range = nothing,

    image_export_filename = "transient_voltage_plot",
)


savehtmlplot(p, "media/transient_voltage_bus5")


p = makeplots(
    df;

    rows="Line Model",
    cols="Injector Setup",
    color=:"ZIPE Parameters",
    trace_names="tracename",
    hovertext="hovertext",

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