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

# pygui(true)
include("sysbuilder.jl")

"""
all things to make functions here backwards compatible with earlier sims we've run, like adding power setpoint column if it's not present
"""
function preprocess!(df::DataFrame)

    # expand ZIPE params and LineModelParams into individual columns
    expand_columns!(df)

    if !("Power Setpoint" in names(df))
        df[!, Symbol("Power Setpoint")] .= 1.0
    end

    return df
end

############## LOAD DATA ###############
## dataset with all cases from paper, no power variation
# df = preprocess!(load_serde_data("data/fineresults_smallertime"))

## all cases from paper, with power variation
df = preprocess!(load_serde_data("data/fineresults_powersetpt"))

## massive dataset with all combos of injectors and zipe gridsearch
# df = preprocess!(load_serde_data("data/results"))


############### CASES ###################
case1(df::DataFrame) = @subset df begin
    :"injector at {Bus1}" .== "SM"
    :"injector at {Bus 2}" .== "GFL"
    :"injector at {Bus 3}" .== "GFM"
end

case2(df::DataFrame) = @subset df begin
    :"injector at {Bus1}" .== "SM"
    :"injector at {Bus 2}" .== "GFM"
    :"injector at {Bus 3}" .== "GFL"
end

case3(df::DataFrame) = @subset df begin
    :"injector at {Bus1}" .== "GFM"
    :"injector at {Bus 2}" .== "GFL"
    :"injector at {Bus 3}" .== "SM"
end

cases(df::DataFrame) = Dict(
    "SM, GFL, GFM" => case1(df),
    "SM, GFM, GFL" => case2(df), 
    "GFM, GFL, SM" => case3(df),
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

# list of colors to use in order
# so it's consistent between all plots
colorlist = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # safety orange
    "#2ca02c",  # cooked asparagus green
    "#d62728",  # brick red
    "#9467bd",  # muted purple
    "#8c564b",  # chestnut brown
    "#e377c2",  # raspberry yogurt pink
    "#7f7f7f",  # middle gray
    "#bcbd22",  # curry yellow-green
    "#17becf"   # blue-teal
]

"""
makes fig with rows [statpi, dynpi] and columns [case1, case2, case3].

just a utility for these two plots since I don't want to copy/paste it.

`x_title` and `y_title` are the axis titles for the subplots. `filename` is the filename the plot saves as when you click the "download plot" button.
"""
function makefig(df::DataFrame, x_title::Union{String, LaTeXString}, y_title::Union{String, LaTeXString}, filename::String="plot", title_text::Union{String, Nothing}=nothing, yaxis_range::Union{Vector{<:Real}, Nothing}=nothing, slider::Bool=false)
    sp = Subplots(
        rows=length(unique(df."Line Model")), cols=length(cases(df)),
        shared_xaxes="all", shared_yaxes="all",
        start_cell = "top-left",
        column_titles=["Case $i ($key)" for (i, key) in enumerate(keys(cases(df)))],
        row_titles=["Static Pi Model", "Dynamic Pi Model"],
        x_title=x_title,
        y_title=y_title,
        horizontal_spacing=0.03,
        vertical_spacing=0.03,
    )

    layout_kwargs = Dict()
    if !isnothing(title_text) layout_kwargs[:title_text] = title_text end
    if !isnothing(yaxis_range) layout_kwargs[:yaxis_range] = yaxis_range end
    if slider
        values = sort(unique(df."Power Setpoint"))
        layout_kwargs[:sliders] = [attr(
            steps=[
                attr(
                    label=string(round(values[i], digits=3)),
                    method="restyle",
                    # args=["visible", [j%2==i%2 for j in 1:48]]
                    # args = [Dict(:visible=>(i%2==0)) for i in 1:48],
                    # args=[attr(z=(zz[i, :, :],), text=(text[i, :, :],)), 0]
                )
                for i in 1:length(values)
            ],
            active=1,
            currentvalue_prefix="Power Setpoint: ",
            pad_t=40,
        )]
    end

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
            filename=filename,
            height=900,
            width=1600,
            scale=1
        ).fields,
        displayModeBar=true,
    ))

    for ax in [:xaxis, :xaxis1, :xaxis2, :xaxis3, :xaxis4, :xaxis5, :xaxis6, :yaxis, :yaxis1, :yaxis2, :yaxis3, :yaxis4, :yaxis5, :yaxis6]
        fig.plot.layout[ax][:showticklabels] = true
    end
    fig.plot.layout.margin[:l] *= 2
    fig.plot.layout.margin[:b] *= 2
    fig.plot.layout.margin[:r] *= 2
    fig.plot.layout.margin[:t] *= 2
    return fig
end


################## PLOTTING FUNCTIONS #################
# Eigenvalue Plot
function eigplot(df::DataFrame)
    selectZIPE(df, η) = @subset(df, :z_percent.≈η[1], :i_percent.≈η[2], :p_percent.≈η[3])
    geteigs(df) = [
        (
            name=name, 
            eigs=Dict(selectZIPE(df, η)."Power Setpoint" .=> selectZIPE(df, η).Eigenvalues),
            # p_factors=Dict(selectZIPE(df, η)."Power Setpoint" .=> selectZIPE(df, η)."Participation Factors"),
        ) 
        for (name, η) in (η_combos)
    ]

    stateigs = [(name=casename, data=geteigs(@subset(casedf, :"Line Model" .== "statpi"))) for (casename, casedf) in cases(df)]
    dyneigs = [(name=casename, data=geteigs(@subset(casedf, :"Line Model" .== "dynpi"))) for (casename, casedf) in cases(df)]

    fig = makefig(
        df,
        L"\mathrm{Re}(\lambda)", 
        L"\mathrm{Im}(\lambda)", 
        "eigenvalue_plot", 
        "System Eigenvalues over Line Model, Injector Configuration, and ZIPE Parameters",
        nothing,
        true
    )
    # pwr_setpts = []
    xs = []
    # p_factors = []
    for η_combo_idx in 1:length(η_combos)
        for inj_case_idx in 1:length(cases(df))
            pwr_setpt = 1.0
            hovertext = "η="*string(η_combos[stateigs[inj_case_idx].data[η_combo_idx].name])
            # println(stateigs[inj_case_idx].name, stateigs[inj_case_idx].data[η_combo_idx].name, length(stateigs[inj_case_idx].data[η_combo_idx].eigs))
            η_case = stateigs[inj_case_idx].data[η_combo_idx].name
            η_case *= ": $(η_combos[η_case])"
            if stateigs[inj_case_idx].data[η_combo_idx].eigs isa Missing
                println("No data for statpi, case $inj_case_idx ($(stateigs[inj_case_idx].name)), $(stateigs[inj_case_idx].data[η_combo_idx].name)")
            else
                add_trace!(fig, PlotlyJS.scatter(
                    # x=real(stateigs[inj_case_idx].data[η_combo_idx].eigs[pwr_setpt]),
                    # y=imag(stateigs[inj_case_idx].data[η_combo_idx].eigs[pwr_setpt]),
                    mode="markers",
                    name="Case $(inj_case_idx), statpi",
                    legendgroup=η_case,
                    legendgrouptitle_text=η_case,
                    marker=attr(size=7, color=colorlist[η_combo_idx]),
                    opacity=0.7,
                    hovertext=hovertext,
                    ), row=1, col=inj_case_idx)
                push!(xs, stateigs[inj_case_idx].data[η_combo_idx].eigs)
                # push!(p_factors, stateigs[inj_case_idx].data[η_combo_idx].p_factors)
            end
                
            if dyneigs[inj_case_idx].data[η_combo_idx].eigs isa Missing
                println("No data for dynpi, case $inj_case_idx ($(dyneigs[inj_case_idx].name)), $(dyneigs[inj_case_idx].data[η_combo_idx].name)")
            else
                add_trace!(fig, PlotlyJS.scatter(
                    # x=real.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs[pwr_setpt])),
                    # y=imag.(reduce(vcat, dyneigs[inj_case_idx].data[η_combo_idx].eigs[pwr_setpt])),
                    mode="markers",
                    name="Case $(inj_case_idx), dynpi",
                    legendgroup=η_case,
                    legendgrouptitle_text=η_case,
                    marker=attr(size=7, color=colorlist[η_combo_idx]),
                    opacity=0.7,
                    hovertext=hovertext,
                    ), row=2, col=inj_case_idx)
                push!(xs, dyneigs[inj_case_idx].data[η_combo_idx].eigs)
                # push!(p_factors, dyneigs[inj_case_idx].data[η_combo_idx].p_factors)
            end
        end
    end
    for (step, pwr_setpt) in zip(fig.plot.layout[:sliders][1].steps, sort(unique(df."Power Setpoint")))
        step.args = [Dict(
            "x"=>[real(x[pwr_setpt]) for x in xs], 
            "y"=>[imag(x[pwr_setpt]) for x in xs],
            # "hovertext"=>[["p=$(round(p, sigdigits=4))" for p in p_vec[pwr_setpt]] for p_vec in p_factors],
        )]
    end

    return fig
end

"""
plot transient results from dataframe. must have columns :z_percent, :i_percent, :p_percent, :e_percent, :"Power Setpoint", "Load Voltage at Bus 5", "Line Model", 
"""
function transient(df::DataFrame)
    selectZIPE(df, η) = @subset(df, :z_percent.≈η[1], :i_percent.≈η[2], :p_percent.≈η[3])
    getvoltage(df) = [(name=name, data=Dict(selectZIPE(df, η)."Power Setpoint" .=> selectZIPE(df, η)."Load Voltage at Bus 5")) for (name, η) in (η_combos)]
    selectLineModel(df, linemodel) = [(name=casename, data=getvoltage(@subset(casedf, :"Line Model" .== linemodel))) for (casename, casedf) in cases(df)]

    # statpi = [(name=casename, data=getvoltage(@subset(casedf, :"Line Model" .== "statpi"))) for (casename, casedf) in cases(df)]
    # dynpi = [(name=casename, data=getvoltage(@subset(casedf, :"Line Model" .== "dynpi"))) for (casename, casedf) in cases(df)]
    line_models = [
        (name="statpi", data=statpi), 
        (name="dynpi", data=dynpi),
    ]
    line_models = [(name=name, data=selectLineModel(df, name)) for name in unique(df."Line Model")]

    fig = makefig(df, L"\mathrm{Time} \:\: \mathrm{[s]}", L"V_\text{rms} \:\: \mathrm{[\ p.u.]}", "transient_plot", "Transient Simulation: Voltage at Bus 5, Branch Trip on Line Bus 5-Bus 4", [0.0, 10.0], true)
    xs = []
    for η_combo_idx in 1:length(η_combos)
        for inj_case_idx in 1:length(cases(df))
            for line_model_idx in 1:length(line_models)
                hovertext = "η="*string(η_combos[line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].name])
                η_case = line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].name
                η_case *= ": $(η_combos[η_case])"
                if line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].data isa Missing
                    println("No data for $(line_models[line_model_idx].name), case $inj_case_idx ($(line_models[line_model_idx].data[inj_case_idx].name)), $(line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].name)")
                else
                    add_trace!(fig, PlotlyJS.scatter(
                        x0=0.48,# collect(0.48:0.00005:1.0),
                        dx=0.00005,
                        # DON'T plot y here, we don't want extra repeat data (all data is included in slider setup)
                        # y=round.(line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].data[1.0], sigdigits=4),
                        # mode="line",
                        name="Case $(inj_case_idx), $(line_models[line_model_idx].name)",
                        legendgroup=η_case,
                        legendgrouptitle_text=η_case,
                        marker=attr(size=7, color=colorlist[η_combo_idx]),
                        opacity=0.7,
                        hovertext=hovertext,
                        ), row=line_model_idx, col=inj_case_idx)
                    push!(xs, line_models[line_model_idx].data[inj_case_idx].data[η_combo_idx].data)

                end
            end
                
            # if dynpi[inj_case_idx].data[η_combo_idx].data isa Missing
            #     println("No data for dynpi, case $inj_case_idx ($(dynpi[inj_case_idx].name)), $(dynpi[inj_case_idx].data[η_combo_idx].name)")
            # else
            #     add_trace!(fig, PlotlyJS.scatter(
            #         x0=0.48,# collect(0.48:0.00005:1.0),
            #         dx=0.00005,
            #         # DON'T plot y here, we don't want extra repeat data (all data is included in slider setup)
            #         # y=round.(dynpi[inj_case_idx].data[η_combo_idx].data[1.0], sigdigits=4),
            #         # mode="lines",
            #         name="Case $(inj_case_idx), dynpi",
            #         legendgroup=η_case,
            #         legendgrouptitle_text=η_case,
            #         marker=attr(size=7, color=colorlist[η_combo_idx]),
            #         opacity=0.7,
            #         hovertext=hovertext,
            #         ), row=2, col=inj_case_idx)
            #     push!(xs, dynpi[inj_case_idx].data[η_combo_idx].data)
            # end
        end
    end
    for (step, pwr_setpt) in zip(fig.plot.layout[:sliders][1].steps, sort(unique(df."Power Setpoint")))
        step.args = ["y", [round.(x[pwr_setpt], sigdigits=6) for x in xs]]
    end
    return fig
end

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

function get_participation_factors(sim, bus_idx)
    sm = small_signal_analysis(sim)
    jacwrapper = get_jacobian(ResidualModel, sim.inputs, sm.operating_point, 0)
    eigen_vals, eigen_vecs = PSID._get_eigenvalues(jacwrapper.Jv, true)
    statemap = deepcopy(PSID.make_global_state_map(sim.inputs))
    statemap["busses"] = Dict([Symbol(string(i))=>i for i in 1:19]...)
    participation_factors = PSID._get_participation_factors(eigen_vecs, PSID._make_reduced_jacobian_index(statemap, [true for _ in 0:length(sim.results.solution.u[1])]))

    # covariance matrix of participation factors shows high correlation between v_i and v_r at each bus.
    # so we'll just take the average. (surely this is the participation factor towards the voltage magnitude)

    ## code to get cov matrix
    # p = participation_factors
    # P = hcat([p[Symbol(string(i))].-(sum(p[Symbol(string(i))])/length(p[Symbol(string(i))])) for i in 1:19]...)
    # for i in 1:19
    #     P[:, i] /= sqrt(sum(P[:, i].^2))
    # end
    # cov_p = P'*P
    return eigen_vals, (participation_factors["busses"][Symbol(string(bus_idx))] .+ participation_factors["busses"][Symbol(string(bus_idx+PSID.get_bus_count(sim.inputs)))]).*0.5
end

function add_eigs!(df)
    x = hcat([[a, b] for (a, b) in get_participation_factors.(df.sim, 5)]...)
    df.Eigenvalues .= x[1, :]
    df[!, :"Participation Factors"] .= real.(x[2, :]) # type is complex but im component is always zero. just from hcat making types funny
end


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
function plot(
    df::DataFrame;

    rows::Union{String, Symbol, Nothing}=nothing,
    cols::Union{String, Symbol, Nothing}=nothing,
    color::Union{String, Symbol, Nothing}=nothing,
    trace_names::Union{String, Symbol, Nothing}=nothing,
    hovertext::Union{String, Symbol, Nothing}=nothing,

    slider::Union{String, Symbol, Nothing}=nothing,
    slider_sort_func::Function = identity,
    slider_label_func::Function = x->round(x, sigdigits=2),
    slider_current_value_prefix::Union{String, LaTeXString, Nothing} = nothing,

    x::Union{String, Symbol, NamedTuple{(:x0, :dx), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    y::Union{String, Symbol, NamedTuple{(:y0, :dy), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    x_title::Union{String, LaTeXString, Nothing} = nothing,
    y_title::Union{String, LaTeXString, Nothing} = nothing,

    col_titles::Union{String, Symbol, Vector{Union{String, LaTeXString}}, Nothing}=nothing,
    row_titles::Union{String, Symbol, Vector{Union{String, LaTeXString}}, Nothing}=nothing,
    supertitle::Union{String, LaTeXString, Nothing} = nothing,

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
                # println([i for i in [x, y, hovertext] if i isa Union{String, Symbol}])
                push!(sliderdata, Dict((isnothing(slider) ? [i for i in 1:length(eachrow(data))] : data[!, slider]) .=> eachrow(select(data, [i for i in [x, y, hovertext, trace_names] if i isa Union{String, Symbol}]))))

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
                        if x isa Union{String, Symbol} scatargs[:x]=i[x] end
                        if y isa Union{String, Symbol} scatargs[:y]=i[y] end
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
        if x isa Union{String, Symbol}
            step.args[1]["x"]=[data[sliderval][x] for data in sliderdata]
        end
        if y isa Union{String, Symbol}
            step.args[1]["y"]=[data[sliderval][y] for data in sliderdata]
        end
        if !isnothing(trace_names) step.args[1]["name"]      = [data[sliderval][trace_names] for data in sliderdata] end
        if !isnothing(hovertext)   step.args[1]["hovertext"] = [data[sliderval][hovertext]   for data in sliderdata] end
    end
    return fig
end




η_names = Dict(
    [0.1, 0.1, 0.1, 0.7]     => "High E Load",
    [0.0, 0.0, 0.0, 1.0]     => "Full E Load",
    [0.3, 0.3, 0.3, 0.1]     => "Low E Load",
    [0.2, 0.2, 0.2, 0.4]     => "Medium E Load",
    [0.15, 0.15, 0.15, 0.55] => "Low P High E Load",
    [0.0, 0.0, 1.0, 0.0]     => "Constant Power",
    [0.15, 0.15, 0.55, 0.15] => "High P Low E Load",
    [1.0, 0.0, 0.0, 0.0]     => "Constant Impedance",
)
inj_case_names = Dict(
    "SM, GFL, GFM" => "Case 1",
    "SM, GFM, GFL" => "Case 2", 
    "GFM, GFL, SM" => "Case 3",
)

# first we add columns for all the data we want to include in the plot
getvec(x::LoadParams) = [x.z_percent, x.i_percent, x.p_percent, x.e_percent]
df[!, "ZIPE Parameters"] = (x->"$(η_names[round.(x, digits=3)]): $(round.(x, digits=3))").(getvec.(df."ZIPE Load Params"))
df[!, "hovertext"] = (x->"η=$x").(getvec.(df."ZIPE Load Params"))
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ")
df[!, "tracename"] = map(x->inj_case_names[x[1]]*", "*x[2], zip(df[!, "Injector Setup"], df[!, "Line Model"]))
df[!, "Injector Setup"] = map(x->inj_case_names[x]*" ($x)", df[!, "Injector Setup"])

# then we plot it!
plot(
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
    y="Load Voltage at Bus 5",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Voltage}\:\:[\mathrm{p.u.}]",

    col_titles="Injector Setup",
    row_titles="Line Model",
    supertitle="Voltage Plot",

    yaxis_home_range = (min=0, max=10),
    xaxis_home_range = nothing,

    image_export_filename = "transient_plot",
)


