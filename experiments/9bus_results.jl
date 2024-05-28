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
        rows=2, cols=3,
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
            eigs=Dict(selectZIPE(df, η)."Power Setpoint" .=> selectZIPE(df, η).Eigenvalues)
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
    # ys = []
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
                
            end
        end
    end
    for (step, pwr_setpt) in zip(fig.plot.layout[:sliders][1].steps, sort(unique(df."Power Setpoint")))
        step.args = [Dict("x"=>[real(x[pwr_setpt]) for x in xs], "y"=>[imag(x[pwr_setpt]) for x in xs])]
    end

    return fig
end

function transient(df::DataFrame)
    selectZIPE(df, η) = @subset(df, :z_percent.≈η[1], :i_percent.≈η[2], :p_percent.≈η[3])
    getvoltage(df) = [(name=name, data=Dict(selectZIPE(df, η)."Power Setpoint" .=> selectZIPE(df, η)."Load Voltage at Bus 5")) for (name, η) in (η_combos)]

    statpi = [(name=casename, data=getvoltage(@subset(casedf, :"Line Model" .== "statpi"))) for (casename, casedf) in cases(df)]
    dynpi = [(name=casename, data=getvoltage(@subset(casedf, :"Line Model" .== "dynpi"))) for (casename, casedf) in cases(df)]

    fig = makefig(df, L"\mathrm{Time} \:\: \mathrm{[s]}", L"V_\text{rms} \:\: \mathrm{[\ p.u.]}", "transient_plot", "Transient Simulation: Voltage at Bus 5, Branch Trip on Line Bus 5-Bus 4", [0.0, 10.0], true)
    xs = []
    for η_combo_idx in 1:length(η_combos)
        for inj_case_idx in 1:length(cases(df))
            hovertext = "η="*string(η_combos[statpi[inj_case_idx].data[η_combo_idx].name])
            η_case = statpi[inj_case_idx].data[η_combo_idx].name
            η_case *= ": $(η_combos[η_case])"
            if statpi[inj_case_idx].data[η_combo_idx].data isa Missing
                println("No data for statpi, case $inj_case_idx ($(statpi[inj_case_idx].name)), $(statpi[inj_case_idx].data[η_combo_idx].name)")
            else
                add_trace!(fig, PlotlyJS.scatter(
                    x0=0.48,# collect(0.48:0.00005:1.0),
                    dx=0.00005,
                    # DON'T plot y here, we don't want extra repeat data (all data is included in slider setup)
                    # y=round.(statpi[inj_case_idx].data[η_combo_idx].data[1.0], sigdigits=4),
                    # mode="line",
                    name="Case $(inj_case_idx), statpi",
                    legendgroup=η_case,
                    legendgrouptitle_text=η_case,
                    marker=attr(size=7, color=colorlist[η_combo_idx]),
                    opacity=0.7,
                    hovertext=hovertext,
                    ), row=1, col=inj_case_idx)
                push!(xs, statpi[inj_case_idx].data[η_combo_idx].data)

            end
                
            if dynpi[inj_case_idx].data[η_combo_idx].data isa Missing
                println("No data for dynpi, case $inj_case_idx ($(dynpi[inj_case_idx].name)), $(dynpi[inj_case_idx].data[η_combo_idx].name)")
            else
                add_trace!(fig, PlotlyJS.scatter(
                    x0=0.48,# collect(0.48:0.00005:1.0),
                    dx=0.00005,
                    # DON'T plot y here, we don't want extra repeat data (all data is included in slider setup)
                    # y=round.(dynpi[inj_case_idx].data[η_combo_idx].data[1.0], sigdigits=4),
                    # mode="lines",
                    name="Case $(inj_case_idx), dynpi",
                    legendgroup=η_case,
                    legendgrouptitle_text=η_case,
                    marker=attr(size=7, color=colorlist[η_combo_idx]),
                    opacity=0.7,
                    hovertext=hovertext,
                    ), row=2, col=inj_case_idx)
                push!(xs, dynpi[inj_case_idx].data[η_combo_idx].data)
            end
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
