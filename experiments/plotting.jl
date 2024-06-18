using PlotlyJS
using DataFrames
using DataFramesMeta
using LaTeXStrings
using OrderedCollections


"""
Plot data from a dataframe. allows high dimensional data through a grid of subplots, trace colors, and a slider.

# Arguments
## Mandatory Arguments 
- `df`: the DataFrame to get data from.

## DATA TO PLOT
`x` and `y` can be:
 - AbstractString, Symbol: get from column of df
 - NamedTuple of `(x0, dx)` or `(y0, dy)`: evenly spaced values with given parameters
 - Vector{<:Real}: this specific array every time

If you choose to get the series from the dataframe, each element of the column will be one series/trace, so it better be a vector!

## Specifying which columns to use for what
 - `rows`: the column of the dataframe specifying which row of subplots this data should be on
 - `cols`: the column of the dataframe specifying which column of subplots this data should be on
 - `color`: the column of the dataframe specifying which color this trace should be.
 - `trace_names`: the column of the dataframe specifying the text name of this trace
 - `hovertext`: the column of the dataframe specifying any additional text that appears for this trace on hover
 - `slider`: the column of the dataframe specifying which slider value this trace corresponds to. See **Slider Config**

## Titles
 - `x_title`: x axis title. must be the same for every plot.
 - `y_title`: y axis title. must be the same for every plot.
 - `col_title_func`: function to get column title from value in column of df passed in `cols`. Default: identity
 - `row_title_func`: function to get row title from value in column of df passed in `row`. Default: identity
 - `supertitle`: Biggest title on the graph.

## Other graph options
 - `yaxis_home_range` and `xaxis_home_range` can be NamedTuple{(:min, :max), <:Tuple{Real, Real}}
 - `scatterplot_args`: Dict to be passed to PlotlyJS.scatter() with any additional user-desired args.
 - `image_export_filename`: default filename when image is exported using button on plot.
 - `colorlist`: list of colors to use. If you have more than 10 different values you want represented by color, set the `color` array.
 
## Slider Config
 - `slider_sort_func`: a function to get keys for the sort of the slider labels. For example, (x -> -x) would reverse the slider.
 - `slider_label_func`: gets slider tick label from value. defaults to `x->round(x, sigdigits=2)`
 - `slider_current_value_prefix`: prefix to put before the value on the slider's label. For example, "val: " would make the label "val: 0.5" or something
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

    # col_titles::Union{AbstractString, Symbol, Vector{AbstractString}, Nothing}=nothing,
    col_title_func::Function=identity,
    row_title_func::Function=identity,
    # row_titles::Union{AbstractString, Symbol, Vector{AbstractString}, Nothing}=nothing,
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
    if !isnothing(cols) spkwargs[:column_titles] = col_title_func.(unique(df[!, cols])) end
    if !isnothing(rows) spkwargs[:row_titles] = row_title_func.(unique(df[!, rows])) end

    # if !isnothing(col_titles) spkwargs[:column_titles] = col_titles isa Vector ? col_titles : unique(df[!, col_titles]) end
    # if !isnothing(row_titles) spkwargs[:row_titles] = row_titles isa Vector ? row_titles : unique(df[!, row_titles]) end
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
                    scatargs[:legendgrouptitle_text]=colorval
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