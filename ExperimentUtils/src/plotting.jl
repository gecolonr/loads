using PlotlyJS
using DataFrames
using DataFramesMeta
using LaTeXStrings
using Colors

"""
```julia
makeplots(
    df::DataFrame;

    rows::Union{AbstractString, Symbol, Nothing}=nothing,
    cols::Union{AbstractString, Symbol, Nothing}=nothing,
    legendgroup::Union{AbstractString, Symbol, Nothing}=nothing,
    markershape::Union{AbstractString, Symbol, Nothing}=nothing,
    trace_names::Union{AbstractString, Symbol, Nothing}=nothing,
    color::Union{AbstractString, Symbol, Nothing}=nothing,
    opacity::Union{AbstractString, Symbol, Real}=1.0,
    markersize::Union{AbstractString, Symbol, Real}=8,
    hovertext::Union{AbstractString, Symbol, Nothing}=nothing,
    scattermode::Union{AbstractString, Symbol}="markers",
    scattertext::Union{AbstractString, Symbol, Nothing}=nothing,
    
    colorbar::Bool=false,
    map_to_colorbar::Bool=true,
    marker_line_width::Union{Real, Nothing} = nothing,

    slider_trace_id::Union{AbstractString, Symbol, Nothing}=nothing,

    row_sort_func::Function = identity,
    col_sort_func::Function = identity,
    color_sort_func::Function = identity,
    slider_sort_func::Function = identity,


    slider::Union{AbstractString, Symbol, Nothing}=nothing,
    slider_label_func::Function = x->round(x, sigdigits=2),
    slider_current_value_prefix::Union{AbstractString, Nothing} = nothing,
    x::Union{AbstractString, Symbol, NamedTuple{(:x0, :dx), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    y::Union{AbstractString, Symbol, NamedTuple{(:y0, :dy), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    x_title::Union{AbstractString, Nothing} = nothing,
    y_title::Union{AbstractString, Nothing} = nothing,

    col_title_func::Function=identity,
    row_title_func::Function=identity,
    legendgroup_title_func::Function=identity,
    supertitle::Union{AbstractString, Nothing} = nothing,

    yaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,
    xaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,

    image_export_filename::String = "saved_plot",
    image_export_size::@NamedTuple{height::Int64, width::Int64}=(height=1200, width=1200),
    colorlist::Union{Vector{String}, Vector{<:Colors.Colorant}} = Colors.distinguishable_colors(10),
    symbollist::Vector{<:AbstractString} = ["circle", "square", "diamond", "cross", "triangle-up", "star", "circle-cross", "y-up", "circle-open", "square-open", "diamond-open", "cross-open", "triangle-up-open", "star-open"],
    colorbar_args::Dict = Dict(attr(autocolorscale=true, colorbar=attr(outlinecolor=colorant"black", outlinewidth=1))),
    use_webgl::Bool = false,
    hide_legend_duplicates::Bool = true,
    legend_location::@NamedTuple{x::Float64, y::Float64} = (x=0.9, y=0.9),
    legend_visible::Bool=true,
    shared_xaxes::Union{String, Bool}="all",
    shared_yaxes::Union{String, Bool}="all",
    fontsize::Real=18,
    offsets::NamedTuple{(:xtitle, :ytitle, :coltitle, :rowtitle), <:Tuple{Real, Real, Real, Real}} = (xtitle = -0.03, ytitle=-0.05, coltitle=0.01, rowtitle=0.01),
    subplotgaps::NamedTuple{(:horizontal_spacing, :vertical_spacing), <:Tuple{Real, Real}}=(horizontal_spacing=0.1, vertical_spacing=0.06),
)
```

Plot data from a dataframe. allows high dimensional data through a grid of subplots, trace colors, and a slider.

+ 2 variables for (x, y) on plot
+ 2 variables for row and column in subplot grid
+ 1 variable for color
+ 1 variable for marker symbol
+ 1 variable on a slider
+ 1 variable on hovertext
+ 1 variable on marker text
+ 1 variable on opacity
+ 1 variable on trace name
+ 1 variable for trace or datapoint size/thickness
+ 0 variables for legend group (typically can't work independently; pair with color or marker symbol)
= 12 max variables

In reality, you should use less. For example, it's impossible to read the opacity correctly, so you should make hovertext and opacity the same.
Also, if you don't need that many variables, you can just leave them out. If you don't specify what to separate along subplot rows and columns, `makeplots` will just give you a single subplot.

## Notes
 - colorbar doesn't work very well with line plots. If you want to vary color by trace, that's ok, just set the `colorlist` to make the legend effectively a discrete colorbar.
 - for slider:
    - make sure that each slider value has the exact same number of traces
    - either provide the `slider_trace_id` or sort the dataframe such that subsetting to each slider value gives the remaining traces in the same order every time

# Arguments
## Mandatory Arguments 
- `df`: the DataFrame to get data from.

## DATA TO PLOT
`x` and `y` can be:
 - AbstractString, Symbol: get from column of df
 - NamedTuple of `(x0, dx)` or `(y0, dy)`: evenly spaced values with given parameters
 - Vector{<:Real}: this specific array every time

If you choose to get the series from the dataframe, each element of the column will be one series/trace, so it better be a vector!
If you have scalar or categorical values, try the default `PlotlyJS.plot(df, ...)` - it's pretty good for that.

## Specifying which columns to use for what
 - `rows`: the column of the dataframe specifying which row of subplots this data should be on
 - `cols`: the column of the dataframe specifying which column of subplots this data should be on
 - `legendgroup`: the column of the dataframe specifying which legendgroup this trace should be. Value of this column is the legendgroup name.
 - `color`: the column of the dataframe specifying the color of this trace (or the color of each datapoint in this trace). Mapped to `colorlist` after being sorted by `color_sort_func` if `colorbar==false`, otherwise used raw for colorbar.
 - `opacity` (default: 1.0): the opacity of all datapoints OR the column of the dataframe specifying the opacity of each datapoint
 - `markershape`: the column of the dataframe to map to marker shape. Uses shapes provided in `symbollist`.
 - `markersize` (default: 8): the size of all markers OR the column of the dataframe specifying the size of each trace or datapoint.
 - `trace_names`: the column of the dataframe specifying the text name of this trace. Visible in legend and on hover. If `hide_legend_duplicates==true`, duplicates of `trace_names` will not be seen in each legend group.
 - `hovertext`: the column of the dataframe specifying any additional text that appears for this trace on hover
 - `slider`: the column of the dataframe specifying which slider value this trace corresponds to. See **Slider Config**
 - `scattertext`: the column of the dataframe specifying the text to appear next to datapoints. Only has effect if `scattermode` contains the `text` flag.

## Titles
 - `x_title`: x axis title for every plot.
 - `y_title`: y axis title for every plot.
 - `col_title_func`: function to get column title from value in column of df passed in `cols`. Default: identity
 - `row_title_func`: function to get row title from value in column of df passed in `row`. Default: identity
 - `legendgroup_title_func`: function to get the legendgroup title text from the legendgroup name. Default: identity
 - `supertitle`: Biggest title on the graph.

## Legend & Axes
 - `colorbar` (default: false): whether to display the colorbar. Can be finnicky, only really works for non-line plots.
 - `colorbar_args`: Dictionary of parameters for the colorbar. See Julia PlotlyJS docs for Layout.coloraxis.
 - `shared_xaxes`, `shared_yaxes`: how to share axes. see `Subplots`.
 - `hide_legend_duplicates` (default: true): hide duplicate trace names for each legend group.
 - `legend_visible` (default: true): show the legend

## Markers & Colors
 - `scattermode` (default: "markers"): mode for scatterplot. Some combination of "markers", "lines", and "text", joined with "+". ie, "markers+text".
 - `colorlist` (default: 10 ok colors): list of colors to use. If you have more than 10 different values you want represented by color, or you want a specific set of colors, set this array.
 - `color_sort_func` (default: identity): function by which to sort values of `color` before associating with `colorlist` (if no colorbar)
 - `symbollist` (default: 14 different symbols): list of symbols to use for the different values of `markershape`. If you have more than 14 different values, set this array.
 - `use_webgl` (default: `true`): Whether to use scattergl or plain scatter.

## Sizes
 - `fontsize` (default 18): what do you think this option does lmao. take a guess
 - `margin` (default: 200): Int, just the margin size (px)
 - `offsets` (default: `(xtitle = -0.03, ytitle=-0.05, coltitle=0.01, rowtitle=0.01)`): where to put the titles, as a fraction of the plot
 - `yaxis_home_range` and `xaxis_home_range` can be `NamedTuple{(:min, :max), <:Tuple{Real, Real}}`
 - `legend_location` (default: (x=0.9, y=0.9)): NamedTuple with Float64 `x` and `y` defining the legend location relative to the paper.
 - `image_export_size` (default: (height=1200, width=1200)): NamedTuple with Integers `height` and `width` specifying the size in pixels of the SVG plots saved with the "download plot" button.
     - `image_export_filename` (default: `"saved_plot"`): default filename when image is exported using button on plot.

## Slider Config
 - `slider_sort_func` (default: `identity`): a function to get keys for the sort of the slider labels. For example, `(x -> -x)` would reverse the slider.
 - `slider_trace_id`: the column in the dataframe corresponding to a trace identifier. There must be exactly one row for each unique combination of this column and the slider column. If not passed, will use the order of the dataframe.
 - `slider_label_func` (default: `x->round(x, sigdigits=2)`): gets slider tick label from value.
 - `slider_current_value_prefix`: prefix to put before the value on the slider's label. For example, "val: " would make the label "val: 0.5" or something
"""
function makeplots(
    df::DataFrame;

    rows::Union{AbstractString, Symbol, Nothing}=nothing,
    cols::Union{AbstractString, Symbol, Nothing}=nothing,
    legendgroup::Union{AbstractString, Symbol, Nothing}=nothing,
    markershape::Union{AbstractString, Symbol, Nothing}=nothing,
    trace_names::Union{AbstractString, Symbol, Nothing}=nothing,
    color::Union{AbstractString, Symbol, Nothing}=nothing,
    opacity::Union{AbstractString, Symbol, Real}=1.0,
    markersize::Union{AbstractString, Symbol, Real}=8,
    hovertext::Union{AbstractString, Symbol, Nothing}=nothing,
    scattermode::Union{AbstractString, Symbol}="markers",
    scattertext::Union{AbstractString, Symbol, Nothing}=nothing,
    
    colorbar::Bool=false,
    map_to_colorbar::Bool=true,
    marker_line_width::Union{Real, Nothing} = nothing,

    slider_trace_id::Union{AbstractString, Symbol, Nothing}=nothing,

    row_sort_func::Function = identity,
    col_sort_func::Function = identity,
    color_sort_func::Function = identity,
    slider_sort_func::Function = identity,


    slider::Union{AbstractString, Symbol, Nothing}=nothing,
    slider_label_func::Function = x->round(x, sigdigits=2),
    slider_current_value_prefix::Union{AbstractString, Nothing} = nothing,
    x::Union{AbstractString, Symbol, NamedTuple{(:x0, :dx), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    y::Union{AbstractString, Symbol, NamedTuple{(:y0, :dy), <:Tuple{Real, Real}}, Vector{<:Real}}=nothing,
    x_title::Union{AbstractString, Nothing} = nothing,
    y_title::Union{AbstractString, Nothing} = nothing,

    col_title_func::Function=identity,
    row_title_func::Function=identity,
    legendgroup_title_func::Function=identity,
    supertitle::Union{AbstractString, Nothing} = nothing,

    yaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,
    xaxis_home_range::Union{NamedTuple{(:min, :max), <:Tuple{Real, Real}}, Nothing} = nothing,

    image_export_filename::String = "saved_plot",
    image_export_size::@NamedTuple{height::Int64, width::Int64}=(height=1200, width=1200),
    margin::Int = 200,
    colorlist::Union{Vector{String}, Vector{<:Colors.Colorant}} = Colors.distinguishable_colors(10),
    symbollist::Vector{<:AbstractString} = ["circle", "square", "diamond", "cross", "triangle-up", "star", "circle-cross", "y-up", "circle-open", "square-open", "diamond-open", "cross-open", "triangle-up-open", "star-open"],
    colorbar_args::Dict = Dict(attr(autocolorscale=true, colorbar=attr(outlinecolor=colorant"black", outlinewidth=1))),
    use_webgl::Bool = false,
    hide_legend_duplicates::Bool = true,
    legend_location::@NamedTuple{x::Float64, y::Float64} = (x=0.9, y=0.9),
    legend_visible::Bool=true,
    shared_xaxes::Union{String, Bool}="all",
    shared_yaxes::Union{String, Bool}="all",
    fontsize::Real=18,
    offsets::NamedTuple{(:xtitle, :ytitle, :coltitle, :rowtitle), <:Tuple{Real, Real, Real, Real}} = (xtitle = -0.03, ytitle=-0.05, coltitle=0.01, rowtitle=0.01),
    subplotgaps::NamedTuple{(:horizontal_spacing, :vertical_spacing), <:Tuple{Real, Real}}=(horizontal_spacing=0.1, vertical_spacing=0.06),
)

    if cols ∈ names(df)
        coldict = Dict(reverse.(enumerate(sort(unique(df[!, cols]), by=col_sort_func))))
        to_col_idx = c->coldict[c]
    else
        to_col_idx = c->1
    end
    if rows ∈ names(df)
        rowdict = Dict(reverse.(enumerate(sort(unique(df[!, rows]), by=row_sort_func))))
        to_row_idx = r->rowdict[r]
    else
        to_row_idx = r->1
    end

    to_color = identity
    to_marker = identity
    if color ∈ names(df)
        if occursin("lines", scattermode) || colorbar==false
            if eltype(df[!, color]) <: AbstractVector
                colorvals = unique(Iterators.flatten(df[!, color]))
                @warn "it is not known whether per-element color works when not using colorbar"
            else
                colorvals = unique(df[!, color])
            end
            @assert (length(colorvals) <= length(colorlist)) "need at least as many colors as unique color values!!"
            colorvals = sort(colorvals, by=color_sort_func)
            colordict = Dict(colorvals .=> first(colorlist, length(colorvals)))
            to_color(c) = colordict[c]
            # println(colordict)
            # return
        else
            to_color = identity
        end
    end

    if markershape ∈ names(df)
        if eltype(df[!, markershape]) <: AbstractVector
            markervals = unique(Iterators.flatten(df[!, markershape]))
            @assert length(markervals) <= length(symbollist)
        else
            markervals = unique(df[!, markershape])
            @assert length(markervals) <= length(symbollist)
        end
        markerdict = Dict(markervals .=> first(symbollist, length(markervals)))
        to_marker(m) = markerdict[m]
        println(markerdict)
    end
    # return
    scatterfunc = use_webgl ? PlotlyJS.scattergl : PlotlyJS.scatter
    # if !isnothing(color) && length(unique(df[!, color]))>length(colorlist)
    #     AssertionError("Not enough colors to give each series a unique color. Need at least $(length(unique(df[!, color]))).")
    # end

    spkwargs = Dict()
    if !isnothing(cols) spkwargs[:column_titles] = col_title_func.(sort(unique(df[!, cols]), by=col_sort_func)) end
    if !isnothing(rows) spkwargs[:row_titles] = row_title_func.(sort(unique(df[!, rows]), by=row_sort_func)) end

    # if !isnothing(col_titles) spkwargs[:column_titles] = col_titles isa Vector ? col_titles : unique(df[!, col_titles]) end
    # if !isnothing(row_titles) spkwargs[:row_titles] = row_titles isa Vector ? row_titles : unique(df[!, row_titles]) end
    if !isnothing(x_title) spkwargs[:x_title] = x_title end
    if !isnothing(y_title) spkwargs[:y_title] = y_title end

    sp = Subplots(
        rows=isnothing(rows) ? 1 : length(unique(df[!, rows])), 
        cols=isnothing(cols) ? 1 : length(unique(df[!, cols])),
        shared_xaxes=shared_xaxes, shared_yaxes=shared_yaxes,
        start_cell = "top-left",
        horizontal_spacing=subplotgaps.horizontal_spacing,
        vertical_spacing=subplotgaps.vertical_spacing;
        spkwargs...
    )
    

    layout_kwargs = Dict()
    layout_kwargs[:font]=attr(
        family="Computer Modern",
        size=fontsize,
    )
    # layout_kwargs[:editable] = true
    layout_kwargs[:paper_bgcolor] = "white"
    layout_kwargs[:plot_bgcolor] = "white"
    layout_kwargs[:yaxis] = attr(automargin=true, showline=true, linewidth=1, linecolor="black", mirror=true, gridcolor="#eeeeee", griddash="dot", zeroline=false)
    layout_kwargs[:xaxis] = attr(automargin=true, showline=true, linewidth=1, linecolor="black", mirror=true, gridcolor="#eeeeee", griddash="dot", zeroline=false)
    # layout_kwargs[:coloraxis] = attr(autocolorscale=true, colorbar=attr(outlinecolor=colorant"black", outlinewidth=1, xref="container", yref="container"))
    layout_kwargs[:coloraxis] = colorbar_args
    # layout_kwargs.coloraxis.colorbar.labelalias= attr(; layout_kwargs.coloraxis.colorbar.labelalias...)
    if colorbar 
        layout_kwargs[:legend] = attr(
            orientation="h", 
            yref="paper", 
            xref="paper", 
            x=legend_location.x, 
            y=legend_location.y,
            xanchor="right",
            yanchor="bottom",
        ) 
    end
    layout_kwargs[:showlegend] = legend_visible

    if !isnothing(supertitle) 
        layout_kwargs[:title_text] = supertitle 
    end
    # layout_kwargs[:xaxis_automargin]=true
    # layout_kwargs[:yaxis_automargin]=true
    if !isnothing(xaxis_home_range) layout_kwargs[:xaxis_range] = [xaxis_home_range.min, xaxis_home_range.max] end
    if !isnothing(yaxis_home_range) layout_kwargs[:yaxis_range] = [yaxis_home_range.min, yaxis_home_range.max] end

    slider_values = isnothing(slider) ? [0.0] : sort(unique(df[!, slider]), by=slider_sort_func)
    layout_kwargs[:sliders] = [attr(
        steps=[
            attr(
                label=string(slider_label_func(slider_values[i])),
                method="restyle",
            )
            for i in 1:length(slider_values)
        ],
        active=1,
        currentvalue_prefix=slider_current_value_prefix,
        # pad_t=40,
        pad = attr(b=30, l=0, r=0, t=120),
    )]
    layout_kwargs[:margin] = attr(
        b=margin,
        l=margin,
        pad=0,
        r=margin,
        t=margin,
    )
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
            height=image_export_size.height,
            width=image_export_size.width,
            scale=1
        ).fields,
        displayModeBar=true,
        # queuelength=5,
        # editable=true
    ))
    update_annotations!(fig, font_size=fontsize)
    for i in fig.plot.layout.annotations
        if i.text == x_title
            println("here1, $(i.x)")
            i.y += offsets.xtitle
        elseif i.text == y_title
            println("here2, $(i.y)")
            i.x += offsets.ytitle
        elseif rows ∈ names(df) && i.text ∈ row_title_func.(unique(df[!, rows]))
            println("here3, $(i.x)")
            i.x += offsets.rowtitle
        elseif cols ∈ names(df) && i.text ∈ col_title_func.(unique(df[!, cols]))
            println("here4, $(i.y)")
            i.y += offsets.coltitle
        end
    end
    legend_seen_already = []
    layoutkeys = keys(fig.plot.layout)
    layoutkeys = collect(layoutkeys)[(x->(occursin("axis", string(x)))).(layoutkeys)]
    for ax in layoutkeys
        fig.plot.layout[ax][:showticklabels] = true
    end
    # fig.plot.layout.margin[:l] *= 4
    # fig.plot.layout.margin[:b] *= 4
    # fig.plot.layout.margin[:r] *= 4
    # fig.plot.layout.margin[:t] *= 4

    if slider ∈ names(df)
        slidervals = sort(unique(df[!, slider]), by=slider_sort_func)
        df_sliderfilt = df[df[!, slider].==last(slidervals), :]

        if slider_trace_id ∈ names(df) @assert nrow(df_sliderfilt)==length(unique(df_sliderfilt[!, slider_trace_id])) "conflict detected: slider invariant has duplicates for a slider value" end
    else
        df_sliderfilt = df
    end
    sliderdata = []

    for (idx, trace) in enumerate(eachrow(df_sliderfilt))
        if slider ∈ names(df)
            push!(sliderdata, 
                if isnothing(slider_trace_id)
                    DataFrames.DataFrame([eachrow(df[df[!, slider].==val, :])[idx] for val in slidervals])
                else
                    df[df[!, slider_trace_id].==trace[!, slider_trace_id], :]
                end
            )
        end

        scatargs = Dict()
        scatargs[:marker] = attr()
        if marker_line_width isa Real
            scatargs[:marker].line_width=marker_line_width
        end
        if (x isa Vector) 
            scatargs[:x] = x
        elseif (x isa NamedTuple) 
            scatargs[:x0] = x.x0
            scatargs[:dx] = x.dx
        else
            scatargs[:x] = trace[x]
        end

        if (y isa Vector) 
            scatargs[:y] = y
        elseif (y isa NamedTuple)
            scatargs[:y0] = y.y0
            scatargs[:dy] = y.dy
        else
            scatargs[:y] = trace[y]
        end
        
        if legendgroup ∈ names(df)
            scatargs[:legendgroup] = trace[legendgroup]
            scatargs[:legendgrouptitle_text] = legendgroup_title_func(trace[legendgroup])
            if hide_legend_duplicates
                if (trace[legendgroup], if trace_names in names(trace) trace[trace_names] else nothing end) in legend_seen_already
                    scatargs[:showlegend] = false
                else
                    push!(legend_seen_already, (trace[legendgroup], if trace_names in names(trace) trace[trace_names] else nothing end))
                end
            end
        end
        if markershape ∈ names(df) scatargs[:marker].symbol = to_marker(trace[markershape]) end
        if trace_names ∈ names(df) scatargs[:name] = trace[trace_names] end
        if color ∈ names(df)       
            scatargs[:marker].color = to_color.(trace[color])
            scatargs[:marker].line = attr(color=to_color.(trace[color]))
            # scatargs[:marker].line_color = to_color.(trace[color])
            # if marker_line_width isa Real
            #     scatargs[:marker].line.width = marker_line_width
            # end
            if colorbar
                # scatargs[:showlegend]=false
                if occursin("lines", scattermode)
                    # scatargs[:line_color]=trace[color]
                    scatargs[:line_coloraxis]="coloraxis"
                    scatargs[:line] = attr(
                        color=to_color(first(trace[color]))
                    )
                end
                if trace[color] isa Real && occursin("markers", scattermode)
                    scatargs[:marker].color = [trace[color] for _ in 1:length(if :x in keys(scatargs) scatargs[:x] else scatargs[:y] end)]
                end
                if map_to_colorbar && !occursin("lines", scattermode)
                    scatargs[:marker].coloraxis = "coloraxis"
                    scatargs[:marker].showscale = true
                end
                # scatargs[:marker].line.coloraxis = "coloraxis"
            else
                scatargs[:legendrank]=color_sort_func(trace[color])
            end
        end
        scatargs[:opacity] = if opacity ∈ names(df)     
            trace[opacity] 
        else
            opacity
        end
        scatargs[:marker].size = if markersize ∈ names(df)  
           trace[markersize] 
        else
            markersize
        end
        if hovertext ∈ names(df)   scatargs[:hovertext] = trace[hovertext] end
        # if visible ∈ names(df)     scatargs[:visible] = trace[visible] end
        scatargs[:mode] = if scattermode ∈ names(df) 
            trace[scattermode]
        else
            scattermode
        end
        if scattertext ∈ names(df) scatargs[:text] = trace[scattertext] end

        if markershape ∈ names(df) scatargs[:marker].symbol = to_marker(trace[markershape]) end
        # println(scatargs[:marker])
        PlotlyJS.add_trace!(fig, scatterfunc(; scatargs...), 
            row = if (rows in names(df)) to_row_idx(trace[rows]) else 1 end, 
            col = if (cols in names(df)) to_col_idx(trace[cols]) else 1 end,
        )
    end
    if slider ∉ names(df)
        if colorbar && occursin("lines", scattermode)
            PlotlyJS.add_trace!(fig, scatterfunc(x=[0,0,0], y=[0,0,0], visible="legendonly", mode="markers", marker=attr(coloraxis="coloraxis", showscale=true)))
        end
        return fig
    end
    # for data in sliderdata
    #     Base.show(data)
    # end
    col_updates_on_slider(col) = col ∈ names(df) && any((length(unique(data[!, col])) != 1) for data in sliderdata)


    for (step, sliderval) in zip(fig.plot.layout[:sliders][1].steps, sort(unique(df[!, slider]), by=slider_sort_func))
        step.args = [Dict()]
        args = step.args[1]
        getcol = colname->[first(data[data[!, slider].==sliderval, :][!, colname]) for data in sliderdata]
        # return getcol(x)
        if col_updates_on_slider(x)           args["x"] = map(z->round.(z, sigdigits=6), getcol(x)); println("reached x") end
        if col_updates_on_slider(y)           args["y"] = map(z->round.(z, sigdigits=6), getcol(y)); println("reached y") end
        if col_updates_on_slider(legendgroup) args["legendgroup"] = getcol(legendgroup); println("reached legendgroup") end
        if col_updates_on_slider(legendgroup) args["legendgrouptitle_text"] = getcol(legendgroup); println("reached legendgroup") end
        if col_updates_on_slider(markershape) args["marker.symbol"] = to_marker.(getcol(markershape)); println("reached markershape") end
        if col_updates_on_slider(trace_names) args["name"] = getcol(trace_names); println("reached trace_names") end
        if col_updates_on_slider(color)       args["marker.color"] = getcol(color); println("reached color")
                                              args["line.color"] = getcol(color);
                                            #   args["marker.line.color"] = getcol(color);
        end
        if col_updates_on_slider(opacity)     args["opacity"] = getcol(opacity); println("reached opacity") end
        if col_updates_on_slider(markersize)  args["marker.size"] = getcol(markersize); println("reached markersize") end
        if col_updates_on_slider(hovertext)   args["hovertext"] = getcol(hovertext); println("reached hovertext") end
        # if col_updates_on_slider(visible)     args["visible"] = getcol(visible); println("reached visible") end
        if col_updates_on_slider(scattermode) args["mode"] = getcol(scattermode); println("reached scattermode") end
        if col_updates_on_slider(scattertext) args["text"] = getcol(scattertext); println("reached scattertext") end
        # args["marker.line.color"] = getcol(color)
    end
    if colorbar && occursin("lines", scattermode)
        PlotlyJS.add_trace!(fig, scatterfunc(x=[0,0,0], y=[0,0,0], visible="legendonly", mode="markers", marker=attr(coloraxis="coloraxis", showscale=true)))
    end
    return fig
end



"""
    savehtmlplot(plot, filename::String=nothing)

utility to save `plot` to `filename`.html. Default filename is the plot's image download filename, if it exists, or worst case just "plot.html"
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
        PlotlyBase.to_html(io, plot, full_html=true)
    end
end

