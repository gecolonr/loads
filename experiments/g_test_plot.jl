using PlotlyJS
using DataFrames
using DataFramesMeta
using LaTeXStrings
using OrderedCollections
using Colors

path = "data/gab_tests/results0.jls"
df = load_serde_data(path::String)
df[!, "real_eigs"] = map(real, df[!, "Eigs"])
df[!, "imag_eigs"] = map(imag, df[!, "Eigs"])

include("plotting.jl")

p = makeplots(
    df,
    rows= "Line Model",
    cols="Power Setpoint",
    color="Line impedance increase",
    # legendgroup="hovertext",
    # opacity=0.7,
    # markershape="markershape",
    trace_names="Line Model",
    colorbar=true,
    # slider="Power Setpoint",
    # hovertext="hovertext_eig",
    # opacity="colorbar",
    # scattertext="scattertext",
    scattermode="markers",
    # markersize=7,
    # color_sort_func=x->(-x),

    # slider_current_value_prefix="Load Scale: ",

    x="real_eigs",
    y="imag_eigs",
    # x="time",
    # y="Bus 3 Injector Current",

    x_title=L"\Large{\mathrm{Re}(\lambda)}",
    y_title=L"\Large{\mathrm{Im}(\lambda)}",
    # col_title_func=x->LaTeXString("\$\\Large{\\text{Load Scale: }$(string(x))}\$"),
    col_title_func=x->"Load Scale: "*string(x),
    row_title_func=x->(x*" lines"),
    # col_title_func=x->Lastring(x)

    supertitle=L"\Large{\text{System Eigenvalues}}",
    image_export_filename="eigenvalue_plot",
    
    # colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11)),
    legend_visible=false,
    colorlist = vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))),
    fontsize=24,
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax= 1.0,
        # cmid=0.5,
        cmin= 1.0,
        title="",
        colorbar=attr(
            # bordercolor=colorant"black",
            # borderwidth=1,
            # title=attr(text="Load Scale"),
            # blahblahblah=[1,2,3,4,5],
            outlinecolor=colorant"black",
            outlinewidth=1,
            tickmode="array",
            ticktext=(["Z=1.0", "P=1.0", "P=0.8", "P=0.6", "P=0.4", "P=0.2", "P=E=0", "E=0.2", "E=0.4", "E=0.6", "E=0.8", "E=1.0"]),
            # tickprefix="\$",
            # ticksuffix="\$",

            tickvals=(collect(range(-1.2, 1, 13)).+(2.2/(2*13)))[1:end - 1],
            xref="paper",
            yref="paper",
            x=1.02,
            thickness=40,
            # xpad=20,
            
        ),
        # colorscale=[[0.0, RGB(colorant"yellow")], [1.0, RGB(colorant"purple")]],
        colorscale=vcat([[[(idx-1)/12, color], [idx/12, color]] for (idx, color) in enumerate(vcat(RGB(0, 1, 0), RGB.(range(colorant"red", colorant"blue", length=11))))]...),
    )),
    image_export_size=(height=800, width=1600),
)
savehtmlplot(p, "eigenvalues2")

# p = makeplots(
#     df;
	
#     rows="Line Model",
#     color="Line impedance increase",
#     slider="Power Setpoint",
#     slider_current_value_prefix="Power Setpoint: ",
	
#     x=(x0=0.48, dx=0.00005),
#     y="Load Voltage at Bus 5",
#     x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
#     y_title = L"\mathrm{Current}\:\:[\mathrm{p.u.}]",
	
#     supertitle="Bus 3 Injector Current",
#     yaxis_home_range = (min=0, max=10),
# )

savehtmlplot(p, "data/gab_tests/plot0.html")