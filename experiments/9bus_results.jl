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

include("sysbuilder.jl")
include("plotting.jl")



getvec(x::LoadParams) = round.([x.z_percent, x.i_percent, x.p_percent, x.e_percent], digits=3)

df = load_serde_data("/data/reiddye/loads/forplot_fine_twocase").df
# if neccessary, filter injector setup to only get the ones you want

df[!, "ZIPE Load Params"] = getvec.(df[!, "ZIPE Load Params"]);
df[!, "markershape"] = map(x-> if x[1] ≈ 1.0 L"Constant Impedance" elseif x[3] > 0 L"\Large{\eta_P > 0}" elseif x[4] > 0 L"\Large{\eta_E > 0}" else "E=P=0" end, df[!, "ZIPE Load Params"]);
df[!, "legendgroup"] = map(x-> if x[1] ≈ 1.0 L"{}{}{}" elseif x[3] > 0 L"{}" elseif x[4] > 0 L"{}{}" else L"{{}}" end, df[!, "ZIPE Load Params"]);
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"],
                                             df[!, "injector at {Bus 2}"],
                                             df[!, "injector at {Bus 3}"])), ", ");
df[!, "colorbar"] = map(x->((x[1]≈1.0) ? 1.2 : (x[3]-x[4])), df[!, "ZIPE Load Params"]);
df[!, "real"] = real.(df.Eigenvalues);
df[!, "imag"] = imag.(df.Eigenvalues);
df[!, "hovertext"] = "η=" .* string.(df[!, "ZIPE Load Params"]);
df[!, "e_percent"] = map(last, df[!, "ZIPE Load Params"]);
df[!, "transient_legendgroup_hack"] = map(x->LaTeXString("\$\\iffalse $(x) \\fi\$"), df[!, "hovertext"]);
df[!, "power_percent"] = map(x->(1.0-2*first(x)), df[!, "ZIPE Load Params"]);
function pfactors(sm)
    pf = summary_participation_factors(sm)
    return Dict(pf.Name .=> eachrow(select(pf, Not(:Name))))
end
df[!, "participation_factors"] = map(sm->get(pfactors(sm), "generator-3-1 ir_cnv", missing), df.sm);
df[!, "participation_factors"] = map(x->x isa Missing ? "" : "\nPF: ".*string.(round.((collect(x)), sigdigits=5)), df.participation_factors);
df[!, "hovertext_eig"] = [i.*j for (i, j) in zip(df[!, "hovertext"], df[!, "participation_factors"])];
df[!, "reverse_colorbar"] = -df[!, "colorbar"];
df = sort(df, ["Power Setpoint", "colorbar", "Line Model", "hovertext", "e_percent", "power_percent", "Injector Setup"]);
p = makeplots(
    df;
    rows="Line Model",
    cols="Injector Setup",
    color="reverse_colorbar",
    legendgroup="transient_legendgroup_hack",
    legendgroup_title_func=x->nothing,
    # markershape="markershape",
    trace_names="hovertext",
    colorbar=true,
    legend_visible=false,
    slider="Power Setpoint",
    hovertext="hovertext",
    # opacity=0.7,
    scattermode="lines",
    markersize=7,
    color_sort_func=identity,

    slider_current_value_prefix="Load Scale: ",

    # x="real",
    # y="imag",
    x="time",
    y="Bus 3 Injector Current",

    x_title=L"\Large{\mathrm{Time}\:\: [\mathrm{s}]}",
    y_title = L"\Large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",
    # row_title_func=x->LaTeXString("\$\\Large{\\text{$x lines}}\$"),
    row_title_func = x->(x*" lines"),

    supertitle=L"\Large{\text{Transient Current at Bus 3}}",
    image_export_filename="transient_plot",
    colorlist = vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))),
    shared_yaxes="rows",
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax= 1.0,
        # cmid=0.5,
        cmin= -1.2,
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

            tickvals=(collect(range(-1.2, 1, 13)).+(2.2/(2*12)))[1:end - 1],
            xref="paper",
            yref="paper",
            x=1.02,
            thickness=45
            # xpad=20,
            
        ),
        # colorscale=[[0.0, RGB(colorant"yellow")], [1.0, RGB(colorant"purple")]],
        colorscale=vcat([[[(idx-1)/12, color], [idx/12, color]] for (idx, color) in enumerate(vcat(RGB(0, 1, 0), RGB.(range(colorant"red", colorant"blue", length=11))))]...),
    )),
    fontsize=24,
)
savehtmlplot(p, "transient_niceified_fine_bus3")

p = makeplots(
    df;
    rows="Line Model",
    cols="Injector Setup",
    color="reverse_colorbar",
    legendgroup="hovertext",
    # markershape="markershape",
    trace_names="Line Model",
    colorbar=true,
    # opacity=1.0,
    # marker_line_width=1,
    slider="Power Setpoint",
    hovertext="hovertext_eig",
    # opacity="colorbar",
    # scattertext="scattertext",
    scattermode="markers",
    # markersize=7,
    color_sort_func=x->(-x),

    slider_current_value_prefix="Load Scale: ",

    x="real",
    y="imag",
    # x="time",
    # y="Bus 3 Injector Current",

    x_title=L"\Large{\mathrm{Re}(\lambda)}",
    y_title=L"\Large{\mathrm{Im}(\lambda)}",
    # row_title_func=x->LaTeXString("\$\\Large{\\text{$x lines}}\$"),
    row_title_func=x->(x*" lines"),

    supertitle=L"\Large{\text{System Eigenvalues}}",
    image_export_filename="eigenvalue_plot",
    
    # colorlist = "#" .* hex.(range(colorant"blue", colorant"red", length=11)),
    legend_visible=false,
    # map_to_colorbar=false,
    colorlist = vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))),
    fontsize=24,
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax= 1.0,
        # cmid=0.5,
        cmin= -1.2,
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
            width=45,
            # xpad=20,
            
        ),
        # colorscale=[[0.0, RGB(colorant"yellow")], [1.0, RGB(colorant"purple")]],
        colorscale=vcat([[[(idx-1)/12, color], [idx/12, color]] for (idx, color) in enumerate(vcat(RGB(0, 1, 0), RGB.(range(colorant"red", colorant"blue", length=11))))]...),
    )),
)
savehtmlplot(p, "eigplot_niceified_fine")


# p = makeplots(
#     @subset df begin
#         :"power_percent" .> 0.0
#         #* UNCOMMENT NEXT LINE TO EXCLUDE HIGH LOAD SCALES
#         # :"Power Setpoint" .< 0.37
#     end;
#     rows="Line Model",
#     # cols="markershape",
#     color="Power Setpoint",
#     legendgroup="legendgroup",
#     markershape="markershape",
#     trace_names="markershape",
#     colorbar=true,
#     slider="power_percent",
#     hovertext="hovertext_eig",
#     # opacity=0.7,
#     scattermode="markers",
#     markersize=10,

#     slider_current_value_prefix="P/E Fraction: ",

#     x="real",
#     y="imag",
#     # x="time",
#     # y="Bus 3 Injector Current",

#     x_title=L"\Large{\mathrm{Re}(\lambda)}",
#     y_title=L"\Large{\mathrm{Im}(\lambda)}",
#     row_title_func=x->(x*" lines"),
#     # row_title_func=x->LaTeXString("\$\\Large{\\text{$x lines}}\$"),

#     supertitle=L"\Large{\text{System Eigenvalues}}",
#     image_export_filename="eigenvalue_plot",
#     # colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11)),
#     legend_location=(x=0.97, y=1.02),
#     colorbar_args=Dict(attr(
#         autocolorscale=false, 
#         cmax=1.0,
#         # cmid=0.5,
#         cmin=0.05,
#         # title="Load Scale",
#         colorbar=attr(
#             # bordercolor=colorant"black",
#             # borderwidth=1,
#             title=attr(text="Load Scale"),
#             outlinecolor=colorant"black",
#             outlinewidth=1,
#             # tickmode="array",
#             # ticktext=["P=1.0", "P=0.8", "P=0.6", "P=0.4", "P=0.2", "P=E=0", "E=0.2", "E=0.4", "E=0.6", "E=0.8", "E=1.0"],
#             # tickvals=collect(-1.0:0.2:1.0),
#             xref="paper",
#             yref="paper",
#             x=1.02
#             # xpad=20,
            
#         ),
#         # colorscale=[[-1.0, RGB(colorant"red")], [1.0, RGB(colorant"blue")]],
#         colorscale=zip(0:0.1:1, RGB.(range(colorant"red", colorant"blue", length=11))),

#     )),
# )
# savehtmlplot(p, "eigplot_vsloadscale_fine")


p = makeplots(
    @subset df begin
        (:"Power Setpoint".≈0.2) .|| (:"Power Setpoint".≈0.5) .|| (:"Power Setpoint".≈0.8)
        :"Line Model" .== "dynpi"
    end;
    # rows="Line Model",
    cols="Power Setpoint",
    color="reverse_colorbar",
    legendgroup="hovertext",
    # opacity=0.7,
    # markershape="markershape",
    trace_names="Line Model",
    colorbar=true,
    # slider="Power Setpoint",
    hovertext="hovertext_eig",
    # opacity="colorbar",
    # scattertext="scattertext",
    scattermode="markers",
    # markersize=7,
    # color_sort_func=x->(-x),

    slider_current_value_prefix="Load Scale: ",

    x="real",
    y="imag",
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
        cmin= -1.2,
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
            width=45,
            # xpad=20,
            
        ),
        # colorscale=[[0.0, RGB(colorant"yellow")], [1.0, RGB(colorant"purple")]],
        colorscale=vcat([[[(idx-1)/12, color], [idx/12, color]] for (idx, color) in enumerate(vcat(RGB(0, 1, 0), RGB.(range(colorant"red", colorant"blue", length=11))))]...),
    )),
    image_export_size=(height=800, width=1600),
)
savehtmlplot(p, "eigplot1")

p = makeplots(
    @subset df begin
        (:"Power Setpoint".≈0.2) .|| (:"Power Setpoint".≈0.5) .|| (:"Power Setpoint".≈0.8)
        :"Line Model" .== "dynpi"
    end;
    # rows="Line Model"
    cols="Power Setpoint",
    # cols="Injector Setup",
    color="reverse_colorbar",
    legendgroup="transient_legendgroup_hack",
    legendgroup_title_func=x->nothing,
    # markershape="markershape",
    trace_names="hovertext",
    colorbar=true,
    legend_visible=false,
    # slider="Power Setpoint",
    hovertext="hovertext",
    # opacity=0.7,
    scattermode="lines",
    markersize=7,
    color_sort_func=identity,

    slider_current_value_prefix="Load Scale: ",

    # x="real",
    # y="imag",
    x="time",
    y="Bus 3 Injector Current",

    x_title=L"\Large{\mathrm{Time}\:\: [\mathrm{s}]}",
    y_title = L"\Large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",
    # row_title_func=x->LaTeXString("\$\\Large{\\text{$x lines}}\$"),
    row_title_func = x->(x*" lines"),
    col_title_func=x->"Load Scale: "*string(x),

    supertitle=L"\Large{\text{Transient Current at Bus 3}}",
    image_export_filename="transient_plot",
    # colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11)),
    shared_yaxes="columns",
    colorlist = vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))),
    fontsize=24,
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax= 1.0,
        # cmid=0.5,
        cmin= -1.2,
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
            width=45,
            # xpad=20,
            
        ),
        # colorscale=[[0.0, RGB(colorant"yellow")], [1.0, RGB(colorant"purple")]],
        colorscale=vcat([[[(idx-1)/12, color], [idx/12, color]] for (idx, color) in enumerate(vcat(RGB(0, 1, 0), RGB.(range(colorant"red", colorant"blue", length=11))))]...),
    )),
    image_export_size=(height=800, width=1600),
)
savehtmlplot(p, "transient1")


df2 = @subset df begin
    :reverse_colorbar .> -1.1
end
p = PlotlyJS.plot(PlotlyJS.scatter(;
    x=df2[!, "reverse_colorbar"], 
    y=df2[!, "dt"]/1e9, 
    marker=attr(color=map(x->((x=="SM, SM, SM") ? hex(colorant"red") : hex(colorant"blue")), df2[!, "Injector Setup"])),
    mode="markers",
),
Layout(title=attr(text="red=SM/SM/SM, blue=SM/GFM/GFL"), xaxis=attr(title=L"\eta_E-\eta_P"), yaxis=attr(title="Simulation Runtime (s)")))

df2[!, L"\eta_E-\eta_P"] = df2[!, "reverse_colorbar"]
df2[!, "Runtime (s)"] = df2[!, "dt"]/1e9
p = PlotlyJS.plot(df2; 
    x=Symbol(L"\eta_E-\eta_P"),
    y=Symbol("Runtime (s)"),
    color=Symbol("Injector Setup"),
    symbol=Symbol("Line Model"),
    mode="markers"
)

savehtmlplot(p, "runtime")