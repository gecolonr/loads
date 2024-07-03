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

df = load_serde_data("data/forplot_fine").df
df = @subset df begin
    @byrow getvec(:"ZIPE Load Params")[1] != 1.0
end
df[!, "ZIPE Load Params"] = getvec.(df[!, "ZIPE Load Params"])
df[!, "markershape"] = map(x-> if x[3] > 0 L"\large{\eta_P > 0}" elseif x[4] > 0 L"\large{\eta_E > 0}" else "E=P=0" end, df[!, "ZIPE Load Params"])
df[!, "legendgroup"] = map(x-> if x[3] > 0 L"{}" elseif x[4] > 0 L"{}{}" else L"{{}}" end, df[!, "ZIPE Load Params"])
df[!, "Injector Setup"] = join.(eachrow(hcat(df[!, "injector at {Bus1}"], 
df[!, "injector at {Bus 2}"],
df[!, "injector at {Bus 3}"])), ", ")
df[!, "colorbar"] = map(x->x[3]-x[4], df[!, "ZIPE Load Params"])
df[!, "real"] = real.(df.Eigenvalues)
df[!, "imag"] = imag.(df.Eigenvalues)
df[!, "hovertext"] = "η=" .* string.(df[!, "ZIPE Load Params"])
df[!, "e_percent"] = map(last, df[!, "ZIPE Load Params"])
df[!, "transient_legendgroup_hack"] = map(x->LaTeXString("\$\\iffalse $(x) \\fi\$"), df[!, "hovertext"])
df[!, "power_percent"] = map(x->(1.0-2*first(x)), df[!, "ZIPE Load Params"])
function pfactors(sm)
    pf = summary_participation_factors(sm)
    return Dict(pf.Name .=> eachrow(select(pf, Not(:Name))))
end
df[!, "participation_factors"] = map(sm->get(pfactors(sm), "generator-3-1 ir_cnv", missing), df.sm)
df[!, "participation_factors"] = map(x->x isa Missing ? "" : "\nPF: ".*string.(round.((collect(x)), sigdigits=5)), df.participation_factors)
df[!, "hovertext_eig"] = [i.*j for (i, j) in zip(df[!, "hovertext"], df[!, "participation_factors"])]
df = sort(df, ["Power Setpoint", "Line Model", "colorbar", "hovertext", "e_percent", "power_percent"])

p = makeplots(
    df;
    rows="Line Model",
    # cols="Injector Setup",
    color="colorbar",
    legendgroup="transient_legendgroup_hack",
    legendgroup_title_func=x->nothing,
    # markershape="markershape",
    trace_names="hovertext",
    colorbar=false,
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
    y="Bus 2 Injector Current",

    x_title=L"\large{\mathrm{Time}\:\: [\mathrm{s}]}",
    y_title = L"\large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",
    # row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
    row_title_func = x->(x*" lines"),

    supertitle=L"\large{\text{Transient Current at Bus 2}}",
    image_export_filename="transient_plot",
    colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11))
)
savehtmlplot(p, "transient_niceified_fine_bus2")

p = makeplots(
    df;
    rows="Line Model",
    # cols="Injector Setup",
    color="colorbar",
    legendgroup="hovertext",
    # markershape="markershape",
    trace_names="Line Model",
    colorbar=true,
    slider="Power Setpoint",
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

    x_title=L"\large{\mathrm{Re}(\lambda)}",
    y_title=L"\large{\mathrm{Im}(\lambda)}",
    # row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
    row_title_func=x->(x*" lines"),

    supertitle=L"\large{\text{System Eigenvalues}}",
    image_export_filename="eigenvalue_plot",
    
    colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11)),
    legend_visible=false,
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax=1.0,
        # cmid=0.5,
        cmin=-1.0,
        title="",
        colorbar=attr(
            # bordercolor=colorant"black",
            # borderwidth=1,
            # title=attr(text="Load Scale"),
            outlinecolor=colorant"black",
            outlinewidth=1,
            tickmode="array",
            ticktext=reverse(["P=1.0", "P=0.8", "P=0.6", "P=0.4", "P=0.2", "P=E=0", "E=0.2", "E=0.4", "E=0.6", "E=0.8", "E=1.0"]),
            tickvals=collect(-1.0:0.2:1.0),
            xref="paper",
            yref="paper",
            x=1.02
            # xpad=20,
            
        ),
        colorscale=[[-1.0, RGB(colorant"red")], [1.0, RGB(colorant"blue")]],
    )),
)
savehtmlplot(p, "eigplot_niceified_fine")


p = makeplots(
    @subset df begin
        :"power_percent" .> 0.0
        #* UNCOMMENT NEXT LINE TO EXCLUDE HIGH LOAD SCALES
        # :"Power Setpoint" .< 0.37
    end;
    rows="Line Model",
    # cols="markershape",
    color="Power Setpoint",
    legendgroup="legendgroup",
    markershape="markershape",
    trace_names="markershape",
    colorbar=true,
    slider="power_percent",
    hovertext="hovertext_eig",
    # opacity=0.7,
    scattermode="markers",
    markersize=10,
    # color_sort_func=identity,

    slider_current_value_prefix="P/E Fraction: ",

    x="real",
    y="imag",
    # x="time",
    # y="Bus 3 Injector Current",

    x_title=L"\large{\mathrm{Re}(\lambda)}",
    y_title=L"\large{\mathrm{Im}(\lambda)}",
    row_title_func=x->(x*" lines"),
    # row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),

    supertitle=L"\large{\text{System Eigenvalues}}",
    image_export_filename="eigenvalue_plot",
    colorlist = "#" .* hex.(range(colorant"red", colorant"blue", length=11)),
    legend_location=(x=0.97, y=1.02),
    colorbar_args=Dict(attr(
        autocolorscale=false, 
        cmax=1.0,
        # cmid=0.5,
        cmin=0.05,
        # title="Load Scale",
        colorbar=attr(
            # bordercolor=colorant"black",
            # borderwidth=1,
            title=attr(text="Load Scale"),
            outlinecolor=colorant"black",
            outlinewidth=1,
            # tickmode="array",
            # ticktext=["P=1.0", "P=0.8", "P=0.6", "P=0.4", "P=0.2", "P=E=0", "E=0.2", "E=0.4", "E=0.6", "E=0.8", "E=1.0"],
            # tickvals=collect(-1.0:0.2:1.0),
            xref="paper",
            yref="paper",
            x=1.02
            # xpad=20,
            
        ),
        colorscale=[[-1.0, RGB(colorant"red")], [1.0, RGB(colorant"blue")]],
    )),
)
savehtmlplot(p, "eigplot_vsloadscale_fine")
