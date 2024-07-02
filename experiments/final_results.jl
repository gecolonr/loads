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
using Colors

include("sysbuilder.jl")
include("plotting.jl")

getvec(x::LoadParams) = round.([x.z_percent, x.i_percent, x.p_percent, x.e_percent], digits=3)

function make_plots_from_paper(;save_to_html=true, gss=nothing)
    if isnothing(gss)
        gss = load_serde_data("data/forplot")
    end
    map!(getvec, gss.df[!, "ZIPE Load Params"], gss.df[!, "ZIPE Load Params"])
    gss.df[!, "hovertext"] = map(x->"η=$x", gss.df[!, "ZIPE Load Params"])
    gss.df[!, "real_eigs"] = map(real, gss.df[!, "Eigenvalues"])
    gss.df[!, "imag_eigs"] = map(imag, gss.df[!, "Eigenvalues"])
    function pfactors(sm)
        pf = summary_participation_factors(sm)
        if "generator-3-1 ir_filter" ∉ pf.Name
            return missing
        end
        return ", PF: " .* string.(
            round.(
                collect(
                    first(select(
                        pf[pf.Name .== "generator-3-1 ir_filter", :], 
                        Not(:Name)
                    ))
                ), sigdigits=5
            )
        )
    end
    gss.df[!, "pfactors"] = map(pfactors, gss.df[!, "sm"])
    gss.df[!, "hovertext_eig"] = reduce.(.*, eachrow(gss.df[!, ["hovertext", "pfactors"]]))
    gss.df[!, "colorbar"] = map(x->(x[4]-x[3]), gss.df[!, "ZIPE Load Params"])
    eigplot_slider = makeplots(
        gss.df;
        rows="Line Model",
        color="hovertext", #ZIPE Load Params
        trace_names="Line Model",
        hovertext="hovertext_eig",

        slider="Power Setpoint",
        slider_current_value_prefix="Power Setpoint: ",

        x="real_eigs",
        y="imag_eigs",
        x_title=L"\large{\mathrm{Re}(\lambda)}",
        y_title=L"\large{\mathrm{Im}(\lambda)}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->-Dict(
            "η=[1.0, 0.0, 0.0, 0.0]" => 1,
            "η=[0.0, 0.0, 1.0, 0.0]" => 2,
            "η=[0.1, 0.1, 0.8, 0.0]" => 3,
            "η=[0.2, 0.2, 0.6, 0.0]" => 4,
            "η=[0.3, 0.3, 0.4, 0.0]" => 5,
            "η=[0.4, 0.4, 0.2, 0.0]" => 6,
            "η=[0.5, 0.5, 0.0, 0.0]" => 7,
            "η=[0.4, 0.4, 0.0, 0.2]" => 8,
            "η=[0.3, 0.3, 0.0, 0.4]" => 9,
            "η=[0.2, 0.2, 0.0, 0.6]" => 10,
            "η=[0.1, 0.1, 0.0, 0.8]" => 11,
            "η=[0.0, 0.0, 0.0, 1.0]" => 12,
         )[x],
        supertitle=L"\large{\text{System Eigenvalues}}",

        image_export_filename="eigenvalue_plot",
        scatterplot_args=Dict(:mode=>"markers"),
        colorlist=reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11)))),
    )
    eigplot_noslider = makeplots(
        @subset gss.df :"Power Setpoint" .≈ 0.2;
        rows="Line Model",
        color="hovertext", # ZIPE Load Params
        colorbar="colorbar",
        trace_names="Line Model",
        hovertext="hovertext_eig",

        x="real_eigs",
        y="imag_eigs",
        x_title=L"\large{\mathrm{Re}(\lambda)}",
        y_title=L"\large{\mathrm{Im}(\lambda)}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->-Dict(
            "η=[1.0, 0.0, 0.0, 0.0]" => 1,
            "η=[0.0, 0.0, 1.0, 0.0]" => 2,
            "η=[0.1, 0.1, 0.8, 0.0]" => 3,
            "η=[0.2, 0.2, 0.6, 0.0]" => 4,
            "η=[0.3, 0.3, 0.4, 0.0]" => 5,
            "η=[0.4, 0.4, 0.2, 0.0]" => 6,
            "η=[0.5, 0.5, 0.0, 0.0]" => 7,
            "η=[0.4, 0.4, 0.0, 0.2]" => 8,
            "η=[0.3, 0.3, 0.0, 0.4]" => 9,
            "η=[0.2, 0.2, 0.0, 0.6]" => 10,
            "η=[0.1, 0.1, 0.0, 0.8]" => 11,
            "η=[0.0, 0.0, 0.0, 1.0]" => 12,
         )[x],
        supertitle=L"\large{\text{System Eigenvalues}}",

        image_export_filename="eigenvalue_plot",
        scatterplot_args=Dict(:mode=>"markers"),
        # colorlist=reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11)))),
    )
    transient_plot_slider = makeplots(
        gss.df;
        rows="Line Model",
        color="hovertext", # ZIPE Load Params
        trace_names="Line Model",
        hovertext="hovertext",

        slider="Power Setpoint",
        slider_current_value_prefix="Power Setpoint: ",

        x="time",
        y="Bus 3 Injector Current",
        x_title=L"\large{\mathrm{Time}\:\: [\mathrm{s}]}",
        y_title = L"\large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->-Dict(
            "η=[1.0, 0.0, 0.0, 0.0]" => 1,
            "η=[0.0, 0.0, 1.0, 0.0]" => 2,
            "η=[0.1, 0.1, 0.8, 0.0]" => 3,
            "η=[0.2, 0.2, 0.6, 0.0]" => 4,
            "η=[0.3, 0.3, 0.4, 0.0]" => 5,
            "η=[0.4, 0.4, 0.2, 0.0]" => 6,
            "η=[0.5, 0.5, 0.0, 0.0]" => 7,
            "η=[0.4, 0.4, 0.0, 0.2]" => 8,
            "η=[0.3, 0.3, 0.0, 0.4]" => 9,
            "η=[0.2, 0.2, 0.0, 0.6]" => 10,
            "η=[0.1, 0.1, 0.0, 0.8]" => 11,
            "η=[0.0, 0.0, 0.0, 1.0]" => 12,
         )[x],
        supertitle=L"\large{\text{Bus 3 (GFL) Current Magnitude}}",

        yaxis_home_range = (min=0.0, max=2.0),
        xaxis_home_range = (min=0.48, max=1.0),

        image_export_filename = "transient_current_plot",
        colorlist = reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))))
    )
    transient_plot_noslider = makeplots(
        @subset gss.df :"Power Setpoint" .≈ 0.2;
        rows="Line Model",
        color="hovertext", # ZIPE Load Params
        trace_names="Line Model",
        hovertext="hovertext",

        x="time",
        y="Bus 3 Injector Current",
        x_title=L"\large{\mathrm{Time}\:\: [\mathrm{s}]}",
        y_title = L"\large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->-Dict(
            "η=[1.0, 0.0, 0.0, 0.0]" => 1,
            "η=[0.0, 0.0, 1.0, 0.0]" => 2,
            "η=[0.1, 0.1, 0.8, 0.0]" => 3,
            "η=[0.2, 0.2, 0.6, 0.0]" => 4,
            "η=[0.3, 0.3, 0.4, 0.0]" => 5,
            "η=[0.4, 0.4, 0.2, 0.0]" => 6,
            "η=[0.5, 0.5, 0.0, 0.0]" => 7,
            "η=[0.4, 0.4, 0.0, 0.2]" => 8,
            "η=[0.3, 0.3, 0.0, 0.4]" => 9,
            "η=[0.2, 0.2, 0.0, 0.6]" => 10,
            "η=[0.1, 0.1, 0.0, 0.8]" => 11,
            "η=[0.0, 0.0, 0.0, 1.0]" => 12,
         )[x],
        supertitle=L"\large{\text{Bus 3 (GFL) Current Magnitude}}",

        yaxis_home_range = (min=0.317, max=0.365),
        xaxis_home_range = (min=0.498, max=0.53),

        image_export_filename = "transient_current_plot",
        colorlist = reverse(vcat("#00FF00", "#" .* hex.(range(colorant"red", colorant"blue", length=11))))
    )
    if save_to_html
        savehtmlplot(eigplot_slider, "eigplot_slider_from_paper")
        savehtmlplot(eigplot_noslider, "eigplot_noslider_from_paper")
        savehtmlplot(transient_plot_slider, "transient_plot_slider_from_paper")
        savehtmlplot(transient_plot_noslider, "transient_plot_noslider_from_paper")
    end
    return (eigplot_slider, eigplot_noslider, transient_plot_slider, transient_plot_noslider)
end

function make_plots_more_data(;save_to_html=true)
    gss = load_serde_data("data/fineresults_powersetpt")
    add_result!(gss, ["Bus 3 Injector Current", "Bus 1 Injector Current", "Bus 2 Injector Current"], get_injector_currents)
    inj_case_names = Dict(
        "SM, GFL, GFM" => "Case 1",
        "SM, GFM, GFL" => "Case 2", 
        "SM, GFM, GFL" => "Case 2", 
        "GFM, GFL, SM" => "Case 3",
        "SM, SM, SM" => "Case 4",
        "GFM, GFM, SM" => "Case 5",
        "GFM, SM, SM" => "Case 6",
    )
    gss.df[!, "Injector Setup"] = map(
        x->inj_case_names[x]*" ($x)",  
        join.(eachrow(hcat(gss.df[!, "injector at {Bus1}"], 
                           gss.df[!, "injector at {Bus 2}"], 
                           gss.df[!, "injector at {Bus 3}"])), ", ")
    )
    gss.df[!, "ZIPE Load Params"] = map(getvec, gss.df[!, "ZIPE Load Params"])
    gss.df[!, "hovertext"] = map(x->"η=$x", gss.df[!, "ZIPE Load Params"])
    gss.df[!, "real_eigs"] = map(real, gss.df[!, "Eigenvalues"])
    gss.df[!, "imag_eigs"] = map(imag, gss.df[!, "Eigenvalues"])
    gss.df[!, "tracename"] = map(x->first(x[1], 6)*", "*x[2], zip(gss.df[!, "Injector Setup"], gss.df[!, "Line Model"]))

    function pfactors(sm)
        pf = summary_participation_factors(sm)
        if "generator-3-1 ir_filter" ∉ pf.Name
            return missing # if this is a synchronous machine, output current isn't a dynamic state, so we can't get the participation factor.
        end
        return ", PF: " .* string.(
            round.(
                collect(
                    first(select(
                        pf[pf.Name .== "generator-3-1 ir_filter", :], 
                        Not(:Name)
                    ))
                ), sigdigits=5
            )
        )
    end
    gss.df[!, "pfactors"] = map(pfactors, gss.df[!, "sm"])
    gss.df[!, "hovertext_eig"] = reduce.(.*, eachrow(gss.df[!, ["hovertext", "pfactors"]]))


    eigplot = makeplots(
        gss.df;
        rows="Line Model",
        cols="Injector Setup",
        color="hovertext", # ZIPE Load Params
        trace_names="tracename",
        hovertext="hovertext_eig",

        slider="Power Setpoint",
        slider_current_value_prefix="Power Setpoint: ",

        x="real_eigs",
        y="imag_eigs",
        x_title=L"\large{\mathrm{Re}(\lambda)}",
        y_title=L"\large{\mathrm{Im}(\lambda)}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->(parse(Float64,first(last(x, 4), 3)) - parse(Float64, first(last(x, 9), 3))),

        supertitle=L"\large{\text{System Eigenvalues}}",

        image_export_filename="eigenvalue_plot",
        scatterplot_args=Dict(:mode=>"markers"),
        colorlist=reverse("#" .* hex.(range(colorant"red", colorant"blue", length=21))),
    )
    transient_plot = makeplots(
        gss.df;
        rows="Line Model",
        cols="Injector Setup",
        color="hovertext", # ZIPE Load Params
        trace_names="tracename",
        hovertext="hovertext",

        slider="Power Setpoint",
        slider_current_value_prefix="Power Setpoint: ",

        x="time" in names(gss.df) ? "time" : (x0=0.48, dx=0.00005),
        y="Bus 3 Injector Current",
        x_title=L"\large{\mathrm{Time}\:\: [\mathrm{s}]}",
        y_title = L"\large{||i_\text{filter}||\:\:[\mathrm{p.u.}]\:\:\text{(Bus 3)}}",

        row_title_func=x->LaTeXString("\$\\large{\\text{$x lines}}\$"),
        color_sort_func=x->(parse(Float64,first(last(x, 4), 3)) - parse(Float64, first(last(x, 9), 3))),
        supertitle=L"\large{\text{Bus 3 (GFL) Current Magnitude}}",

        yaxis_home_range = (min=0.0, max=2.0),
        xaxis_home_range = (min=0.48, max=1.0),

        image_export_filename = "transient_current_plot",
        colorlist = reverse("#" .* hex.(range(colorant"red", colorant"blue", length=21)))
    )
    if save_to_html
        savehtmlplot(eigplot, "eigplot_more_data")
        savehtmlplot(transient_plot, "transient_plot_more_data")
    end
    return (eigplot, transient_plot)

end