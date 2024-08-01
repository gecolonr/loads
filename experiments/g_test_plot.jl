using TLmodels
path = "data/gab_tests/results0.jls"
df = load_serde_data(path::String)

include("plotting.jl")

p = makeplots(
    df;
	
    rows="Line Model",
    # color="ZIPE Parameters",
    slider="Power Setpoint",
    slider_current_value_prefix="Power Setpoint: ",
	
    x=(x0=0.48, dx=0.00005),
    y="Load Voltage at Bus 5",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Current}\:\:[\mathrm{p.u.}]",
	
    supertitle="Bus 3 Injector Current",
    yaxis_home_range = (min=0, max=10),
)

savehtmlplot(p, "data/gab_tests/plot0.html")

x = df[!,"Load Voltage at Bus 5"]
using Plots
plot(t, x)