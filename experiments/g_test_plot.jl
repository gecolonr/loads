path = "data/gab_tests/results0.jls"

df = load_serde_data(path::String)

p = makeplots(
    df;
	
    rows="Line Model",
    # cols="Injector Setup",
    color="ZIPE Load Params",
    slider="Power Setpoint",
    slider_current_value_prefix="Power Setpoint: ",
	
    x=(x0=0.48, dx=0.00005),
    y="Bus 3 Inj Current",
    x_title=L"\mathrm{Time}\:\: [\mathrm{s}]",
    y_title = L"\mathrm{Current}\:\:[\mathrm{p.u.}]",
	
    supertitle="Bus 3 Injector Current",
    yaxis_home_range = (min=0, max=10),
)