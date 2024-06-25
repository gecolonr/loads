using DataFrames
using DataFramesMeta
using PyCall
pd = pyimport("pandas")
include("sysbuilder.jl")

gss = load_serde_data("data/fineresults_powersetpt")


function get_zipe_load_voltages_sm(gss::GridSearchSys, sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing})
    if isnothing(sim) return missing end
    sys = deepcopy(sim.sys)
    remove_component!(Line, sys, sim.perturbations[1].branch_name)
    newsim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 1.0), disable_timer_outputs=true)
    sm = small_signal_analysis(newsim)
    out = zeros(length(get_components(StandardLoad, gss.base)))
    for (idx, i) in enumerate(get_number.(get_bus.(get_components(StandardLoad, gss.base))))
        bus_lookup = PSID.get_bus_lookup(sim.results)
        bus_ix = bus_lookup[i]
        V_R = sm.operating_point[bus_ix]
        V_I = sm.operating_point[bus_ix+PSID.get_bus_count(sim.results)]
        out[idx] = sqrt(V_R^2 + V_I^2)
    end

    return out
end
getvec(x::LoadParams) = round.([x.z_percent, x.i_percent, x.p_percent, x.e_percent], digits=3)

add_result!(gss, ["Bus 3 Inverter Current", "Bus 1 Inverter Current", "Bus 2 Inverter Current"], get_injector_currents)
add_result!(gss, "posttripsm", small_signal_tripped)
add_result!(gss, ["Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages)
add_result!(gss, ["Final Load Voltage at $busname" for busname in get_name.(get_bus.(get_components(StandardLoad, gss.base)))], get_zipe_load_voltages_sm)

gss.df[!, "initial_eigenvalues"] = gss.df[!, "Eigenvalues"]
gss.df[!, "final_eigenvalues"] = map(x->x.eigenvalues, gss.df[!, "posttripsm"])
gss.df[!, "initial_op_pt"] = map(x->x.operating_point, gss.df[!, "sm"])
gss.df[!, "final_op_pt"] = map(x->x.operating_point, gss.df[!, "posttripsm"])

gss.df[!, "ZIPE Load Params"] = map(getvec, gss.df[!, "ZIPE Load Params"])



df = select(gss.df, [
    "injector at {Bus1}",
    "injector at {Bus 2}",
    "injector at {Bus 3}",
    "Bus 1 Inverter Current",
    "Bus 2 Inverter Current",
    "Bus 3 Inverter Current",
    "Load Voltage at Bus 5",
    "Load Voltage at Bus 6",
    "Load Voltage at Bus 8",
    "Power Setpoint",
    "Line Model",
    "ZIPE Load Params",
    "initial_eigenvalues",
    "initial_op_pt",
    "final_eigenvalues",
    "final_op_pt",
])

df[!, 1] = Vector{String}(df[!, 1])
df[!, 2] = Vector{String}(df[!, 2])
df[!, 3] = Vector{String}(df[!, 3])
df[!, 10] = Vector{Float64}(df[!, 10])
df[!, 11] = Vector{String}(df[!, 11])
df[!, 13] = Vector{Vector{ComplexF64}}(df[!, 13])

using Pandas
df_jlpd = Pandas.DataFrame(df)

df_pd = pd.DataFrame(df_jlpd)

df_pd.to_pickle("/data/reiddye/loads/fineresults_powersetpt_python.pkl")

