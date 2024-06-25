#%%
import pandas as pd
import numpy as np
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
#%%
df = pd.read_pickle("../data/fineresults_powersetpt_python.pkl")
# %%
df["Injector Setup"] = list(map(
    lambda x: ",".join(x), 
    zip(df["injector at {Bus1}"], 
        df["injector at {Bus 2}"], 
        df["injector at {Bus 3}"])
))
df.drop(["injector at {Bus1}", "injector at {Bus 2}", "injector at {Bus 3}"], axis=1, inplace=True)
# %%
df["ZIPE Load Params"] = list(map(lambda x: tuple(x.tolist()), df["ZIPE Load Params"]))
def categorize(name): df[name] = df[name].astype("category")
_ = list(map(categorize, ["Power Setpoint", "Line Model", "ZIPE Load Params", "Injector Setup"]))
# %%

df["initial_eigs_real"] = list(map(np.real, df["initial_eigenvalues"]))
df["initial_eigs_imag"] = list(map(np.imag, df["initial_eigenvalues"]))
df["final_eigs_real"] = list(map(np.real, df["final_eigenvalues"]))
df["final_eigs_imag"] = list(map(np.imag, df["final_eigenvalues"]))
df.drop(["initial_eigenvalues", "final_eigenvalues"], axis=1, inplace=True)
#%%
df_eigs = df.explode(["initial_eigs_real", "initial_eigs_imag", "final_eigs_real", "final_eigs_imag"], ignore_index=True)
df_eigs.drop(["Bus 1 Inverter Current", "Bus 2 Inverter Current", "Bus 3 Inverter Current", "Load Voltage at Bus 5", "Load Voltage at Bus 6", "Load Voltage at Bus 8", "initial_op_pt", "final_op_pt"], axis=1, inplace=True)

df_eigs["Index"] = sum([list(range(i)) for i in map(len, df["initial_eigs_real"])], start=[])

df_eigs_melted = df_eigs.melt(
    id_vars=["Injector Setup", "ZIPE Load Params", "Power Setpoint", "Line Model", "Index"],
    value_vars=["initial_eigs_real", "final_eigs_real"],
    var_name="initialorfinal",
    value_name="real"
)
df_eigs_melted["imag"] = df_eigs.melt(value_vars=["initial_eigs_imag", "final_eigs_imag"], value_name="imag")["imag"]
df_eigs_melted["initialorfinal"]=list(map(lambda x: {"f":"final", "i":"initial"}[x[0]], df_eigs_melted["initialorfinal"]))
df_eigs_melted.sort_values("Power Setpoint", ascending=True, inplace=True)
# df_eigs_melted["index"]=list(range(len(df_eigs_melted)))
# %%
fig = px.scatter(df_eigs_melted, 
           x="real",
           y="imag",
           color="ZIPE Load Params",
           facet_col="Injector Setup",
           facet_row="Line Model",
           animation_frame="Power Setpoint",
           symbol="initialorfinal",
        #    hover_data="Index"
           labels={
               "real":r"$\mathrm{Re}[\lambda]$", 
               "imag":r"$\mathrm{Im}[\lambda]$",
            },
            title="System Eigenvalues Before and After Branch Trip",
            size_max=7)
fig["layout"].pop("updatemenus") # drop animation buttons
fig
# %%
fig.write_html("../media/plotly_plot.html", auto_play=False, include_plotlyjs="cdn", include_mathjax="cdn")
# %%
df_trans = df.drop([
       'Bus 1 Inverter Current', 'Bus 2 Inverter Current',
       'Bus 3 Inverter Current', #'Load Voltage at Bus 5',
       'Load Voltage at Bus 6', 'Load Voltage at Bus 8', 'initial_eigenvalues',
       'initial_op_pt', 'final_eigenvalues', 'final_op_pt'
       ], axis=1).explode(["time", "Load Voltage at Bus 5"])
df_trans["Load Voltage at Bus 5"] = np.round(df_trans["Load Voltage at Bus 5"].astype(np.float32), decimals=5)
# %%
fig = px.scatter(df_trans,
           x="time",
           y="Load Voltage at Bus 5",
           color="ZIPE Load Params",
           facet_row="Line Model",
           facet_col="Injector Setup",
           animation_frame="Power Setpoint",
           title="Transient Simulation: Load Voltage at Bus 5")
fig["layout"].pop("updatemenus") # drop animation button
# %%
fig.write_html("../media/plotly_transient_plot.html", auto_play=False, include_plotlyjs="cdn", include_mathjax="cdn")
# %%
