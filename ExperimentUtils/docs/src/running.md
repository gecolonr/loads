# Running Simulations

```@setup 1
# this is code from the previous page
using Pkg; Pkg.activate("../..")
using PowerSystems
using PowerSimulationsDynamics
using Downloads
using ZIPE_loads
using TLmodels
using ExperimentUtils
import Logging
Logging.disable_logging(Logging.Error)
system_raw_file = Downloads.download("https://raw.githubusercontent.com/gecolonr/loads/main/data/raw_data/WSCC_9bus.raw");
mv(system_raw_file, system_raw_file*".raw")
system_raw_file *= ".raw"
sys = System(system_raw_file, time_series_in_memory=true);
set_runchecks!(sys, false)
rm(system_raw_file)
const eu = ExperimentUtils
gfm_inj() = DynamicInverter(
    "GFM", # Grid forming control
    1.0, # ω_ref,
    eu.converter_high_power(), #converter
    eu.VSM_outer_control(), #outer control
    eu.GFM_inner_control(), #inner control voltage source
    eu.dc_source_lv(), #dc source
    eu.pll(), #pll
    eu.filt(), #filter
)
sm_inj() = DynamicGenerator(
    "SM", # Synchronous Machine
    1.0, # ω_ref,
    eu.AF_machine(), #machine
    eu.shaft_no_damping(), #shaft
    eu.avr_type1(), #avr
    eu.tg_none(), #tg
    eu.pss_none(), #pss
)
gss = GridSearchSys(sys, [gfm_inj(), sm_inj()])

function set_power_setpt!(sys::System, scale::Real)
    for load in get_components(StandardLoad, sys)
        set_impedance_active_power!(load, get_impedance_active_power(load)*scale)
        set_current_active_power!(load, get_current_active_power(load)*scale)
        set_constant_active_power!(load, get_constant_active_power(load)*scale)
        
        set_impedance_reactive_power!(load, get_impedance_reactive_power(load)*scale)
        set_current_reactive_power!(load, get_current_reactive_power(load)*scale)
        set_constant_reactive_power!(load, get_constant_reactive_power(load)*scale)
    end
    for gen in get_components(Generator, sys)
        if gen.bus.bustype == ACBusTypes.PV
            set_active_power!(gen, get_active_power(gen) * scale)
            set_reactive_power!(gen, get_reactive_power(gen) * scale)
        end
    end
    return sys
end

load_scales_to_try = [0.4, 1.0]
add_generic_sweep!(gss, "Load Scale", set_power_setpt!, load_scales_to_try)

ZIPE_params_to_try = (x->ZIPE_loads.LoadParams(x...)).([
#     Z    I    P    E
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    # [0.3, 0.3, 0.3, 0.1],
    # [0.2, 0.2, 0.2, 0.4],
    # [0.1, 0.1, 0.1, 0.7],
    # [0.0, 0.0, 0.0, 1.0],
])

function add_zipe_load(sys::System, params::ZIPE_loads.LoadParams)
    ZIPE_loads.create_ZIPE_load(sys, params)
    return sys
end

add_generic_sweep!(gss, "ZIPE Params", add_zipe_load, ZIPE_params_to_try)

add_result!(gss, "Eigenvalues", eu.get_eigenvalues)
add_result!(gss, ["Bus 3 Injector Current", "Bus 1 Injector Current", "Bus 2 Injector Current"], eu.get_injector_currents)

```

Now that we have everything set up, we can run our simulations! This is fairly simple with the [`execute_sims!`](@ref) method. The stuff about saving to file and `chunksize` will be explained on the [next page](saving.md)

```@docs; canonical=false
execute_sims!
```

For our example, let's run short transient simulations with a BranchTrip:

```@example 1
execute_sims!(
    gss, 
    BranchTrip(0.5, ACBranch, "Bus 5-Bus 4-i_1"), 
    tspan=(0.49, 0.55), 
    dtmax=0.05, 
    run_transient=true, 
    log_path="example_sims"
)
```

Now we can inspect the results!
```@example 1
gss
```
```@repl 1
gss.df
```

## Notes
 - You must only use one perturbation if you want to save the results to file. For some reason, passing a `Vector{Perturbation}` changes something about the typing and breaks serialization.
 - if your perturbation does not occur at $$t=0.5$$, make sure to pass the `tstops` argument and include at least the time of your perturbation. This significantly helps numerical stability.
 - `execute_sims!` is fully parallelized, so use as many cores as you can! sometimes the REPL is limited to one thread, so I've found more success running sims from the command line and passing an explicit thread count with `julia -t 69 my_experiments.jl`
    - be sure to mind your memory usage - these systems are big and therefore solving them can be costly.
    - use `tmux` so you don't have to sit there all day remoted into the server.
 - don't spend a ton of time trying to tune `ida_opts`. It doesn't help much.
 - mind the ratio of the time span to `dtmax`. The data accumulates quickly, and RAM is limited.