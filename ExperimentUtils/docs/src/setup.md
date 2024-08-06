```@setup 1
using Pkg; Pkg.activate("../..")
using PowerSystems
using PowerSimulationsDynamics
using Downloads
using ZIPE_loads
using TLmodels
using ExperimentUtils
```
# Setup

The best way to see how to use this package is a detailed example.
We'll start with how to set up the simulations you want to run.

## Constructor
In order to create a `GridSearchSys`, we can use the constructor. This creates all configurations of the given system and injectors, or specific given configurations.

```@docs; canonical=false
GridSearchSys(::System, ::Union{AbstractArray{DynamicInjection}, AbstractArray{DynamicInjection, 2}}, ::Union{AbstractArray{Vector{String}}, AbstractArray{String}, Nothing})
```

!!! note "Future Changes Possible"
    Admittedly, this is quite unclean. Ideally, injector permutations are separate from construction and are performed through the same interface as all the other sweeps. However, this was the way I began, so it's been stuck here. Hopefully in the future this can change.


Let's start with the IEEE 9-bus test case:

![IEEE WSCC 9-Bus Test Case](https://icseg.iti.illinois.edu/files/2013/10/WSCC9.png)

We need it to be all static components. We can use a .raw file to load it to a `PowerSystems.System` object like this:

```@example 1
system_raw_file = Downloads.download("https://raw.githubusercontent.com/gecolonr/loads/main/data/raw_data/WSCC_9bus.raw");
mv(system_raw_file, system_raw_file*".raw")
system_raw_file *= ".raw"
sys = System(system_raw_file, time_series_in_memory=true);
rm(system_raw_file)
sys
```

I'd also recommend you do this, as it helps remove clutter from the output.
```@example 1
set_runchecks!(sys, false);
```

Using `time_series_in_memory=true` I *believe* makes serialization work better but I'm not sure that it's necessary.

Now, we want to choose our injectors. Let's create a dynamic inverter and a dynamic generator. The functions we use here are from 

```@example 1
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
```

Now we're ready to call the constructor! For this example, we'll just get all combinations of the two injectors at each of the busses.

!!! warn "grid following inverters"
    Grid following inverters CANNOT be attached to the reference bus. If you need to do this, you'll likely have to switch your reference frame to `ConstantFrequency`. 
    If you try to run simulations with a grid following inverter attached to the reference bus, you will get
    ```
    ┌ Error: ResidualModel failed to build
    │   exception =
    │    KeyError: key :ω_oc not found
    │    Stacktrace:
    │      [1] getindex
    ...
    ```

!!! details "implementation details about memory usage"
    The injector variations created by the constructor are the only things that are **not** done lazily. This means that for every new injector setup, the base system will be deepcopied.  
    Other variations (those added with `add_generic_sweep!` or similar) will not instantiate new systems, so the total number of systems held in memory at once stays small enough for most use cases.


```@example 1
gss = GridSearchSys(sys, [gfm_inj(), sm_inj()])
```

Here are some possible alternatives that showcase the full behavior of the constructor:

#### Two specific test cases
Here, notice that the `injectors` array is now two dimensional.
We'll try these combinations of injectors:

| Bus 1 | Bus 2 | Bus 3 |
| ----- | ----- | ----- |
| GFM   | SM    | SM    |
| GFM   | GFM   | SM    |

```@example 1
GridSearchSys(sys, [gfm_inj() sm_inj() sm_inj()
                    gfm_inj() gfm_inj() sm_inj()],
                    ["bus1", "bus 2", "Bus 3"])
```

#### All combinations, with bus groups

Here, we group busses one and two together so they always have the same injector. Note that they aren't electrically connected; they simply always have identical machines attached.

```@example 1
#                                             bus group 1     bus group 2
GridSearchSys(sys, [gfm_inj(), sm_inj()], [["Bus1", "Bus 2"], ["Bus 3"]])
```

#### Bus groups and specific combinations
You can combine these two approaches, too. Let's try two specific test cases with those bus groups:

| Bus 1, Bus 2 | Bus 3 |
| ------------ | ----- |
| GFM          | SM    |
| SM           | SM    |

```@example 1
GridSearchSys(sys, [gfm_inj() sm_inj()
                    sm_inj()  sm_inj()],
        [["Bus1", "Bus 2"], ["Bus 3"]])
```

!!! note "a few notes"
    - if passing specific combinations, it is recommended to also pass `busgroups` to explicitly set the order of the busses. In the second example, `busgroups` is technically redundant, but it makes it clear what injector is being assigned to what bus.
    - When using `busgroups`, the strings are the bus names. Make sure they match precisely! You'll notice here that `"Bus1"` is missing a space. That's because that is the actual name of the bus.



## Sweep Config

Now that we have our system all set up, we can add sweeps of other parameters. This happens with the `add_generic_sweep!` method.

```@docs; canonical=false
add_generic_sweep!
```

Let's say we want to vary the load scale of the grid by scaling the power draw and supply setpoints of each generator and load. We can write a function that does this for one system:

```@example 1
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
```
Notice that although this function is marked as in-place and modifies `sys`, it still returns a system. This is important! The system that will actually be used in simulation is whatever is *returned*.

Now we can provide our function to `add_generic_sweep!` to apply for all systems:

```@example 1
load_scales_to_try = [0.4, 1.0]
add_generic_sweep!(gss, "Load Scale", set_power_setpt!, load_scales_to_try)
gss
```

This operation is incredibly lightweight; although this output seems to indicate that there are now 16 systems (8 initially times 2 load scale options), these have not actually been instantiated. This will happen lazily to ensure we don't run out of memory.

We can continue to add as many sweeps as we want. Let's add one for ZIPE load parameters.

```@example 1
ZIPE_params_to_try = (x->ZIPE_loads.LoadParams(x...)).([
#     Z    I    P    E
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    # [0.3, 0.3, 0.3, 0.1], # let's only run
    # [0.2, 0.2, 0.2, 0.4], # two sets of params
    # [0.1, 0.1, 0.1, 0.7], # so docs don't take
    # [0.0, 0.0, 0.0, 1.0], # years to build
])

function add_zipe_load(sys::System, params::ZIPE_loads.LoadParams)::System
    ZIPE_loads.create_ZIPE_load(sys, params)
    return sys
end

add_generic_sweep!(gss, "ZIPE Params", add_zipe_load, ZIPE_params_to_try)
gss
```

## Results

Results columns can be added before or after the simulations are run using the `add_result!` method. These are the things recorded for each simulation that is run.

```@docs; canonical=false
add_result!
```

One of these methods adds one result and the other adds multiple at a time, one for each element of the output of the `getter` function.

signature of thhe `getter` function is very important. Thankfully, there are many many examples! See [Results Getters](@ref) for all of the builtin getters you can try.

For now, let's just grab the eigenvalues of the system at the initial operating point and the transient current magnitude at each injector.

```@example 1
add_result!(gss, "Eigenvalues", eu.get_eigenvalues)
add_result!(gss, ["Bus 3 Injector Current", "Bus 1 Injector Current", "Bus 2 Injector Current"], eu.get_injector_currents)
```
Note that the vector of titles has a strange order - this is specific to match the order that `get_injector_currents` returns!