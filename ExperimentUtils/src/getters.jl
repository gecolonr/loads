"""
WARNING: mostly deprecated. just use `get_sim`.

Saves the `sim` object to file, then returns the filename to be added to the results dataframe.
"""
function get_serialized_sim_filename_function_builder(path::String)
    try
        mkdir(path)
    catch
        println("Failed to create directory $path (maybe it already exists)")
    end

    function get_serialized_sim_filename(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing})
        if sim isa Missing
            return missing
        end
        filename = "$(path)/SimObject$(time_ns()).jls"
        Serialization.serialize(filename, sim)
        return [filename]
    end
    return get_serialized_sim_filename
end

"""
gets the time taken to run this simulation.
"""
function get_dt(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return dt
end

"""
gets system eigenvalues from small signal analysis.
"""
function get_eigenvalues(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return sm isa Missing ? missing : [sm.eigenvalues]
end


"""
gets system eigenvectors from small signal analysis.
"""
function get_eigenvectors(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return sm isa Missing ? missing : [sm.eigenvectors]
end

"""
gets the small signal analysis object `sm`
"""
function get_sm(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return sm
end

"""
**RETURNS A VECTOR!** use a vector of column titles if you want each bus to be a separate column.

gets voltage time series at every bus in the system.

returns them in the same order as `get_components(Generator, gss.base)`.
"""
function get_bus_voltages(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Bus, gss.base)))
    else
        return [get_voltage_magnitude_series(sim.results, i)[2] for i in get_number.(get_components(Bus, gss.base))]
    end
end

# """
# extract array of current magnitude over time at the object with name `name` from the simulation results `res`.

# Timestamps are thrown away. If you're using the GridSearchSystem, the returned array will correspond to an even time grid with spacing `dtmax` (the argument you passed to [`execute_sims`](@ref), default 0.02s).
# """
function current_magnitude(res::SimulationResults, name::String)
    _, ir = PSID.post_proc_real_current_series(res, name, nothing)
    _, ii = PSID.post_proc_imaginary_current_series(res, name, nothing)
    return sqrt.(ir.^2 .+ ii.^2)
end

"""
**RETURNS A VECTOR!** use a vector of column titles if you want each inverter to be a separate column.

Returns a vector with an entry for each `Generator` in the base system. Each entry is a time series of the current magnitude.

returns them in the same order as `get_components(Generator, gss.base)`.
"""
function get_injector_currents(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    gens = get_components(Generator, gss.base)
    if sim isa Missing
        return Array{Missing}(missing, length(gens))
    end

    gen_dict = Dict(get_name.(gens) .=> get_number.(get_bus.(gens)))
    # get dynamic injectors. filter by those in gen_dict (ie, in gss.base) to exclude ZIPE inverters
    injectors = [i for i in PSID.get_dynamic_injectors(sim.inputs) if get_name(i) in keys(gen_dict)]
    injectors = Dict(map(x->gen_dict[x], get_name.(injectors)) .=> injectors)

    return [(i in keys(injectors) ? current_magnitude(read_results(sim), get_name(injectors[i])) : missing) for i in get_number.(get_bus.(gens))]
end

"""
**RETURNS A VECTOR!** use a vector of column titles if you want each inverter to be a separate column.

Returns a vector with an entry for each `Generator` in the base system. If there is an inverter there, the entry is a time series of the current magnitude. Otherwise, it's `missing`.

returns them in the same order as `get_components(Generator, gss.base)`.
"""
function get_inverter_currents(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    gens = get_components(Generator, gss.base)
    if sim isa Missing
        return Array{Missing}(missing, length(gens))
    end
    gen_dict = Dict(get_name.(gens) .=> get_number.(get_bus.(gens)))
    # get dynamic inverters. filter by those in gen_dict (ie, in gss.base) to exclude ZIPE inverters
    inverters = [i for i in get_components(DynamicInverter, sim.sys) if get_name(i) in keys(gen_dict)]
    inverters = Dict(map(x->gen_dict[x], get_name.(inverters)) .=> inverters)
    return [(i in keys(inverters) ? current_magnitude(read_results(sim), get_name(inverters[i])) : missing) for i in get_number.(get_bus.(gens))]
end


"""
**RETURNS A VECTOR!** use a vector of column titles if you want each generator to be a separate column.

gets speed (in rad/s) of all generators in the system.

returns them in the same order as `get_components(Generator, gss.base)`.
"""
function get_generator_speeds(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    throw("unimplemented: not tested yet")
    if sim isa Missing
        return Array{Missing}(missing, length(get_components(Generator, gss.base)))
    end
    gen_busses = collect(get_bus.(get_components(Generator, gss.base)))
    gen_dict = Dict(get_name.(get_components(Generator, gss.base)) .=> get_number.(gen_busses))
    generators = [i for i in get_components(DynamicGenerator, sys) if i.name in keys(gen_dict)]
    generators = Dict(map(x->gen_dict[x], get_name.(generators)) .=> generators)
    return [(i in keys(generators) ? generator_speed(res, generators[i].name) : missing) for i in get_number.(gen_busses)]
end

"""
**RETURNS A VECTOR!** use a vector of column titles if you want each load to be a separate column.

Gets the voltage magnitude time series at all ZIPE loads. 

returns them in the same order as `get_components(Generator, gss.base)`.
"""
function get_zipe_load_voltages(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    if (sim isa Missing) || (sim.results isa Nothing)
        return Array{Missing}(missing, length(get_components(StandardLoad, gss.base)))
    end
    # the [2] here just gets the voltage instead of a tuple of (time, voltage) since we know the time divisions
    return [get_voltage_magnitude_series(sim.results, i)[2] for i in get_number.(get_bus.(get_components(StandardLoad, gss.base)))]
end

"""
Gets the timestamps for transient results. 
"""
function get_time(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    if sim isa Missing
        return [missing]
    end
    return [get_voltage_magnitude_series(sim.results, get_number(get_bus(first(get_components(StandardLoad, gss.base)))))[1]]
end
"""
**RETURNS A VECTOR!** use a vector of column titles if you want each inverter to be a separate column.

Gets the currrent magnitude time series at each ZIPE load.
"""
function get_zipe_load_current_magnitudes(gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    throw("unimplemented: not tested yet.")
    if (sim isa Missing) || (sim.results isa Nothing)
        return Array{Missing}(missing, length(get_components(StandardLoad, gss.base)))
    end
    
    return ([current_magnitude(sim.results, "load_GFL_inverter"*string(get_number(get_bus(load)))).+ current_magnitude(sim.results, load.name) for load in get_components(StandardLoad, gss.base)])
end

"""
gets the value `sim.status` from the Simulation object (or missing).
"""
function get_sim_status(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return (sim isa Missing) ? missing : sim.status
end

"""
gets the string representation of the error raised during simulation (or missing)
"""
function get_error(_gss::GridSearchSys, _sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, error::Union{String, Missing}, dt::Real)
    return error
end

"""
gets the whole Simulation object.
"""
function get_sim(_gss::GridSearchSys, sim::Union{Simulation, Missing}, _sm::Union{PSID.SmallSignalOutput, Missing}, _error::Union{String, Missing}, dt::Real)
    return sim
end
