push!(LOAD_PATH,"../src/")
using Pkg; Pkg.activate("../..")
using PowerSystems
using PowerSimulationsDynamics


using Documenter, ExperimentUtils
DocMeta.setdocmeta!(ExperimentUtils, :DocTestSetup, :(using ExperimentUtils; using Pkg; Pkg.activate(".."); using PowerSystems; using PowerSimulationsDynamics); recursive=true)
makedocs(
    sitename="ExperimentUtils Documentation", 
    modules=[ExperimentUtils],
    pages = [
        "index.md",
        "Example" => [
            "setup.md",
            "running.md",
            "saving.md",
            "plotting.md",
        ],
        "api_reference.md"
    ],
    format = Documenter.HTML(size_threshold=10000000),
    # draft = true,
)