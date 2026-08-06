using ThermalResponseDeconvolution
using Documenter

# Make `using ThermalResponseDeconvolution` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    ThermalResponseDeconvolution,
    :DocTestSetup,
    :(using ThermalResponseDeconvolution);
    recursive = true,
)

makedocs(;
    modules = [ThermalResponseDeconvolution],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "ThermalResponseDeconvolution.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/ThermalResponseDeconvolution.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => "theory.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/ThermalResponseDeconvolution.jl",
    devbranch = "main",
)
