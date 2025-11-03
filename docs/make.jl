using PeriodicARModels
using Documenter

DocMeta.setdocmeta!(PeriodicARModels, :DocTestSetup, :(using PeriodicARModels); recursive=true)

makedocs(;
    modules=[PeriodicARModels],
    authors="ArnaudG0649 <arnaudcmc@hotmail.com> and contributors",
    sitename="PeriodicARModels.jl",
    format=Documenter.HTML(;
        canonical="https://ArnaudG0649.github.io/PeriodicARModels.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArnaudG0649/PeriodicARModels.jl",
    devbranch="master",
)
