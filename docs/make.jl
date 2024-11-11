using ddHodge
using Documenter

DocMeta.setdocmeta!(ddHodge, :DocTestSetup, :(using ddHodge); recursive=true)

makedocs(;
    modules=[ddHodge],
    authors="Kazumitsu Maehara <kazumits@gmail.com> and contributors",
    sitename="ddHodge.jl",
    format=Documenter.HTML(;
        canonical="https://kazumits.github.io/ddHodge.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kazumits/ddHodge.jl",
    devbranch="main",
)
