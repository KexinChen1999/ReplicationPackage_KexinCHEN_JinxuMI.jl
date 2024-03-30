using ReplicationPackage_KexinCHEN_JinxuMI
using Documenter

DocMeta.setdocmeta!(ReplicationPackage_KexinCHEN_JinxuMI, :DocTestSetup, :(using ReplicationPackage_KexinCHEN_JinxuMI); recursive=true)

makedocs(;
    modules=[ReplicationPackage_KexinCHEN_JinxuMI],
    authors="KexinChen1999 <chenkexin.vincent@gmail.com> and contributors",
    sitename="ReplicationPackage_KexinCHEN_JinxuMI.jl",
    format=Documenter.HTML(;
        canonical="https://KexinChen1999.github.io/ReplicationPackage_KexinCHEN_JinxuMI.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KexinChen1999/ReplicationPackage_KexinCHEN_JinxuMI.jl",
    devbranch="main",
)
