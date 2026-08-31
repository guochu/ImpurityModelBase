using ImpurityModelBase
using Documenter

DocMeta.setdocmeta!(ImpurityModelBase, :DocTestSetup, :(using ImpurityModelBase); recursive=true)

makedocs(;
    modules=[ImpurityModelBase],
    authors="Guo Chu <guochu604b@gmail.com> and contributors",
    sitename="ImpurityModelBase",
    format=Documenter.HTML(;
        mathengine=Documenter.MathJax(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/guochu/ImpurityModelBase",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "简介" => "index.md",
        "谱密度函数" => "spectrumfuncs.md",
        "粒子浴" => "baths.md",
        "精确对角化" => "exactdiagonalizations.md",
        "解析解" => "analyticsolutions.md",
        "工具函数" => "utilities.md",
    ],
)

deploydocs(;
    repo="github.com/guochu/ImpurityModelBase.git",
    devbranch="master",
)
