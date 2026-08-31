using ImpurityModelBase
using Documenter

DocMeta.setdocmeta!(ImpurityModelBase, :DocTestSetup, :(using ImpurityModelBase); recursive=true)

makedocs(;
    modules=[ImpurityModelBase],
    authors="Guo Chu <guochu604b@gmail.com> and contributors",
    sitename="ImpurityModelBase",
    format=Documenter.HTML(;
        mathengine=Documenter.MathJax3(),
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

# Point MathJax at the locally vendored copy (docs/src/assets/mathjax-tex-svg-full.js)
# instead of the CDN, so that formulas render even without network access.
# documenterBaseURL is set per-page by Documenter ("." at the root, ".." in subdirs).
let documenter_js = joinpath(@__DIR__, "build", "assets", "documenter.js")
    s = read(documenter_js, String)
    s2 = replace(
        s,
        "script.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg-full.js';" =>
            "script.src = documenterBaseURL + '/assets/mathjax-tex-svg-full.js';",
    )
    (s2 == s) || write(documenter_js, s2)
end

deploydocs(;
    repo="github.com/guochu/ImpurityModelBase.git",
    devbranch="master",
)
