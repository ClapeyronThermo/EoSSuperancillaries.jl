using Documenter,EoSSuperancillaries

makedocs(sitename = "EoSSuperancillaries.jl",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        canonical = "https://ClapeyronThermo.github.io/EoSSuperancillaries.jl/",
        ),
    warnonly = Documenter.except(),
    authors = "Pierre J. Walker and Andrés Riedemann.",
    pages = ["Home" => "index.md",
            "API" =>["van der Wals" => "vdw.md",
            "Redlich-Kwong" => "rk.md",        
            "Peng-Robinson" => "pr.md",
            "PC-SAFT" => "pcsaft.md"
            "REFPROP" => [
                "REFPROP superancillaries" => "refprop.md",
                "list of fluids" => "refprop_fluids.md"]
            ]]
    )

deploydocs(repo="github.com/ClapeyronThermo/EoSSuperancillaries.jl.git")
