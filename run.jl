
using GLMakie
using ColorSchemes
using CSV
using MAT
using JLD2

include("functions.jl")

@mtkcompile sys = SingleCell(Qcell=4.96, params=OKane2022())

parameters = []
cases = [
    # Calendar conditions
    [sys.Iin=>-5, sys.Tin=>273.15+0, sys.soc_min=>0.3, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+0, sys.soc_min=>0.5, sys.soc_max=>0.5],
    [sys.Iin=>-5, sys.Tin=>273.15+0, sys.soc_min=>0.8, sys.soc_max=>0.8],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.3, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.5, sys.soc_max=>0.5],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.8, sys.soc_max=>0.8],
    [sys.Iin=>-5, sys.Tin=>273.15+45, sys.soc_min=>0.3, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+45, sys.soc_min=>0.5, sys.soc_max=>0.5],
    [sys.Iin=>-5, sys.Tin=>273.15+45, sys.soc_min=>0.8, sys.soc_max=>0.8],
    # Cycle 0-30% expt1
    [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.0, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.0, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.0, sys.soc_max=>0.3],
    # Cycle 70-85%
    [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.7, sys.soc_max=>0.85],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.7, sys.soc_max=>0.85],
    [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.7, sys.soc_max=>0.85],
    #cycle 85-100%
    [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.85, sys.soc_max=>1],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.85, sys.soc_max=>1],
    [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.85, sys.soc_max=>1],
]

sols = test_all_cases(sys, parameters, cases)
@save "output.jld2" sols

plots = [
    (
        title="Calendar 0degC",
        filename="cal_0.png",
        datasets=[1,2,3],
        labels=["30% SoC", "50% SoC", "80% SoC"]
    ),
    (
        title="Calendar 25degC",
        filename="cal_25.png",
        datasets=[4,5,6],
        labels=["30% SoC", "50% SoC", "80% SoC"]
    ),    
    (
        title="Calendar 45degC",
        filename="cal_45.png",
        datasets=[7,8,9],
        labels=["30% SoC", "50% SoC", "80% SoC"]
    ),
    (
        title="Cyclic 0-30%",
        filename="cyc_0-30.png",
        datasets=[10,11,12],
        labels=["10degC", "25degC", "40degC"]
    ),
    (
        title="Cyclic 70-85%",
        filename="cyc_70-85.png",
        datasets=[13,14,15],
        labels=["10degC", "25degC", "40degC"]
    ),
        (
        title="Cyclic 85-100%",
        filename="cyc_85-100.png",
        datasets=[16,17,18],
        labels=["10degC", "25degC", "40degC"]
    ),
]

real_data = load_datasets()

## Plot
GLMakie.set_theme!(theme_latexfonts(), fontsize=24)

for plt in plots
    
    f = Figure()

    ax = Axis(f[1,1],
        title=plt.title,
        xlabel="Time [days]",
        ylabel="Capacity [A.h]",
    )

    colors = get(ColorSchemes.matter, range(0,1, length=length(plt.datasets)))
    
    for (i,p) in enumerate(plt.datasets)
        GLMakie.lines!(ax, sols[p][sys.t]/24/3600, sols[p][sys.cell.C_cell], label=plt.labels[i], color=colors[i])
        GLMakie.scatter!(ax, real_data[p].t, real_data[p].q/1000, color=colors[i])
    end

    axislegend(ax)
    GLMakie.xlims!(0, 365)

    save(joinpath(@__DIR__,"plots",plt.filename), f)
end