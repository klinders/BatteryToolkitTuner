
using GLMakie
using ColorSchemes
using CSV
using MAT
using JLD2
using OrdinaryDiffEq
using Revise

include("src/functions.jl")


## Calendar ageing
@mtkcompile sys = SingleCell(calendar_ageing=false)

parameters = [
    sys.cell.α=>0.3847905431966652, 
    sys.cell.k_sei=>4.552377035512908e-16, 
    sys.cell.D_ec=>2.926256584834017e-22, 
    sys.cell.cathode_diss.i₀=>1.0e-6, 
    sys.cell.cathode_diss.Ediss=>3.5,
    sys.cell.plating.k_plating=>1e-10,
    sys.cell.plating.γ₀=>1e-7,
    sys.cell.plating.α_plating=>0.65,  
    sys.cell.cracking_n.k_cr=>3.9e-20,
    sys.cell.cracking_n.m_cr=>2.2,
    sys.cell.cracking_n.b_cr=>1.12,
    sys.cell.lam_n.β_LAM=>1.47e-7,
    sys.cell.lam_p.β_LAM=>3.43e-7,
    sys.cell.lam_n.m_LAM=>1.02,
    sys.cell.lam_p.m_LAM=>1.02
]


cases = [
    # Calendar conditions
    # [sys.Iin=>-4.89, sys.Tin=>273.15+0,sys.soc_min=>0.0, sys.soc_max=>0.3],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+0,sys.soc_min=>0.0, sys.soc_max=>0.5],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+0,sys.soc_min=>0.0, sys.soc_max=>0.8],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+25, sys.soc_min=>0.3, sys.soc_max=>0.3],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+25, sys.soc_min=>0.5, sys.soc_max=>0.5],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+25, sys.soc_min=>0.8, sys.soc_max=>0.8],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+45, sys.soc_min=>0.3, sys.soc_max=>0.3],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+45, sys.soc_min=>0.5, sys.soc_max=>0.5],
    # [sys.Iin=>-4.89, sys.Tin=>273.15+45, sys.soc_min=>0.8, sys.soc_max=>0.8],
    # Cycle 0-30% expt1
    # [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.0, sys.soc_max=>0.3],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.0, sys.soc_max=>0.3],
    # [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.0, sys.soc_max=>0.3],
    # Cycle 70-85%
    # [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.7, sys.soc_max=>0.85],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.7, sys.soc_max=>0.85],
    # [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.7, sys.soc_max=>0.85],
    #cycle 85-100%
    # [sys.Iin=>-5, sys.Tin=>273.15+10, sys.soc_min=>0.85, sys.soc_max=>1],
    [sys.Iin=>-5, sys.Tin=>273.15+25, sys.soc_min=>0.85, sys.soc_max=>1],
    # [sys.Iin=>-5, sys.Tin=>273.15+40, sys.soc_min=>0.85, sys.soc_max=>1],
]

sols = test_all_cases(sys, cases, QNDF(); parameters=parameters, period=(0, 24*3600*150), maxiters=1e8)

# mkdir("output/unoptimized")
# for (i,sol) in enumerate(sols)
#     @save "output/unoptimized/Case_$(i).jld2" sol
# end

# plots = [
#     (
#         title="Calendar 0degC",
#         filename="cal_0.png",
#         datasets=[1,2,3],
#         labels=["30% SoC", "50% SoC", "80% SoC"]
#     ),
#     (
#         title="Calendar 25degC",
#         filename="cal_25.png",
#         datasets=[4,5,6],
#         labels=["30% SoC", "50% SoC", "80% SoC"]
#     ),    
#     (
#         title="Calendar 45degC",
#         filename="cal_45.png",
#         datasets=[7,8,9],
#         labels=["30% SoC", "50% SoC", "80% SoC"]
#     ),
#     (
#         title="Cyclic 0-30%",
#         filename="cyc_0-30.png",
#         datasets=[10,11,12],
#         labels=["10degC", "25degC", "40degC"]
#     ),
#     (
#         title="Cyclic 70-85%",
#         filename="cyc_70-85.png",
#         datasets=[13,14,15],
#         labels=["10degC", "25degC", "40degC"]
#     ),
#         (
#         title="Cyclic 85-100%",
#         filename="cyc_85-100.png",
#         datasets=[16,17,18],
#         labels=["10degC", "25degC", "40degC"]
#     ),
# ]

real_data = load_datasets()[[10,13,16]]

# ## Plot
# GLMakie.set_theme!(theme_latexfonts(), fontsize=24)

# for plt in plots
    
#     f = Figure()

#     ax = Axis(f[1,1],
#         title=plt.title,
#         xlabel="Time [days]",
#         ylabel="Capacity [A.h]",
#     )

#     colors = get(ColorSchemes.matter, range(0,1, length=length(plt.datasets)))
    
#     for (i,p) in enumerate(plt.datasets)
#         GLMakie.lines!(ax, sols[p][sys.t]/24/3600, sols[p][sys.cell.C_cell], label=plt.labels[i], color=colors[i])
#         GLMakie.scatter!(ax, real_data[p].t, real_data[p].q/1000, color=colors[i])
#     end

#     axislegend(ax)
#     GLMakie.xlims!(0, 365)

#     save(joinpath(@__DIR__,"plots",plt.filename), f)
# end

## Plot los
colors = get(ColorSchemes.managua, range(0,1, length=3))

fig = Figure(size=(1200,800))
q_init = 4.89

for (i,sol) in enumerate(sols)
    @info sol.retcode
    ax = Axis(fig[1, i], title="Experiment $(i)", xlabel="Time [days]", ylabel="Capacity loss [%]")
    band!(ax, sol[sys.t]/3600/24, zeros(length(sol[sys.t])), sol[sys.cell.sei.Q_loss]./q_init*100, label="SEI")
    band!(ax, sol[sys.t]/3600/24, (sol[sys.cell.sei.Q_loss])./q_init*100, (sol[sys.cell.sei.Q_loss].+sol[sys.cell.plating.Q_dead])./q_init*100, label="Plating")
    band!(ax, sol[sys.t]/3600/24, (sol[sys.cell.sei.Q_loss].+sol[sys.cell.plating.Q_dead])./q_init*100, (sol[sys.cell.sei.Q_loss]+sol[sys.cell.plating.Q_dead].+sol[sys.cell.cracking_n.Q_sei])./q_init*100, label="Cracking")
    band!(ax, sol[sys.t]/3600/24, (sol[sys.cell.sei.Q_loss]+sol[sys.cell.plating.Q_dead].+sol[sys.cell.cracking_n.Q_sei])./q_init*100, (sol[sys.cell.sei.Q_loss]+sol[sys.cell.plating.Q_dead].+sol[sys.cell.cracking_n.Q_sei].+sol[sys.cell.lam_p.Q_loss])./q_init*100, label="LAM+")
    band!(ax, sol[sys.t]/3600/24, (sol[sys.cell.sei.Q_loss]+sol[sys.cell.plating.Q_dead].+sol[sys.cell.cracking_n.Q_sei].+sol[sys.cell.lam_p.Q_loss])./q_init*100, (sol[sys.cell.sei.Q_loss]+sol[sys.cell.plating.Q_dead].+sol[sys.cell.cracking_n.Q_sei].+sol[sys.cell.lam_p.Q_loss].+sol[sys.cell.lam_n.Q_loss])./q_init*100, label="LAM-")
    # lines!(ax, df.t/3600/24, (df.q_loss)/q_init*100, label="Total")
    scatter!(ax, real_data[i].t, (1 .- real_data[i].q/real_data[i].q[1])*100)
    xlims!(ax, (0, 150))
    
    ax2 = Axis(fig[2,i], xlabel="Time [days]", ylabel="Capacity [A.h]")
    lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.LAMₚ], label="LAMₚ", color=colors[1])
    lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.LAMₙ], label="LAMₙ", color=colors[2])
    lines!(ax2, sol[sys.t]/3600/24, sol[sys.cell.Q_loss]/q_init*100, label="Total", color=colors[3])
    scatter!(ax2, real_data[i].t, (1 .- real_data[i].q_p./real_data[i].q_p[1])*100, color=colors[1])
    scatter!(ax2, real_data[i].t, (1 .- real_data[i].q_n./real_data[i].q_n[1])*100, color=colors[2])
    scatter!(ax2, real_data[i].t, (1 .- real_data[i].q./real_data[i].q[1])*100, color=colors[3])
    xlims!(ax2, (0, 150))
    
    if i == length(sols)
        axislegend(ax, position=:lt)
        axislegend(ax2, position=:lt)
    end
end
# place a legend on the first subplot (adjust as desired)

fig
@show(fig)