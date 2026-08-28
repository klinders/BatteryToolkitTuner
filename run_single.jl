## 
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using BatteryToolkit
using Setfield
using SciMLStructures
import Logging: global_logger
import TerminalLoggers: TerminalLogger

global_logger(TerminalLogger())

t_end = 150
times = collect(1:1:t_end)
day = 3600 * 24

function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    t = round(int.t,digits=2)
    @warn "Simulation step terminated at t=$t"
    return (;)
end

function charge!(mod, obs, ctx, int)
    # First charge needs to be after V=2.5 so only trigger after one discharge
    if int.t > 3400
        @set! mod.Iin = 4.89*0.3
    end
end

function discharge!(mod, obs, ctx, int)
    @set! mod.Iin = -4.89
end


function SingleCell(;name, params=Chen2020(), Qcell=4.89, T=298.15, kargs...)

    D = Differential(t)

    ModelingToolkit.@variables begin 
        # Pin(t)=0, [input=true]
        Iin(t)=0, [input=true]
        Tin(t)=T, [input=true]
        soc_min(t) = 0.0, [input=true]
        soc_max(t) = 0.3, [input=true]
        V(t)
        I(t)
    end

    @named cell = SPMe(params=params, Q=Qcell) # Battery cell model
    @named source = Current()
    @named ground = Ground()

    eqs = [
        V ~ cell.v
        I ~ cell.i
        D(Tin) ~ 0
        D(Iin) ~ 0
        D(soc_min) ~ 0
        D(soc_max) ~ 0

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        cell.T.u ~ Tin

        source.I.u ~  Iin
    ]

    events = [
        [cell.v ~ 2.5]=>(charge!,(;Iin)),
        [cell.v ~ 4.2]=>(discharge!, (;Iin)),
        [cell.soc ~ soc_min]=>(charge!,(;Iin)),
        [cell.soc ~ soc_max]=>(discharge!,(;Iin))
    ]

    return System(eqs, t; systems=[cell, source, ground],continuous_events=events, name=name, kargs...)
end

@mtkcompile sys = SingleCell(Qcell=4.96, params=OKane2022())

sys = subset_tunables(sys, [
    sys.cell.cracking_n.k_cr, 
    sys.cell.cracking_n.m_cr, 
    sys.cell.cracking_n.b_cr, 
    sys.cell.plating.k_plating, 
    sys.cell.plating.α_plating,
    sys.cell.plating.γ₀,
    sys.cell.lam_p.β_LAM,
    sys.cell.lam_p.m_LAM,
    sys.cell.lam_n.β_LAM,
    sys.cell.lam_n.m_LAM,
])

using OrdinaryDiffEq

calendar_results = [
    [sys.cell.α=>0.2000000000000000, sys.cell.k_sei=>2.141364510787677e-15, sys.cell.D_ec=>7.271349692756622e-21, sys.cell.cathode_diss.i₀=>0.007631492283260097, sys.cell.cathode_diss.Ediss=>4.5],
    [sys.cell.α=>0.3847905431966652, sys.cell.k_sei=>4.552377035512908e-16, sys.cell.D_ec=>2.926256584834017e-22, sys.cell.cathode_diss.i₀=>1.0e-6, sys.cell.cathode_diss.Ediss=>3.5],
    [sys.cell.α=>0.7645442059411685, sys.cell.k_sei=>6.141416974111299e-15, sys.cell.D_ec=>3.0989558838167383e-22, sys.cell.cathode_diss.i₀=>1.0e-6, sys.cell.cathode_diss.Ediss=>3.5],
]

sys = debug_system(sys)

prob = ODEProblem(sys, [sys.Iin=>-5, sys.Tin=>298.15, sys.soc_min=>0.0, sys.soc_max=>0.3], (0, t_end*day))
sols=[]

for (i,T) in enumerate([10, 25, 40])
    nprob = remake(prob, u0=[sys.Tin=>273.15+T],p=calendar_results[2], build_initializeprob = false)

    sol = solve(nprob, DefaultODEAlgorithm(), dtmax=1, saveat=day, progress=true, verbose=true)
    @show sol.retcode
    push!(sols,sol)
end


## Plot los
using CairoMakie
using ColorSchemes

f = Figure()

ax = Axis(f[1,1],
    title="Capacity",
    xlabel="Time [days]",
    ylabel="Capacity [A.h]",
)

colors = get(ColorSchemes.matter, range(0,1, length=3))

for (i,sol) in enumerate(sols)
    # val = minimum.(sol[sys.cell.el.cₑ])
    CairoMakie.lines!(ax, sol[sys.t]/24/3600, sol[sys.cell.ne.c_surf], label="Exp $(i)", color=colors[i])
    # CairoMakie.scatter!(ax, real_data[p].t, real_data[p].q/1000, label="Exp $(i)", color=colors[i])
end

axislegend(ax)
# CairoMakie.xlims!(0, t_end)
@show f