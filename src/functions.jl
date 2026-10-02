using CSV
using DataFrames
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using BatteryToolkit
using Setfield
using SciMLStructures
using OrdinaryDiffEq
using SymbolicIndexingInterface
import Logging: global_logger
import TerminalLoggers: TerminalLogger
using DiffEqGPU
using oneAPI

global_logger(TerminalLogger())
# using ProgressLogging

function rest!(mod, obs, ctx, int)
    if mod.Iin > 0
        @set! mod.Iin = 0
    end
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

function SingleCell(;name, params=OKane2022(), Qcell=4.89, T=298.15, calendar_ageing=false, kargs...)

    D = Differential(t)

    ModelingToolkit.@variables begin 
        # Pin(t)=0, [input=true]
        Iin(t)=0, [input=true]
        Tin(t)=T, [input=true]
        soc_max(t) = 0.3, [input=true]
        soc_min(t) = 0.0, [input=true]
        V(t)
        I(t)
    end

    @named cell = SPMe(params=params, Q=Qcell) # Battery cell model
    @named source = Current()
    @named ground = Ground()

    D = Differential(t)

    eqs = [
        V ~ cell.v
        I ~ cell.i
        
        D(Tin) ~ 0
        D(Iin) ~ 0
        D(soc_max) ~ 0
        D(soc_min) ~ 0

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        cell.T.u ~ Tin
        source.I.u ~  Iin
    ]

    
    if calendar_ageing
        events = [
            [cell.v ~ 2.5]=>(charge!, (;Iin)),
            [cell.soc ~ soc_max]=>(rest!, (;Iin))
        ]
    else
        events = [
            [cell.v ~ 2.5]=>(charge!,(;Iin)),
            [cell.v ~ 4.2]=>(discharge!, (;Iin)),
            [cell.soc ~ soc_min]=>(charge!,(;Iin)),
            [cell.soc ~ soc_max]=>(discharge!,(;Iin))
        ]
    end

    return System(eqs, t; name=name, continuous_events=events,systems=[cell,source,ground], kargs...)
end

function test_all_cases(sys, cases, args...; parameters=[], period=(0, 24*3600*365), kargs...)
    N = length(cases)
    u0 = [sys.Iin=>-4.89, sys.Tin=>298.15, sys.soc_min=>0.0, sys.soc_max=>0.3]
    u0p0 = length(parameters)>0 ? [u0;parameters] : u0

    prob = ODEProblem(sys, u0p0, period)
    
    sols = []

    function prob_func(prob, i, repeat)
        # Iin, Tin, and end_soc are constant input variables. In the compiled
        # system they are initial parameters, not entries in the state vector.
        return remake(prob, u0 = cases[i], build_initializeprob = false)
    end

    eprob = EnsembleProblem(prob, prob_func=prob_func)
    sols = solve(eprob, args...; trajectories=N, saveat=24*3600, progress=true, kargs...)
    return sols
end

using MAT
using CSV

function load_datasets()
    real_data = Vector{DataFrame}()
    
    # first calendar
    data = matread(joinpath(@__DIR__,"../data/Kuzhiyil/RPT_analysis_data.mat"))
    
    for T in [0, 25, 45]
        for soc in [30, 50, 80]
            t_vec = [a[1] for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Days"])]
            q_vec = [a[1]*1000 for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Capacity_Mean"])]

            push!(real_data, DataFrame((t=t_vec, q=q_vec)))
        end
    end

    # Cyclic
    for T in [10, 25, 40]
        for exp in [1, 2, 3]
            data = CSV.read(joinpath(@__DIR__,"../data/Kirkaldy/Expt $(exp) - $(T)degC - Processed Data.csv"), DataFrame)
            rename!(data, [
                "Days of degradation"=>:t,
                "Charge Throughput [A h]"=>:q_ah,
                "NE Capacity [mA h]"=>:q_n, 
                "PE Capacity [mA h]"=>:q_p, 
                "Cell Capacity [mA h]"=>:q, 
                "LAM PE"=>:lam_p,
                "LAM NE_tot"=>:lam_n,
                "LLI"=>:lli,
                "SoH"=>:soc,
                "0.1s Resistance [Ohms]"=>:r
            ])
            push!(real_data, data)
        end
    end

    return real_data

end