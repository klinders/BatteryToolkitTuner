using CSV
using DataFrames
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using BatteryToolkit
using Setfield
using SciMLStructures
using OrdinaryDiffEq, Sundials
using SymbolicIndexingInterface
# import Logging: global_logger
# import TerminalLoggers: TerminalLogger
# global_logger(TerminalLogger())
using ProgressLogging

first = true

function charge!(mod, obs, ctx, int)
    global first
    if !first
        return (; Iin = 4.89*0.3)
    else
        first = false
        return (; Iin = mod.Iin)
    end
end

function discharge!(mod, obs, ctx, int)
    return (; Iin = -4.89)
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
        ModelingToolkit.SymbolicContinuousCallback([
            cell.v ~ 2.5
        ], ModelingToolkit.ImperativeAffect(modified=(;Iin)) do x,o,c,i
            charge!(x, o, c, i)
        end),

        ModelingToolkit.SymbolicContinuousCallback([
            cell.v ~ 4.2
        ], ModelingToolkit.ImperativeAffect(modified=(;Iin)) do x,o,c,i
            discharge!(x, o, c, i)
        end),

        ModelingToolkit.SymbolicContinuousCallback([
            cell.soc ~ soc_min
        ], ModelingToolkit.ImperativeAffect(modified=(;Iin)) do x,o,c,i
            charge!(x, o, c, i)
        end),

        ModelingToolkit.SymbolicContinuousCallback([
            cell.soc ~ soc_max
        ], ModelingToolkit.ImperativeAffect(modified=(;Iin)) do x,o,c,i
            discharge!(x, o, c, i)
        end)
    ]

    return System(eqs, t; systems=[cell, source, ground],continuous_events=events, name=name, kargs...)
end

function test_all_cases(sys, parameters, cases)
    N = length(cases)
    u0 = [sys.Iin=>-4.89, sys.Tin=>298.15, sys.soc_min=>0.0, sys.soc_max=>0.3]
    u0p0 = length(parameters)>0 ? [u0;parameters] : u0

    prob = ODEProblem(sys, u0p0, (0, 24*3600*365))
    
    sols = []

    for case in cases
        nprob = remake(prob, u0 = case, build_initializeprob = false)
        sol = solve(nprob, QNDF(), maxiters=1e8, trajectories=N, saveat=24*3600, progress=true)
        push!(sols,sol)
    end

    return sols
end

using MAT
using CSV

function load_datasets()
    real_data = Vector{DataFrame}()
    
    # first calendar
    data = matread(joinpath(@__DIR__,"../parameter tuning/Kuzhiyil/RPT_analysis_data.mat"))
    
    for T in [0, 25, 45]
        for soc in [30, 50, 80]
            t_vec = [a[1] for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Days"])]
            q_vec = [a[1]*1000 for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Capacity_Mean"])]

            push!(real_data, DataFrame((t=t_vec, q=q_vec)))
        end
    end

    # Cyclic
    for exp in [1, 2, 3]
        for T in [10, 25, 40]
            data = CSV.read(joinpath(@__DIR__,"Kirkaldy/Expt $(exp) - $(T)degC - Processed Data.csv"), DataFrame)
            rename!(data, ["Days of degradation"=>:t,  "NE Capacity [mA h]"=>:q_n, "PE Capacity [mA h]"=>:q_p, "Cell Capacity [mA h]"=>:q])

            push!(real_data, data)
        end
    end

    return real_data

end