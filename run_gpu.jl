using CSV
using DataFrames
using Logging, LoggingExtras

global_logger(MinLevelLogger(FileLogger("output/log.txt"), Info))

temperatures = [10, 25, 40]
socs = [[0,30], [70,85], [85,100]]

labels = ["T$(temp)SOC$(soc[1])-$(soc[2])" for temp in temperatures for soc in socs]

# data["CalResults"][
#   "rptDates" => Dates of when the RPT is done
#   "Resistance"=> [
#       "DisRes"=>
#       "ChaRes"=>"Pulse1CR10"=>[
#           "ResChaErr"=>
#           "ResChaMean"=>
#       ]
#    ]
#    "Capacity"=>[
#         "meanCapStdErr"
#         "meanTotalCap"
#         "meanCap"
#         "indTotalCap"
#         "meanTotalCapStdErr"
#         "indCap"
#    ]
# ]

using Interpolations
# save date every week
t_end = 150
times = collect(1:1:t_end)
day = 3600 * 24


## 
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using BatteryToolkit
using Setfield
using SciMLStructures
using DiffEqGPU
using CUDA

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

@info "Compiling"
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


prob = ODEProblem(sys, [sys.Iin=>-5, sys.Tin=>298.15, sys.soc_min=>0.0, sys.soc_max=>0.3], (0, day*times[end]))

##
using OrdinaryDiffEq
using SymbolicIndexingInterface
using CairoMakie
using ColorSchemes
using Optimization
using ForwardDiff
using DiffEqParamEstim
using OptimizationOptimJL
using OptimizationMetaheuristics

# ["α", "k_sei", "D_ec", "i₀", "U_diss"]
calendar_results = [
    [sys.cell.α=>0.2000000000000000, sys.cell.k_sei=>2.141364510787677e-15, sys.cell.D_ec=>7.271349692756622e-21, sys.cell.cathode_diss.i₀=>0.007631492283260097, sys.cell.cathode_diss.Ediss=>4.5],
    [sys.cell.α=>0.3847905431966652, sys.cell.k_sei=>4.552377035512908e-16, sys.cell.D_ec=>2.926256584834017e-22, sys.cell.cathode_diss.i₀=>1.0e-6, sys.cell.cathode_diss.Ediss=>3.5],
    [sys.cell.α=>0.7645442059411685, sys.cell.k_sei=>6.141416974111299e-15, sys.cell.D_ec=>3.0989558838167383e-22, sys.cell.cathode_diss.i₀=>1.0e-6, sys.cell.cathode_diss.Ediss=>3.5],
]
results = []

for (i_T,T) in enumerate(temperatures)
    @info "Solving T=$(T)"

    real_data = Vector{DataFrame}()
    
    for exp in eachindex(socs)
        data = CSV.read(joinpath(@__DIR__,"../data/Kirkaldy/Expt $(exp) - $(T)degC - Processed Data.csv"), DataFrame)
        rename!(data, ["Days of degradation"=>:t,  "NE Capacity [mA h]"=>:q_n, "PE Capacity [mA h]"=>:q_p, "Cell Capacity [mA h]"=>:q])

        push!(real_data, data)
    end

    function plot(p::EnsembleSolution, name; title="")
        CairoMakie.set_theme!(theme_latexfonts(), fontsize=24)

        f = Figure()

        ax = Axis(f[1,1],
            title=title,
            xlabel="Time [days]",
            ylabel="Capacity [A.h]",
        )
        
        colors = get(ColorSchemes.matter, range(0,1, length=length(p)))
        
        for i in eachindex(p)
            s = p[i]
            CairoMakie.lines!(ax, s[sys.t]/day, s[sys.cell.C_cell], label="SoC=$(socs[i][1])-$(socs[i][2])", color=colors[i])
            CairoMakie.scatter!(ax, real_data[i].t, real_data[i].q/1000, color=colors[i])
            # fake the data to 0 error
            # data[i] = s[sys.cell.C_cell]
        end

        axislegend(ax)
        CairoMakie.xlims!(0, t_end)

        save(joinpath(@__DIR__,"../plots/cyclic",name), f)
    end

    u0_arr = [
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.soc_min=>0.0, sys.soc_max=>0.3],
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.soc_min=>0.7, sys.soc_max=>0.85],
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.soc_min=>0.85, sys.soc_max=>1.0],
    ]

    N = length(u0_arr)

    function prob_func(prob, i, repeat)
        # Iin, Tin, and end_soc are constant input variables. In the compiled
        # system they are initial parameters, not entries in the state vector.
        return remake(prob, u0 = u0_arr[i], p=calendar_results[i_T], build_initializeprob = false)
    end

    function output_func(sol, i)
        if !SciMLBase.successful_retcode(sol)
            return (sol, true)
        end
        return (sol, false)
        
    end

    function eprob_func(eprob, p)
        ps = parameter_values(eprob.prob)

        # Replace the tunable portion with new values
        new_ps = SciMLStructures.replace(SciMLStructures.Tunable(), ps, p)

        return EnsembleProblem(remake(eprob.prob, p=new_ps), prob_func=prob_func)
    end


    eprob = EnsembleProblem(prob, prob_func=prob_func)
    
    @info "Initial solve started"
    sol_t = solve(eprob, trajectories=N, saveat=day)
    plot(sol_t, "not_optimized_$(T).png", title="No optimization at $(T)°C")

    get_q = getsym(sys, sys.cell.C_cell)
    get_qn = getsym(sys, sys.cell.C_neg)
    get_qp = getsym(sys, sys.cell.C_pos)

    function cost(solution)
        err = 0.0
        
        for i in eachindex(solution)
            if solution[i].retcode != SciMLBase.ReturnCode.Success
                err += Inf # Penalize failed simulations
                continue
            end
            q = get_q(solution[i])
            q_n = get_qn(solution[i])
            q_p = get_qp(solution[i])
            
            for d in eachrow(real_data[i])
                # 6 comes from the saves during the three events
                if d.t > t_end
                    continue
                end
                index = round(Int, d.t + 6)
                err += abs2.(q[index] .- d.q./1000) + abs2.(q_n[index] .- d.q_n./1000) + abs2.(q_p[index] .- d.q_p./1000)
            end
        end
        print("\e[A\e[2K")
        @info "Optimizing... Loss: $err"
        return err
    end

    @info "Initial cost: $(cost(sol_t))"

    obj = build_loss_objective(eprob, QNDF(), cost, Optimization.AutoFiniteDiff(),
                               prob_generator = eprob_func,
                               trajectories = N,
                               maxiters=1e8,
                               saveat = day);
    lb = [1e-12, 0.3, 1e-9, 1.5, 1.0, 1e-23, 1e-10, 1.0, 1e-10, 1.0]
    ub = [1e-6,  0.7, 1e-3, 3.0, 1.5, 1e-17, 1e-4,  3.0, 1e-4,  3.0]
    
    p0 = [    
        1e-10,     # k_plating
        0.65,     # α_plating
        1e-7,    # γ₀
        2.2,     # m_cr
        1.12,    # b_cr
        3.9e-20, # k_cr
        1.47e-07, # β_LAM neg
        1.02,        # m_LAM neg
        3.43e-07, # β_LAM pos
        1.02,         # m_LAM pos
    ]

    optprob = OptimizationProblem(obj, p0, lb = lb, ub = ub);
    
    function cb(state, loss)
        print("\e[A\e[2K")
        @info "State:$state, Loss: $loss"
    end
    
    result = solve(optprob, PSO(N=15, C1=2.0, C2=2.0, ω=0.8),  maxiters=100, use_initial=true)

    push!(results, result)

    @info "Done optimizing cyclic ageing T=$(T)"
    @info result.u
    
    ## Optimized
    eprob_opt = eprob_func(eprob, result.u)
    
    sol = solve(eprob_opt, trajectories=N, saveat=day)
    
    plot(sol, "with_optimized_$(T).png", title="With optimization at $(T)°C")

end

# T=0 
# result = retcode: Default
# u: [0.2, 2.141364510787677e-15, 7.271349692756622e-21, 0.007631492283260097, 4.5]
# Final objective value:     0.02468261205925673

# T=25
# result = retcode: Default
# u: [0.38479054319666517, 4.552377035512908e-16, 2.926256584834017e-22, 1.0e-6, 3.5]
# Final objective value:     0.012139333538911034

# T=45
# result = retcode: Default
# u: [0.7645442059411685, 6.141416974111299e-15, 3.0989558838167383e-22, 1.0e-6, 3.5]
# Final objective value:     0.02486319054475867

## Plot the results
# results = [
#     [0.2, 2.141, 72.71, 0.007631, 4.5],
#     [0.3847, 0.4552, 2.926, 1.0e-6, 3.5],
#     [0.7645, 6.141, 3.098, 1.0e-6, 3.5]
# ]
@info "Done" 
@info results

# parameters = ["α", "k_sei", "D_ec", "i₀", "U_diss"]
# colors = get(ColorSchemes.matter, range(0,1, length=2))
# GLMakie.set_theme!(theme_latexfonts(), fontsize=24)

# f = Figure(size=(1600, 400))

# for i in eachindex(parameters)
#     ax = Axis(f[1,i], title="$(parameters[i])", xlabel="Temperature [°C]", ylabel="")
#     # ax.ytickformat = "{:.2f}"
    
#     # GLMakie.tight_ticklabel_spacing!(ax)
#     y = [results[j][i] for j in 1:length(temperatures)]
#     GLMakie.scatter!(ax, temperatures, y)
#     # Sample data
#     X = hcat(temperatures, ones(length(y)))
#     coefs = X \ y
#     a, b = coefs[1], coefs[2]
    
#     GLMakie.lines!(ax, temperatures, a .* temperatures.+b, linestyle=:dash)
# end

# Label(f[1, 2, Top()], L"\times 10^-15", halign = :left, valign = :bottom, padding = (0, 0, 5, 0))
# Label(f[1, 3, Top()], L"\times 10^-22", halign = :left, valign = :bottom, padding = (0, 0, 5, 0))

# save(joinpath(@__DIR__,"../plots/cyclic","optimized_parameters.png"), f)

# @show f