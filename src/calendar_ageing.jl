using MAT
using DataFrames

temperatures = [0, 25, 45]
socs = [30, 50, 80]

labels = ["T$(temp)SOC$(soc)" for temp in temperatures for soc in socs]

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
t_end = 365
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

function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    t = round(int.t,digits=2)
    @warn "Simulation step terminated at t=$t"
    return (;)
end

function charge!(mod, obs, ctx, int)
    @set! mod.Iin = 5*0.3
end

function rest!(mod, obs, ctx, int)
    if mod.Iin > 0
        @set! mod.Iin = 0
    end
end


function SingleCell(;name, params=Chen2020(), Qcell=4.89, T=298.15, kargs...)

    D = Differential(t)

    ModelingToolkit.@variables begin 
        # Pin(t)=0, [input=true]
        Iin(t)=0, [input=true]
        Tin(t)=T, [input=true]
        end_soc(t) = 0.5, [input=true]
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
        D(end_soc) ~ 0

        connect(source.n, cell.n)
        connect(source.p, cell.p)
        connect(ground.g, source.n)
        cell.T.u ~ Tin

        source.I.u ~  Iin
    ]

    events = [
        [cell.v ~ 2.5]=>(charge!, (;Iin)),
        [cell.soc ~ end_soc]=>(rest!, (;Iin))
        # ModelingToolkit.SymbolicContinuousCallback([
        #     cell.soc ~ 0.5
        # ], ModelingToolkit.ImperativeAffect(modified=(;Iin)) do x,o,c,i 
        #     if x.Iin > 0
        #         @set! x.Iin = 0
        #     end
        # end)
    ]

    return System(eqs, t; systems=[cell, source, ground],continuous_events=events, name=name, kargs...)
end

@mtkcompile sys = SingleCell(Qcell=4.96, params=OKane2022())

sys = subset_tunables(sys, [sys.cell.D_ec, sys.cell.k_sei, sys.cell.α, sys.cell.cathode_diss.i₀])
prob = ODEProblem(sys, [sys.Iin=>-5, sys.Tin=>298.15, sys.end_soc=>0.3], (0, day*times[end]))

##
using OrdinaryDiffEq
using SymbolicIndexingInterface
using GLMakie
using ColorSchemes
using Optimization
using ForwardDiff
using DiffEqParamEstim
using OptimizationOptimJL
using OptimizationMetaheuristics

results = []
data = matread(joinpath(@__DIR__,"../data/Kuzhiyil/RPT_analysis_data.mat"))

for T in temperatures
    real_data = Vector{DataFrame}()
    
    for soc in socs
        t_vec = [a[1] for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Days"])]
        q_vec = [a[1] for a in eachrow(data["Temperature_$(T)"]["SOC_$(soc)"]["Capacity_Mean"])]

        push!(real_data, DataFrame((t=t_vec, q=q_vec)))
    end

    function plot(p::EnsembleSolution, name; title="")
        GLMakie.set_theme!(theme_latexfonts(), fontsize=24)

        f = Figure()

        ax = Axis(f[1,1],
            title=title,
            xlabel="Time [days]",
            ylabel="Capacity [A.h]",
        )

        
        colors = get(ColorSchemes.matter, range(0,1, length=length(p)))
        
        for i in eachindex(p)
            s = p[i]
            GLMakie.lines!(ax, s[sys.t]/day, s[sys.cell.C_cell], label="SoC=$(socs[i])", color=colors[i])
            GLMakie.scatter!(ax, real_data[i].t, real_data[i].q, color=colors[i])
            # fake the data to 0 error
            # data[i] = s[sys.cell.C_cell]
        end

        axislegend(ax)
        GLMakie.xlims!(0, t_end)

        save(joinpath(@__DIR__,"../plots/calendar",name), f)
    end

    u0_arr = [
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.end_soc=>0.3],
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.end_soc=>0.5],
        [sys.Iin=>-5, sys.Tin=>273.15+T, sys.end_soc=>0.8],
    ]

    N = length(u0_arr)

    function prob_func(prob, i, repeat)
        # Iin, Tin, and end_soc are constant input variables. In the compiled
        # system they are initial parameters, not entries in the state vector.
        return remake(prob, u0 = u0_arr[i], build_initializeprob = false)
    end

    function eprob_func(eprob, p)
        ps = parameter_values(eprob.prob)

        # Replace the tunable portion with new values
        new_ps = SciMLStructures.replace(SciMLStructures.Tunable(), ps, p)

        return EnsembleProblem(remake(eprob.prob, p=new_ps), prob_func=prob_func)
    end

    eprob = EnsembleProblem(prob, prob_func=prob_func)

    sol_t = solve(eprob, trajectories=N, saveat=day)
    plot(sol_t, "not_optimized_$(T).png", title="No optimization at $(T)°C")

    get_q = getsym(sys, sys.cell.C_cell)

    function cost(solution)
        err = 0.0
        
        for i in eachindex(solution)
            if solution[i].retcode != SciMLBase.ReturnCode.Success
                err += Inf # Penalize failed simulations
                continue
            end
            q = get_q(solution[i])
            
            for d in eachrow(real_data[i])
                # 6 comes from the saves during the three events
                if d.t > t_end
                    continue
                end
                index = round(Int, d.t + 6)
                err += abs2.(q[index] .- d.q)
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
                               maxiters = 1e8,
                               saveat = day);
    lb = [0.2, 1e-18, 1e-24, 1e-6, 3.5]
    ub = [0.8, 1e-14, 1e-20, 1e-2, 4.5]
    
    p0 = [0.5, 1e-16, 1e-21, 1e-6, 4]
    
    optprob = OptimizationProblem(obj, p0, lb = lb, ub = ub);
    
    function cb(state, loss)
        print("\e[A\e[2K")
        @info "State:$state, Loss: $loss"
    end
    
    result = solve(optprob, PSO(N=15, C1=2.0, C2=2.0, ω=0.8),  maxiters=50, use_initial=true)

    @show result
    push!(results, result)

    open("../output/results.txt", "a") do io
        write(io, "Calendar ageing T=", T, "degC\n")
        write(io, result,"\n")
    end
    
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
@show results

parameters = ["α", "k_sei", "D_ec", "i₀", "U_diss"]
colors = get(ColorSchemes.matter, range(0,1, length=2))
GLMakie.set_theme!(theme_latexfonts(), fontsize=24)

f = Figure(size=(1600, 400))

for i in eachindex(parameters)
    ax = Axis(f[1,i], title="$(parameters[i])", xlabel="Temperature [°C]", ylabel="")
    # ax.ytickformat = "{:.2f}"
    
    # GLMakie.tight_ticklabel_spacing!(ax)
    y = [results[j][i] for j in 1:length(temperatures)]
    GLMakie.scatter!(ax, temperatures, y)
    # Sample data
    X = hcat(temperatures, ones(length(y)))
    coefs = X \ y
    a, b = coefs[1], coefs[2]
    
    GLMakie.lines!(ax, temperatures, a .* temperatures.+b, linestyle=:dash)
end

Label(f[1, 2, Top()], L"\times 10^-15", halign = :left, valign = :bottom, padding = (0, 0, 5, 0))
Label(f[1, 3, Top()], L"\times 10^-22", halign = :left, valign = :bottom, padding = (0, 0, 5, 0))

save(joinpath(@__DIR__,"../plots/calendar","optimized_parameters.png"), f)

# @show f