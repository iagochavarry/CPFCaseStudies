using Dates
using Statistics
using Plots
using DataFrames
using PowerModels
using CSV
using Ipopt
using JuMP
using Distributions

include("C:\\Users\\iagoc\\Dropbox\\Lamps PUC\\Energisa Must II\\Packages\\BrazilianPowerModels.jl\\src/ControlPowerFlow.jl")

h = 24*7

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-8, "print_level" => 0)

network = ControlPowerFlow.ParserPWF.parse_file("data/24bus.pwf", add_control_data = true)

profile = Matrix(CSV.read("data/profile.csv", DataFrame))
profile = 1 .+ 0.2*(profile .- 1)

mixed_seasonality(t) = 1 # .+ rand(Uniform(-0.05, 0.05))*sin.(2*pi*collect(1:t)/24) + rand(Uniform(-0.05, 0.05))*cos.(2*pi*collect(1:t)/24)

noise() = rand(Uniform(0.98, 1.02))

steps_ahead, num_scen = size(profile)
contingency_length = 24*2

function handle_contingency_range(steps_ahead, contingency_length)
    start_cont = rand(collect(1:steps_ahead-contingency_length))
    end_cont = start_cont+contingency_length-1
    return start_cont:end_cont
end

contingency_range = [201:200+24*2]
contingency_range = [contingency_range[1] for i = 1:steps_ahead]
# ========================================================
#  Load and Generation Data
# ========================================================

int_br = ["33", "35", "36", "38", "39"]
int_load = ["4", "5", "6", "8"]
int_gen = ["22", "13"]
load_scenarios = Dict()
gen_scenarios = Dict()

for (l, v) in network["load"]
    load_scenarios[l] = network["load"][l]["pd"] 
    load_scenarios[l] *= l in int_load ? profile.*mixed_seasonality(h) : ones(size(profile))
end

for (g, v) in network["gen"]
    gen_scenarios[g] = network["gen"][g]["pg"]
    gen_scenarios[g] *= g in int_gen ? profile.*mixed_seasonality(h) : ones(size(profile))
end

# ========================================================
#  Bus and Branches of Interest
# ========================================================

buses = [4, 5, 6, 8]


contingency = "23"

# =====================================================
#  Results Initialization
# =====================================================

initialize_results() = Dict(
    "vm"       => Dict((b, zeros(steps_ahead,num_scen)) for (b, bus) in network["bus"]),
    "pf"       => Dict((b, zeros(steps_ahead,num_scen)) for (b, br)  in network["branch"]), 
    "pt"       => Dict((b, zeros(steps_ahead,num_scen)) for (b, br)  in network["branch"]), 
    "converge" => 0 .* Matrix{Bool}(undef,steps_ahead,num_scen),
    "data"     => Matrix{Dict}(undef,steps_ahead,num_scen),
    "time"     => zeros(steps_ahead,num_scen)
)

# ===================================
#  Analise dos cenários iniciais 
# ===================================

function apply_load_scenario!(data::Dict, load_scenarios::Dict, t::Int, s::Int)
    fp = 0.2
    for (l, v) in network["load"]
        data["load"][l]["pd"] = load_scenarios[l][t,s]*noise()
        data["load"][l]["qd"] = load_scenarios[l][t,s]*noise()*fp        
    end
end

function apply_gen_scenario!(data::Dict, gen_scenarios::Dict, t::Int, s::Int)
    for (g, v) in network["gen"]
        data["gen"][g]["pg"] *= gen_scenarios[g][t,s]*profile[t, s]*noise()
    end
end


function apply_scenario!(data::Dict, load_scenarios::Dict, gen_scenarios::Dict, t::Int, s::Int)
    apply_load_scenario!(data, load_scenarios, t, s)
    apply_gen_scenario!(data, gen_scenarios, t, s)
end

function apply_contingency!(data::Dict, contingency::String)
    data["branch"][contingency]["br_status"] = 0
end

function handle_results!(iter_data::Dict, results::Dict, t::Int, s::Int)
    for (b, bus) in iter_data["bus"]
        results["vm"][b][t, s] = bus["vm"]
    end
    for (b, branch) in iter_data["branch"]
        results["pf"][b][t, s] = branch["pf"]
        results["pt"][b][t, s] = branch["pt"]
    end
    results["data"][t,s] = iter_data
end

failed_converged_status(t, s) = @warn("AC Power Flow failed to converge at lead time $t of scenario $s")

pm_results = initialize_results()

# for s in 1:n_scenario, t in 1:n_time    
for s in 1:1, t in 1:steps_ahead  
    @info("Scenario $s, Lead Time $t")

    data = deepcopy(network)
    
    apply_scenario!(data, load_scenarios, gen_scenarios, t, s)
    
    t in contingency_range[s] ? apply_contingency!(data, contingency) : nothing

    pm_results["time"][t, s] = @elapsed result = run_ac_pf(data, ipopt)

    result["termination_status"] != MOI.LOCALLY_SOLVED ? failed_converged_status(t, s) : nothing

    update_data!(data, result["solution"])
    update_data!(data, calc_branch_flow_ac(data))
    
    handle_results!(data, pm_results, t, s)
end

lims = (0.85, 0.98)

plot(pm_results["vm"]["4"][:, 1], ylims = lims, label = "Bus 4", size = (900, 600), color = :1)
# plot!(pm_results["vm"]["5"][:, 1], ylims = lims, label = "Bus 5", color = :2)
plot!(pm_results["vm"]["6"][:, 1], ylims = lims, label = "Bus 6", color = :3)
# plot!(pm_results["vm"]["8"][:, 1], ylims = lims, label = "Bus 8", color = :4)

lims = (0.85, 1.05)
plot(pm_results["pf"]["34"][:, 1])
plot!(pm_results["pf"]["36"][:, 1])
plot!(pm_results["pf"]["37"][:, 1])


ct_results = initialize_results()

# for s in 1:n_scenario, t in 1:n_time    
for s in 1:1, t in 1:steps_ahead  
    @info("Scenario $s, Lead Time $t")

    data = deepcopy(network)
    data["info"] = Dict(
        "actions" => Dict(
            "csca" => true
            )
        )
    
    apply_scenario!(data, load_scenarios, gen_scenarios, t, s)
    
    if t in contingency_range[s] 
        apply_contingency!(data, contingency)
    
        ct_results["time"][t, s] = @elapsed begin
            pm = instantiate_model(data, ControlPowerFlow.ControlACPPowerModel, ControlPowerFlow.build_pf)
            set_optimizer(pm.model, ipopt)
            result = optimize_model!(pm)
        end
    else
        ct_results["time"][t, s] = @elapsed begin
            pm = instantiate_model(data, ACPPowerModel, PowerModels.build_pf)
            set_optimizer(pm.model, ipopt)
            result = optimize_model!(pm)
        end
    end    
    result["termination_status"] != MOI.LOCALLY_SOLVED ? failed_converged_status(t, s) : nothing

    update_data!(data, result["solution"])
    update_data!(data, calc_branch_flow_ac(data))
    
    handle_results!(data, ct_results, t, s)
end

lims = (0.85, 0.98)

plot!(ct_results["vm"]["4"][:, 1], ylims = lims, label = "Bus 4", size = (900, 600), color = :1, line = :dash)
# plot!(ct_results["vm"]["5"][:, 1], ylims = lims, label = "Bus 5", color = :2, line = :dash)
plot!(ct_results["vm"]["6"][:, 1], ylims = lims, label = "Bus 6", color = :3, line = :dash)
# plot!(ct_results["vm"]["8"][:, 1], ylims = lims, label = "Bus 8", color = :4, line = :dash)


plot(ct_results["pf"]["34"][:, 1])
plot!(ct_results["pf"]["36"][:, 1])
plot!(ct_results["pf"]["37"][:, 1])


# ==============================================
#  Identificação das divisões do sistema
# ==============================================

internal_bus = [1,2,3,4,5,6,7,8,9,10] 
boundary_bus = [11,12,24]
external_bus = [13,14,15,16,17,17,19,20,21,22,23]
interest_bus = union(internal_bus,boundary_bus)

df_branch = ControlActions.get_df_branch()

interest_br_fr = [(df_branch[h,1],df_branch[h,2],df_branch[h,3]) for h=1:length(df_branch[:,1]) if df_branch[h,1] in internal_bus || df_branch[h,2] in internal_bus]
interest_br_to = [(to,fr,c) for (fr,to,c) in interest_br_fr]
interest_br    = union(interest_br_fr,interest_br_to)
boundary_br    = [(3,24,1),(9,11,1),(9,12,1),(10,11,1),(10,12,1)]


# ===============================================
#  Analise dos cenários com ações corretivas 
# ===============================================

s0 = 0
# for s in 1:n_scenario, t in 1:n_time    
results["time"] = @elapsed for s in 1:1, t in 1:n_time
    OrganonDLL.open_file(@__DIR__, "untracked/" * case_file)   
    s != s0 ? println("Analisando o cenário $s...") : Nothing
    s0 = s

    # if t == 200
    if 181 ≤ t ≤ 230
        OrganonDLL.reload_branchstatus(4,9,1,0)
    #elseif t == 221
    else
        OrganonDLL.reload_branchstatus(4,9,1,1)
    end
    
    for h in 1:length(load_names)
        OrganonDLL.change_activeload(load_names[h], load_values[t,h,s])
        OrganonDLL.change_reactiveload(load_names[h], load_values[t,h,s] * fc[load_names[h]])
    end
    for h in 1:length(gen_names)
        OrganonDLL.change_activegen(gen_names[h], gen_values[t,h,s])
    end

    om, converge       = run_organon_pf()
    approved, hist, om = operational_analysis(violation_list, action_list, om; reverse = false, silence = true)
    results["history"][t,s]  = hist
    results["approved"][t,s] = approved
    results["converge"][t,s] = converge

    for (fr,to,c) in interest_br_fr
        br_id = findall(x->x==(fr,to,c), om.pf[:,1])
        results["p"][(fr,to,c)][t,s] = om.pf[br_id,2][1]
        results["p"][(to,fr,c)][t,s] = om.pf[br_id,4][1]
    end
    for h in 1:length(interest_bus)
        bus = interest_bus[h]
        results["vm"][bus][t,s] = om.voltage[bus,2]
    end
end

# =============================
#  Anásele dos resultados      
# =============================
plot_p = empt_plot00("Ações Corretivas", "Fluxo [Mw]")
plot_P!(24,3 ,1,1,5,1.)
plot_inicial_P!(24,3 ,1,1,5,1.0)
plot_P!(11,9 ,1,1,8,1.)
plot_inicial_P!(11,9 ,1,1,8,1.0)
plot_P!(11,10 ,1,1,7,1.)
plot_inicial_P!(11,10 ,1,1,7,1.0)
plot_P!(12,9 ,1,1,6,1.)
plot_inicial_P!(12,9 ,1,1,6,1.0)
plot_P!(12,10 ,1,1,9,1.)
plot_inicial_P!(12,10 ,1,1,9,1.0)

plot_vm = empt_plot00("", "Tensão [pu]")
plot_VM!(8,1,1,1.)
plot_inicial_VM!(8,1,1,1.0)
plot_VM!(6,1,2,1.)
plot_inicial_VM!(6,1,2,1.0)
plot_VM!(5,1,3,1.)
plot_inicial_VM!(5,1,3,1.0)
plot_VM!(4,1,4,1.)
plot_inicial_VM!(4,1,4,1.0)

plot_sh = plot_shunt(1,10,1.;ls=:solid)

plt_shunt = plot(plot_p, plot_vm, plot_sh;
    titlefontsize=12, xtickfontsize=10,ytickfontsize=10, legendfontsize=10, yguidefontsize=10, xguidefontsize=11,
    layout = grid(3, 1, heights=[0.53, 0.33, 0.14]),
    bottom_margin = 0Plots.mm,
    left_margin = 5Plots.mm,
    size   = (800*1.1,600*0.9),
    dpi    = 300,
    legend = :outertopright
)
# gui(plt)
# savefig("untracked/shunt_control")

#######################################################
