
function run_powermodels(file::String)
    network = PWF.parse_file(file, pm = true)
    result = run_ac_pf(network, IPOPT)
    update_data!(network, result["solution"])
    update_data!(network, calc_branch_flow_ac(network))
    return network
end


## Funções auxiliares para ler os resultados tabulares do ANAREDE e 
## comparar com os resultados do PowerModels.jl
using Statistics

float(str) = parse(Float64, replace(str, ","=>"."))

function filter_bus(bus::Dict)
    bus = deepcopy(bus)
    remain_fields =  ["vm", "va", "area", "bus_i", "bus_type"]
    for (k, v) in bus
        k in remain_fields ? nothing : delete!(bus, k)    
    end
    return bus
end

function filter_branch(branch::Dict)
    branch = deepcopy(branch)
    remain_fields = ["f_bus", "t_bus", "br_status", "pf", "qf", "pt", "qt"]
    for (k, v) in branch
        k in remain_fields ? nothing : delete!(branch, k)    
    end
    return branch
end

function initiate_dict(network::Dict)
    return Dict(
        "bus" => Dict((k, filter_bus(v)) for (k, v) in network["bus"]),
        "branch" => Dict((k, filter_branch(v)) for (k, v) in network["branch"])
    )
end

function _handle_rbar(rbar::DataFrame, anarede_dict::Dict)
    n_bus = size(rbar, 1)
    for i = 1:n_bus
        num = rbar[i, 1]
        bus = anarede_dict["bus"]["$num"]
        bus["bus_i"] = num
        bus["vm"] = float(rbar[i, 3])
        bus["va"] = float(rbar[i, 4]) * pi / 180
    end
end

function get_branch(from::Int, to::Int, circuit::Int, branches::Dict)
    br_id = findall(
        branch -> (branch["f_bus"] == from && branch["t_bus"] == to), branches
    )
    if !isempty(br_id)
        return br_id[circuit], true
    end
    br_id = findall(
        branch -> (branch["f_bus"] == to && branch["t_bus"] == from), branches
    )[circuit]
    return br_id, false
end

function _handle_rlin(rlin::DataFrame, anarede_dict::Dict)
    n_branch = size(rlin, 1)
    for i = 1:n_branch
        from    = rlin[i, 1]
        to      = rlin[i, 18]
        circuit = rlin[i, 21]
        
        br_id, orientation = get_branch(from, to, circuit, anarede_dict["branch"])
        
        branch = anarede_dict["branch"][br_id]
        pf     = float(rlin[i, 22]) / 100
        qf     = float(rlin[i, 23]) / 100
        p_loss = float(rlin[i, 33]) / 100
        q_loss = float(rlin[i, 34]) / 100

        branch["pf"] = orientation ? pf : -pf + p_loss
        branch["qf"] = orientation ? qf : -qf + q_loss
        branch["pt"] = orientation ? -pf + p_loss : pf
        branch["qt"] = orientation ? -qf + q_loss : qf
    end
end

function anarede_results(rbar::DataFrame, rlin::DataFrame, network::Dict)
    anarede_dict = initiate_dict(network)
    _handle_rbar(rbar, anarede_dict)
    _handle_rlin(rlin, anarede_dict)
    return anarede_dict
end

## Generate Statistics from 2 results dictionaries

initiate_field(len) = Dict(
    "error" => Matrix(undef, len, 2),
    "max" => NaN,
    "min" => NaN,
    "mean" => NaN,
    "quantiles" => Dict(),    
)

function initiate_error(network) 
    n_bus    = length(keys(network["bus"]))
    n_branch = length(keys(network["branch"]))
    error = Dict(
        "bus" => Dict(
            "vm" => initiate_field(n_bus),
            "va" => initiate_field(n_bus),
        ),
        "branch" => Dict(
            "pf" => initiate_field(n_branch),
            "qf" => initiate_field(n_branch),
            "pt" => initiate_field(n_branch),
            "qt" => initiate_field(n_branch),
        )
    )
    return error
end

function calc_error(net1::Dict, net2::Dict)
    err = initiate_error(net1)
    for (i, bus) in enumerate(sort(parse.(Int, [k for (k, v) in  net1["bus"]])))
        b = "$bus"
        err["bus"]["vm"]["error"][i, 1] = b
        err["bus"]["vm"]["error"][i, 2] = net1["bus"][b]["vm"] - net2["bus"][b]["vm"]
        err["bus"]["va"]["error"][i, 1] = b
        err["bus"]["va"]["error"][i, 2] = net1["bus"][b]["va"] - net2["bus"][b]["va"]
    end
    for (i, branch) in enumerate(sort(parse.(Int, [k for (k, v) in  net1["branch"]])))
        br = "$branch"
        err["branch"]["pf"]["error"][i, 1] = br
        err["branch"]["pf"]["error"][i, 2] = net1["branch"][br]["pf"] - net2["branch"][br]["pf"]
        err["branch"]["pt"]["error"][i, 1] = br
        err["branch"]["pt"]["error"][i, 2] = net1["branch"][br]["pt"] - net2["branch"][br]["pt"]
        err["branch"]["qf"]["error"][i, 1] = br
        err["branch"]["qf"]["error"][i, 2] = net1["branch"][br]["qf"] - net2["branch"][br]["qf"]
        err["branch"]["qt"]["error"][i, 1] = br
        err["branch"]["qt"]["error"][i, 2] = net1["branch"][br]["qt"] - net2["branch"][br]["qt"]
    end
    return err
end

function handle_statistics(field::Dict)
    quantiles_probs = [0.95, 0.75, 0.5, 0.25, 0.05]
    err = field["error"][:, 2]
    err[isnan.(err)] .= 0
    field["max"] = maximum(err)
    field["min"] = minimum(err)
    field["mean"] = mean(err)
    quantiles = quantile(err, quantiles_probs)
    field["quantiles"] = Dict(quantiles_probs[i] => quantiles[i] for i = 1:length(quantiles_probs))
end

function error_statistics(net1, net2)
    # verify_networks!(net1, net2)
    err = calc_error(net1, net2)
    handle_statistics(err["bus"]["vm"])
    handle_statistics(err["bus"]["va"])
    handle_statistics(err["branch"]["pf"])
    handle_statistics(err["branch"]["qf"])
    handle_statistics(err["branch"]["pt"])
    handle_statistics(err["branch"]["qt"])
    return err
end

get_statistics(field) = [
    field["min"],
    field["quantiles"][0.05],
    field["quantiles"][0.25],
    field["mean"],
    field["quantiles"][0.5],
    field["quantiles"][0.75],
    field["quantiles"][0.95],
    field["max"]
]

function create_statistics_table(err::Dict)
    # name | min | 0.05 | 0.25 | mean | median | 0.75 | 0.95 | max 
    #  vm
    #  va
    #  pf
    #  qf
 
    table = Matrix{Any}(undef, 5, 9)
    table[1, :] = ["Field", "Min", "Q05", "Q25", "̄Mean", "Median", "Q75", "Q95", "Max"]
    table[2, 1] = "vm"
    table[3, 1] = "va"
    table[4, 1] = "pf"
    table[5, 1] = "qf"
    table[2, 2:end] = get_statistics(err["bus"]["vm"])
    table[3, 2:end] = get_statistics(err["bus"]["va"])
    table[4, 2:end] = get_statistics(err["branch"]["pf"])
    table[5, 2:end] = get_statistics(err["branch"]["qf"])
    return DataFrame(table[2:end,:], String.(table[1,:]))
end