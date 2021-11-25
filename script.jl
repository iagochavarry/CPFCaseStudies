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

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 0.01);

network = ControlPowerFlow.ParserPWF.parse_file("data/5bus_csca.pwf"; software = ControlPowerFlow.ParserPWF.ANAREDE);

pm = ControlPowerFlow._PM.instantiate_model(network, ACPPowerModel, ControlPowerFlow.build_pf);

set_optimizer(pm.model, ipopt)
result = optimize_model!(pm)
update_data!(network, result["solution"])
PowerModels.print_summary(network)

network = ControlPowerFlow.ParserPWF.parse_file("data/5bus_csca_cont.pwf"; software = ControlPowerFlow.ParserPWF.ANAREDE);

pm = ControlPowerFlow._PM.instantiate_model(network, ACPPowerModel, ControlPowerFlow.build_pf);

set_optimizer(pm.model, ipopt)
result = optimize_model!(pm)
update_data!(network, result["solution"])
PowerModels.print_summary(network)

network = ControlPowerFlow.ParserPWF.parse_file("data/5bus_csca_cont.pwf"; software = ControlPowerFlow.ParserPWF.ANAREDE, add_control_data = true);
network["info"] = Dict("actions"=>Dict("csca"=>true))

pm = ControlPowerFlow._PM.instantiate_model(network, ControlPowerFlow.ControlACPPowerModel, ControlPowerFlow.build_pf);
set_optimizer(pm.model, ipopt)
result = optimize_model!(pm)
update_data!(network, result["solution"])
PowerModels.print_summary(network)


print(pm.model)

set_optimizer(pm.model, ipopt);
result = ControlPowerFlow._PM.optimize_model!(pm)


result["solution"]["bus"]






