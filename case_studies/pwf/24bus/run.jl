using CSV, DataFrames
using PowerModels, Ipopt, JuMP, StatsPlots
const IPOPT = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-4)

# Include parser
include("C:\\Users\\iagoc\\Dropbox\\Lamps PUC\\Energisa Must II\\Packages\\PWF.jl\\src/PWF.jl")

# Include auxiliar functions
include("../functions.jl")

# reading anarede tabular results
rbar = CSV.read("case_studies/pwf/24bus/data/rbar.csv", DataFrame)
rlin = CSV.read("case_studies/pwf/24bus/data/rlin.csv", DataFrame)

file = "case_studies/pwf/24bus/data/24bus.pwf"

# running PowerModels
network = run_powermodels(file)

# converting tabular results into dictionary
anarede = anarede_results(rbar, rlin, network)

# calculating error statistics
errors = error_statistics(network, anarede)

############################
#    Exporting results     #
############################

# Export densisty of error 
vm_err = errors["bus"]["vm"]["error"][:, 2]
p_dens_vm = density(vm_err, label = "vm mismatch", xaxis = "p.u.", yaxis = "frequency", size = (700, 400))
savefig(p_dens_vm, "case_studies/pwf/24bus/results/vm_density.pdf")

va_err = errors["bus"]["va"]["error"][:, 2]
p_dens_va = density(va_err, label = "va mismatch", xaxis = "p.u.", yaxis = "frequency", size = (700, 400))
savefig(p_dens_va, "case_studies/pwf/24bus/results/va_density.pdf")

pf_err = errors["branch"]["pf"]["error"][:, 2]
p_dens_pf = density(pf_err[isnan.(pf_err) .== false], label = "pf mismatch", xaxis = "p.u.", yaxis = "frequency", size = (700, 400))
savefig(p_dens_pf, "case_studies/pwf/24bus/results/pf_density.pdf")

qf_err = errors["branch"]["qf"]["error"][:, 2]
p_dens_qf = density(qf_err[isnan.(qf_err) .== false], label = "qf mismatch", xaxis = "p.u.", yaxis = "frequency", size = (700, 400))
savefig(p_dens_qf, "case_studies/pwf/24bus/results/qf_density.pdf")

# Export Statistical Table
CSV.write("case_studies/pwf/24bus/results/error.csv", create_statistics_table(errors))
