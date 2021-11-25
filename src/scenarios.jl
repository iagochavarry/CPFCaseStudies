using CSV, Dates, DataFrames, Plots, SeasonalTrendLoess, StateSpaceModels, Pipe

profile_history = CSV.read("data/profile_history.csv", DataFrame; delim = ";")[1:end-1, :]

function treat_profile(profile::DataFrame)
    timestamps = profile[:, 1]
    values = profile[:, 2]
    timestamps = DateTime.(timestamps, "dd/mm/yyyy HH:MM")
    values = [@pipe replace(val, "," => ".") |> parse.(Float64, _) for val in values]

    return DataFrame(timestamps = timestamps, values = values)
end

prof_hist = treat_profile(profile_history)

# historical data
y = prof_hist.values
# Weeks forecasted
w = 2
# Lead Times
h = 24*7*w
# Seasonal Cycle Duration
s = 24*7
# Number of scenarios
ω = 200

# Forecast timestamps
forecast_timestamps = [prof_hist.timestamps[end] + Hour(i) for i = 1:h]

# Seasonal Decomposition by STL
stl_decomp = stl(y, s)
# Seasonal Pattern
week_cycle = stl_decomp.seasonal[end-s+1:end]
# Remaider + Trend
y_tr = stl_decomp.remainder + stl_decomp.trend

# Probabilistic forecast
model_ll = LocalLevel(y_tr)
fit!(model_ll)

# Generating Scenarios
scen_ll = simulate_scenarios(model_ll, h, ω)[:, 1, :]

# Seasonal Forecast -> Repeat last week cycle in the data
seas_fore = []
for i in 1:w 
    seas_fore = vcat(seas_fore, week_cycle)
end

# Generating Profile with Seasonal Pattern
scen_y = scen_ll .+ seas_fore

hist = 24*7*3
plot(prof_hist.timestamps[end-hist+1:end], prof_hist.values[end-hist+1:end], size = (1200, 700))
plot!(forecast_timestamps, scen_y, legend = false)

profile = scen_y
# Load profile forecast
for scenario in 1:ω 
    mean = sum(scen_y[:, scenario])/length(scen_y[:, scenario])
    profile[:, scenario] = scen_y[:, scenario] ./ mean
end
# creating a DataFrame to export
profile = DataFrame(profile, :auto)
# export to .csv
CSV.write("data/profile.csv", profile)


