using PeriodicARModels, OrderedCollections, FileIO, JLD2, Random
using Test

##Seed for random reproducibility
Random.seed!(1234)

##Station
file_TN = joinpath(@__DIR__, "..", "stations", "TN_Nantes.txt")
file_TX = joinpath(@__DIR__, "..", "stations", "TX_Nantes.txt")

##### UNIVARIATE AR MODEL #####

##Model hyperparameters  
p = 2
method_ = "monthlyLL"                 # "mean", "median", "concat", "sumLL", "monthlyLL"
periodicity_model = "trigo"           # "trigo", "smooth", "autotrigo", "stepwise_trigo"
degree_period = 0                     # 0 => default value -> "trigo" : 5, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
Trendtype = "LOESS"                   # "LOESS", "polynomial", "null" (for no additive trend)
trendparam = nothing                  # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1
σ_periodicity_model = "trigo"         # "trigo", "smooth", "autotrigo", "stepwise_trigo", "null" (for no multiplicative periodicity)
σ_degree_period = 0                   # 0 => default value -> "trigo" : 5, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
σ_Trendtype = "LOESS"                 # "LOESS", "polynomial", "null" (for no multiplicative trend)
σ_trendparam = nothing                # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1

##Simulations
n = 3

##All settings
settings_uni = OrderedDict((("file", file_TN),
    ("p", p),
    ("method_", method_),
    ("periodicity_model", periodicity_model),
    ("degree_period", degree_period),
    ("Trendtype", Trendtype),
    ("trendparam", trendparam),
    ("σ_periodicity_model", σ_periodicity_model),
    ("σ_degree_period", σ_degree_period),
    ("σ_Trendtype", σ_Trendtype),
    ("σ_trendparam", σ_trendparam),
    ("n", n)))

series_uni = first(extract_series(file_TN), 2000)

model_uni = fit_AR(series_uni[:, 2], series_uni.DATE,
    p=p,
    method_=method_,
    periodicity_model=periodicity_model,
    degree_period=degree_period,
    Trendtype=Trendtype,
    trendparam=trendparam,
    σ_periodicity_model=σ_periodicity_model,
    σ_degree_period=σ_degree_period,
    σ_Trendtype=σ_Trendtype,
    σ_trendparam=σ_trendparam)

sample_uni = rand(model_uni, n, series_uni.DATE)

##### MULTIVARIATE AR MODEL #####

##Model hyperparameters
p = 2
method_ = "monthly"
periodicity_model = "trigo"       # "trigo", "smooth", "autotrigo", "stepwise_trigo"
degree_period = 2                 # 0 => default value -> "trigo" : 5, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
Trendtype = "LOESS"               # "LOESS", "polynomial", "null" (for no additive trend)
trendparam = 0.16                 # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1
σ_periodicity_model = "trigo"     # "trigo", "smooth", "autotrigo", "stepwise_trigo", "null" (for no multiplicative periodicity)
σ_degree_period = 2               # 0 => default value -> "trigo" : 5, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
σ_Trendtype = "LOESS"             # "LOESS", "polynomial", "null" (for no multiplicative trend)
σ_trendparam = 0.16               # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1

##Simulations
n = 3

settings_multi = OrderedDict((("file", file_TN),
    ("p", p),
    ("method_", method_),
    ("periodicity_model", periodicity_model),
    ("degree_period", degree_period),
    ("Trendtype", Trendtype),
    ("trendparam", trendparam),
    ("σ_periodicity_model", σ_periodicity_model),
    ("σ_degree_period", σ_degree_period),
    ("σ_Trendtype", σ_Trendtype),
    ("σ_trendparam", σ_trendparam),
    ("n", n)))

date_vec_multi, x_multi = Common_indexes(file_TN, file_TX)

model_multi = fit_Multi_AR(x_multi[1:2000,:], date_vec_multi[1:2000],
    p=p,
    method_=method_,
    periodicity_model=periodicity_model,
    degree_period=degree_period,
    Trendtype=Trendtype,
    trendparam=trendparam,
    σ_periodicity_model=σ_periodicity_model,
    σ_degree_period=σ_degree_period,
    σ_Trendtype=σ_Trendtype,
    σ_trendparam=σ_trendparam)

sample_multi = rand(model_multi, n)

save("test/references.jld2", "ref_data", (model_uni, sample_uni, model_multi, sample_multi))

