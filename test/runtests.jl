using PeriodicARModels, OrderedCollections, FileIO, JLD2, Random, Dates
using Test

##Seed for random reproducibility
Random.seed!(1234)

GetAllAttributes(object) = map(field -> getfield(object, field), fieldnames(typeof(object)))
## Source : https://discourse.julialang.org/t/get-the-name-and-the-value-of-every-field-for-an-object/87052/2

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

period = model_uni.period[dayofyear_Leap.(series_uni.DATE)]
σ_period = model_uni.σ_period[dayofyear_Leap.(series_uni.DATE)]
nspart = model_uni.trend .+ period
σ_nspart = model_uni.σ_trend .* σ_period
sample_uni = [rand(model_uni, month.(series_uni.DATE),y₁=model_uni.z[1:p]) for _ in 1:n]
sample_uni = map(sim -> sim .* σ_nspart + nspart, sample_uni)
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

date_vec_multi, x_multi = Common_indexes(file_TN, file_TX)
date_vec_multi, x_multi = date_vec_multi[1:2000], x_multi[1:2000, :]

model_multi = fit_Multi_AR(x_multi, date_vec_multi,
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

period = model_multi.period[dayofyear_Leap.(date_vec_multi), :]
σ_period = model_multi.σ_period[dayofyear_Leap.(date_vec_multi), :]
nspart = model_multi.trend .+ period
σ_nspart = model_multi.σ_trend .* σ_period
using PeriodicARModels
sample_multi = [rand(model_multi, month.(date_vec_multi), y₁=model_multi.z[1:p,:], nspart=nspart, σ_nspart=σ_nspart) for _ in 1:n]

ref_data = load(joinpath(@__DIR__, "references.jld2"))["ref_data"]

@testset "PeriodicARModels.jl" begin
    # @test GetAllAttributes(model_uni) == GetAllAttributes(ref_data[1])
    @test sample_uni == ref_data[2]
    # @test GetAllAttributes(model_multi) == GetAllAttributes(ref_data[3])
    @test sample_multi == ref_data[4]

end

















# x_ = series_uni[:, 2]
# trend = LOESS(x_, trendparam)
# y_ = x_ - trend
# trigo_function = fitted_periodicity_fonc(y_, series_uni.DATE, OrderTrig=degree_period)
# periodicity, period = trigo_function.(series_uni.DATE), trigo_function.(Date(0):(Date(1)-Day(1)))
