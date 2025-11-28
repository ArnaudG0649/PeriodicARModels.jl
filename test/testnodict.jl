##Seed for random reproducibility
Random.seed!(1234)

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

z, trend, period, σ_trend, σ_period = decompose(series_uni[:, 2],
    series_uni.DATE,
    periodicity_model,
    degree_period,
    Trendtype,
    trendparam,
    σ_periodicity_model,
    σ_degree_period,
    σ_Trendtype,
    σ_trendparam)

model_uni = fit_AR(z, series_uni.DATE, p=p, method_=method_)

periodicity = period[dayofyear_Leap.(series_uni.DATE)]
σ_periodicity = σ_period[dayofyear_Leap.(series_uni.DATE)]
nspart = trend .+ periodicity
σ_nspart = σ_trend .* σ_periodicity

sample_uni_nd = [rand(model_uni, month.(series_uni.DATE), y₁=z[1:p]) for _ in 1:n]
sample_uni_nd = map(sim -> sim .* σ_nspart + nspart, sample_uni_nd)

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

z, trend, period, σ_trend, σ_period = decompose(x_multi,
    date_vec_multi,
    periodicity_model,
    degree_period,
    Trendtype,
    trendparam,
    σ_periodicity_model,
    σ_degree_period,
    σ_Trendtype,
    σ_trendparam)

model_multi = fit_Multi_AR(z, date_vec_multi, p=p, method_=method_)

periodicity = period[dayofyear_Leap.(date_vec_multi), :]
σ_periodicity = σ_period[dayofyear_Leap.(date_vec_multi), :]
nspart = trend .+ periodicity
σ_nspart = σ_trend .* σ_periodicity

sample_multi_nd = [rand(model_multi, month.(date_vec_multi), y₁=z[1:p, :], nspart=nspart, σ_nspart=σ_nspart) for _ in 1:n]















# x_ = series_uni[:, 2]
# trend = LOESS(x_, trendparam)
# y_ = x_ - trend
# trigo_function = fitted_periodicity_fonc(y_, series_uni.DATE, OrderTrig=degree_period)
# periodicity, period = trigo_function.(series_uni.DATE), trigo_function.(Date(0):(Date(1)-Day(1)))
