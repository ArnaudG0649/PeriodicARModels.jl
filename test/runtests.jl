using PeriodicARModels, OrderedCollections, FileIO, JLD2, Random, Dates
using Test

##Seed for random reproducibility
Random.seed!(1234)

GetAllAttributes(object) = map(field -> getfield(object, field), fieldnames(typeof(object)))
## Source : https://discourse.julialang.org/t/get-the-name-and-the-value-of-every-field-for-an-object/87052/2

#For evaluation of sample and estimated paramaters
α=1E-4
AreClose(x::AbstractArray{T},y::AbstractArray{T},α=α) where T <: AbstractFloat = maximum(abs.(x .- y)) < α
AreClose(x::AbstractArray{T},y::AbstractArray{T},α=α) where T <: AbstractArray{<:AbstractFloat} = AreClose(vcat(x...),vcat(y...),α)
FlatΦ(Φ,d=2) = vcat([hcat([reshape(v_sub,(1,d^2)) for v_sub in Φ[m]]...) for m in 1:12]...)
AreClose(x::AbstractVector{T},y::AbstractVector{T},α=α) where T <: AbstractVector{<:AbstractMatrix} = AreClose(FlatΦ(x),FlatΦ(y),α)


# AreClose(model_uni.monthlyAR.Φ, ref_data[1].Φ)
# AreClose(model_uni.monthlyAR.σ, ref_data[1].σ)
# AreClose(sample_uni, ref_data[2])

# AreClose(model_multi.monthlyAR.Φ, ref_data[3].Φ)
# AreClose(model_multi.monthlyAR.σ, ref_data[3].σ)
# AreClose(sample_multi, ref_data[4])


##Station
file_TN = joinpath(@__DIR__, "..", "stations", "TN_Nantes.txt")
file_TX = joinpath(@__DIR__, "..", "stations", "TX_Nantes.txt")

##### UNIVARIATE AR MODEL #####

##Model hyperparameters
p = 2
method_ = "monthlyLL"                 # "mean", "median", "concat", "sumLL", "monthlyLL"
periodicity_model = "trigo"           # "trigo", "smooth", "autotrigo", "stepwise_trigo"
degree_period = 0                     # 0 => default value -> "trigo" : 8, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
Trendtype = "LOESS"                   # "LOESS", "polynomial", "null" (for no additive trend)
trendparam = nothing                  # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1
σ_periodicity_model = "trigo"         # "trigo", "smooth", "autotrigo", "stepwise_trigo", "null" (for no multiplicative periodicity)
σ_degree_period = 0                   # 0 => default value -> "trigo" : 8, "smooth" : 9, "autotrigo" : 50, "stepwise_trigo" : 50
σ_Trendtype = "LOESS"                 # "LOESS", "polynomial", "null" (for no multiplicative trend)
σ_trendparam = nothing                # nothing => default value -> "LOESS" : 0.08, "polynomial" : 1

##Simulations
n = 3

series_uni = first(extract_series(file_TN), 2000)

model_uni = MonthlySWG(series_uni[:, 2], series_uni.DATE,
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

sample_uni = rand(model_uni, y₁=model_uni.z[1:p], n_sim=n)


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

model_multi = MonthlySWG(x_multi, date_vec_multi,
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


sample_multi = rand(model_multi, y₁=model_multi.z[1:p, :], n_sim=n)

# include(joinpath(@__DIR__, "testnodict.jl"))

ref_data = load(joinpath(@__DIR__, "references.jld2"))["ref_data"]

@testset "PeriodicARModels.jl" begin
    @test AreClose(model_uni.monthlyAR.Φ, ref_data[1].Φ)
    @test AreClose(model_uni.monthlyAR.σ, ref_data[1].σ)
    @test AreClose(sample_uni, ref_data[2])
    rand(model_uni, n_sim=2)
    rand(model_uni)

    @test AreClose(model_multi.monthlyAR.Φ, ref_data[3].Φ)
    @test AreClose(model_multi.monthlyAR.σ, ref_data[3].σ)
    @test AreClose(sample_multi, ref_data[4])
    rand(model_multi, n_sim=2)
    rand(model_multi)

end





# x_ = series_uni[:, 2]
# trend = LOESS(x_, trendparam)
# y_ = x_ - trend
# trigo_function = fitted_periodicity_fonc(y_, series_uni.DATE, OrderTrig=degree_period)
# periodicity, period = trigo_function.(series_uni.DATE), trigo_function.(Date(0):(Date(1)-Day(1)))
