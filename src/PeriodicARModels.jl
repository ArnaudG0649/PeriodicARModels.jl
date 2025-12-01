module PeriodicARModels

include("structure.jl")

export MonthlySWG, MonthlyAR, rand, save_model, load_model, CaracteristicsSeries, init_CaracteristicsSeries, extract_series, Common_indexes, dayofyear_Leap, decompose

#for testing : 
export LOESS, fitted_periodicity_fonc


end
