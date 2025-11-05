module PeriodicARModels

include("structure.jl")

export MonthlyAR, Multi_MonthlyAR, rand, fit_AR, fit_Multi_AR, save_model, load_model, CaracteristicsSeries, init_CaracteristicsSeries, extract_series, Common_indexes, dayofyear_Leap

#for testing : 
export LOESS, fitted_periodicity_fonc


end
