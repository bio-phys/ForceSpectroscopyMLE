# Code for estimating free-energy landscape parameters and dissociation rate from single-molecule force spectroscopy data
# version 1.0 (28/09/2021)
# Jakob Tómas Bullerjahn (jabuller@biophys.mpg.de)
# Wanhao Cai

# Please read and cite the associated publications: 
# W. Cai, J. T. Bullerjahn, M. Lallemang, K. Kroy, B. N. Balzer, and T. Hugel, "Direction dependence of bond strength and polymer chain elasticity", submitted. 
# O. K. Dudko, G. Hummer, and A. Szabo, "Intrinsic rates and activation free energies from single-molecule pulling experiments", Phys. Rev. Lett. 96, 108101 (2006). 



module ForceSpectroscopyMLE

export
        # MLE:
        MLE_estimator, 
        MLE_errors, 

        # SMFS:
        meanF, 
        varF, 
        D, 

        # Data manipulation:
        read_data, 
        read_data_unsorted, 
        reduce_data, 
        random_reduce_data, 
        data_exclusion, 
        random_data_exclusion

using Base.Threads
using BlackBoxOptim
using DelimitedFiles
using ProgressMeter
using SpecialFunctions
using Statistics
using StatsBase

include("data_functions.jl")
include("smfs_functions.jl")
include("mle_functions.jl")

end
