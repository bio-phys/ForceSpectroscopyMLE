# ForceSpectroscopyMLE

This package provides a robust framework to analyze rupture force data from single-molecule force spectroscopy experiments.  It includes a systematic protocol for trimming unwanted outliers and an efficient maximum likelihood estimator, based on the Dudko-Hummer-Szabo (DHS) bond rupture model, to extract parameters characterizing the free-energy landscape of the bond and the force-free disassociation rate.  

For more details on the theoretical framework, please refer to the associated publication:
> W. Cai, J. T. Bullerjahn, M. Lallemang, K. Kroy, B. N. Balzer, and T. Hugel, "Angle-dependent strength of a single chemical bond by stereographic force spectroscopy", *Chemical Science* **13**, 5734-5740 (2022). https://doi.org/10.1039/D2SC01077A

The code makes use of the DHS model of forcible bond rupture, which was originally published in:
> O. K. Dudko, G. Hummer, and A. Szabo, "Intrinsic rates and activation free energies from single-molecule pulling experiments", *Physical Review Letters* **96**, 108101 (2006). https://doi.org/10.1103/PhysRevLett.96.108101

Please cite the references above if you use `ForceSpectroscopyMLE` to analyze your data.  



## Try without installation

Click this [link](https://bio-phys.pages.mpcdf.de/forcespectroscopymle) to run the `pipeline_example.ipynb` notebook in a cloud environment.  

You can either analyze the data sets found in `examples/mock_data/` or upload your own.  The interactive session is only temporary and files will be deleted after termination (File -> Shut Down).  

**Important:** Please shut down JupyterLab properly after use via the drop-down menu (File -> Shut Down) to free resources for other users.  



## Installation

The package is written in the open-source programming language [Julia](https://github.com/JuliaLang/julia), which can be downloaded from their [webpage](https://julialang.org/downloads/#download_julia).  

Currently, the package is not in a registry.  It must therefore be added by specifying a URL to the repository:
```julia
using Pkg; Pkg.add(url="https://github.com/bio-phys/ForceSpectroscopyMLE")
```
Users of older versions of Julia may need to wrap the contents of the brackets with `PackageSpec()`.  



## Usage

### Importing data

The rupture force data should be of the type `Array{Float64,2}`, where the first column contains the rupture forces `F` (in *pN*) and the second column the associated loading rates `dF` (in *pN/s*).  In principle, users can lump all their measured force spectra into a single file and, e.g., read it in as follows:
```julia
using DelimitedFiles

data = readdlm(file_name)
```
However, in order to make use of our data trimming protocol, we recommend keeping data measured at different pulling speeds in separate files (stored in the directory `rupture_forces`), which can be read in using our specialized function:
```julia
using ForceSpectroscopyMLE

data = read_data("./rupture_forces/")
```
The array `data` is then of the type `Array{Array{Float64,2},1}`.  



### Parameter and error estimation

We can estimate the parameters `βΔG_u`, `x_u` and `k_0` of the DHS model using the `MLE_estimator` function:
```julia
all_data = vcat(data...) # only necessary if 'data' is of the type Array{Array{Float64,2},1}
parameters = MLE_estimator(all_data,ν) # βΔG_u, x_u, k_0
```
The parameter `ν` can be set to `1/2` or `2/3` depending on the shape of the underlying free-energy landscape.  For `ν = 1` the DHS model reduces to the Bell-Evans model, which only depends on the parameters `x_u` and `k_0`.  `MLE_estimator` has various optional arguments, most of which are inputs for the [optimizer](https://github.com/robertfeldt/BlackBoxOptim.jl) except for the absolute temperature `T` (in *K*):
```julia
MLE_estimator(all_data,ν,T=295,βΔE_range=(0.1,100.0),Δx_b_range=(0.001,10.0),msteps=100000,mode=:compact,psize=50,tint=60.0)
```
The `MLE_errors` function provides an estimate of the parameter uncertainties:
```julia
errors = MLE_errors(all_data,ν) # δβΔG_u, δx_u, δk_0
```
with (almost) the same optional arguments as `MLE_estimator`:
```julia
MLE_errors(all_data,ν,N=100,T=295,βΔE_range=(0.1,100.0),Δx_b_range=(0.001,10.0),msteps=100000,mode=:silent,psize=50,tint=60.0)
```
We rely on bootstrapping to gauge the uncertainty of the estimates, by generating `N` new data sets from our sample of rupture forces and analyzing the results.  This can become rather sluggish for large `N`, so it is recommended to run the command `export JULIA_NUM_THREADS=n`, with `n` being the number of available (physical) cores, before launching Julia.  This speeds up the numerics significantly.  

To check the number of available cores for threading, simply run
```julia
using Base.Threads; nthreads()
```
This should print the number `n` if the above-mentioned command was executed properly.  



### Data trimming

The function `read_data` sorts the data sets in ascending order with respect to the rupture forces.  We can therefore use `reduce_data` to trim the last `i` datapoints from each data set, resulting in a reduced data set:
```julia
reduced_data = reduce_data(data,i)
```
For comparison, we can also randomly remove `i` datapoints from each data set:
```julia
randomly_reduced_data = random_reduce_data(data,i)
```
A more detailed example that systematically investigates the effect of data trimming on the parameter estimates can be found in [the examples directory](examples).  
