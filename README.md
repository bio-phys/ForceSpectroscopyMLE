# ForceSpectroscopyMLE

This package provides a robust framework to analyze rupture force data from single-molecule force spectroscopy experiments.  It includes a systematic protocol for trimming unwanted outliers and an efficient maximum likelihood estimator, based on the Dudko-Hummer-Szabo (DHS) bond rupture model, to extract parameters characterizing the free-energy landscape of the bond and the force-free disassociation rate.  

For more details on the theoretical framework, please refer to the associated publication:
> W. Cai, J. T. Bullerjahn, M. Lallemang, K. Kroy, B. N. Balzer, and T. Hugel, "Direction dependence of bond strength and polymer chain elasticity", submitted. 

The code makes use of the DHS model of forcible bond rupture, which was originally published in
> O. K. Dudko, G. Hummer, and A. Szabo, "Intrinsic rates and activation free energies from single-molecule pulling experiments", *Physical Review Letters* **96**, 108101 (2006). https://doi.org/10.1103/PhysRevLett.96.108101

Please cite the references above if you use `ForceSpectroscopyMLE` to analyze your data.  



## Installation

Currently, the package is not in a registry.  It must therefore be added by specifying a URL to the repository:
```julia
using Pkg; Pkg.add(url="https://github.com/bio-phys/ForceSpectroscopyMLE")
```
Users of older software versions may need to wrap the contents of the brackets with `PackageSpec()`.  



## Usage

### Importing data

The rupture force data should be of the type `Array{Float64,2}`, where the first column contains the rupture forces `F` (in pN) and the second column the associated loading rates `dF` (in pN/s).  In principle, users can lump all their measured force spectra into a single file and, e.g., read it in as follows:
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
[βΔG_u, x_u, k_0] = MLE_estimator(data,ν)
```
The parameter `ν` can be set to `1/2` or `2/3` depending on the shape of the underlying free-energy landscape.  For `ν = 1` the DHS model reduces to the Bell-Evans model, which only depends on the parameters `x_u` and `k_0`.  