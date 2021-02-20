# Tau_Transport
Lead Author: [Justin Torok](http://github.com/justin-torok)

Email: jut2008@med.cornell.edu

The following code was developed for running the tau axonal transport model described in Torok, *et al.*, 2021 (***insert bioRxiv link here***), along with all of the auxilliary functions required for plotting the outputs as shown in the manuscript.

## 1. Setup
All code is written in MATLAB and has been tested in versions 2019b and 2020a. However, the developers do not anticipate difficulties with using previous MATLAB versions for any of the functions contained within this package. There are no toolboxes required to run the model, but the Statistics & Machine Learning Toolbox is required for the `corr.m` function.

In order to reproduce the figures within the manuscript, we have provided the model outputs in a separate folder named [SampleFiles](https://drive.google.com/file/d/14KIloSmoMDG1tMROd_orl0s2PzqMbpTT/view?usp=sharing), which contains all of the simulation results required as well as the network bias data from [Mezias & Raj, 2020](https://www.biorxiv.org/content/10.1101/2020.11.06.371625v1).

## 2. Files
Below is a short description of each of the code files contained in the **Tau_Transport** repository, grouped by general functionality in alphabetical order. Scripts that are also functions have their inputs and outputs described, with required inputs in boldface text and optional inputs with their default setting in parentheses.

### Running the Model
- `TauTransportPDE.m`: The core function used to generate all model outputs, which solves a system of coupled 1D PDEs to determine the spatiotemporal profiles of soluble and insoluble pathological tau (refer to Torok, *et al.*, 2021(***insert bioRxiv here***) and the documentation inside of the function itself for full model details). All inputs are optional and specified as keyword arguments using `inputParser`. This function uses the MATLAB builtin `pdepe`, whose documentation can be found [here](https://www.mathworks.com/help/matlab/ref/pdepe.html); we recommend users familiarize themselves with `pdepe` before attempting to modify this function. The wrapper scripts described below are designed to reproduce the output files required for the various analyses explored in the manuscript, and we suggest that users use these as templates for running the model as opposed to directly calling this function. 
    - ***Inputs***:
        - *`alpha`: Linear growth rate term (default = 0)*
        - `beta`: Unimolecular fragmentation rate (default = 1e-06)
        - `gamma`: Bimolecular aggregation rate (default = 2e-05)
        - `delta`: Anterograde velocity enhancement factor (default = 1)
        - `epsilon`: Anterograde velocity reduction factor (default = 0.01)
        - `lambda`: Diffusivity barrier strength between compartments (default = 0.01)
        - `frac`: Fraction of tau undergoing diffusion (default = 0.92)
        - `T`: End simulation time in seconds (default = 5e07)
        - `tsteps`: Number of time steps at which `pdepe` evaluates the model solutions (default = 1000)
        - `L1`: Length of the presynaptic somatodendritic compartment in microns (default = 200)
        - `L2`: Length of the postsynaptic somatodendritic compartment in microns (default = 200)
        - `L_int`: Total length of the internal compartments (axon + AIS + synaptic cleft) in microns (default = 1000)
        - `L_ais`: Length of the axon initial segment in microns (default = 40)
        - `L_syn`: Length of the synaptic cleft in microns (default = 40)
        - `n0`: Initial concentration profile of soluble tau within the internal compartments, specified at each micron (default = zeros(1,L_int))
        - `m0`: Initial concentration profile of insoluble tau within the internal compartments, specified at each micron (default = zeros(1,L_int))
        - `N1_0`: Initial presynaptic somatodendritic soluble tau concentration (default = 0)
        - `N2_0`: Initial postsynaptic somatodendritic soluble tau concentration (default = 0)
        - `M1_0`: Initial presynaptic somatodendritic insoluble tau concentration (default = 0)
         - `M2_0`: Initial postsynaptic somatodendritic insoluble tau concentration (default = 0)    
         - `resmesh`: Coarseness of the inhomogeneous x-mesh used; can take values of 'coarse' or 'fine'. Use 'coarse' for faster, less precise simulations (default = 'fine')   
    - ***Outputs***:
        - `n`: A `tsteps` x `length(xmesh)` numeric array of soluble tau concentrations at each spatial coordinate in xmesh for all time points in `trange`.
        - `m`: A `tsteps` x `length(xmesh)` numeric array of insoluble tau concentrations at each spatial coordinate in xmesh for all time points in `trange`.
        - `xmesh`: A vector of spatial coordinates where the model solutions `n` and `m` are evaluated. Length depends on the `resmesh` input and spacing is inhomogeneous.
        - `trange`: A 1 x `tsteps` vector of time points at which the model solutions `n` and `m` are evaluated. Time scale is logarithmic after the first 1000s of model time, with final time point `T`.
        - `jn`: A `tsteps` x `length(xmesh)` numeric array of the estimated flux of `n` at each spatial coordinate in xmesh for all time points in `trange`.
        - `jm`: A `tsteps` x `length(xmesh)` numeric array of the estimated flux of `m` at each spatial coordinate in xmesh for all time points in `trange`.
- `TauTransportPDE_Wrapper_EqHeatmap.m`: Wrapper script for generating output files from `TauTransportPDE.m` for generating equilibrium heatmaps with respect to `delta` and `epsilon`.
- `TauTransportPDE_Wrapper_RandomIC.m`: Wrapper script for generating output files from `TauTransportPDE.m` using random initial conditions.
- `TauTransportPDE_Wrapper_Single.m`: Wrapper script for generating output files from `TauTransportPDE.m` for single parameter sets.

### Plotting the Outputs
While several of these functions are hard-coded to load specific model output files from the "SampleFiles" folder (see above), they can be easily modified to take in different model output files that a user generates.
- `BiasComparisonPlotter.m`: Produces plots of model bias vs. observed network bias from [Mezias & Raj, 2020](https://www.biorxiv.org/content/10.1101/2020.11.06.371625v1). Model outputs used by this function were generated using the `TauTransportPDE_Wrapper_Single.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - `savenclose`: Binary flag that, if set to 1, will save all figures as 300 DPI .png files in a folder called "OutputFigures" within the `basepath` directory (default = 0)
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)
- `ConvergencePlotter.m`: Produces plots of model convergence with respect to initial conditions for three combinations of `delta` and `epsilon` values. Model outputs used by this function were generated using the `TauTransportPDE_Wrapper_RandomIC.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - `savenclose`: Binary flag that, if set to 1, will save all figures as 300 DPI .png files in a folder called "OutputFigures" within the `basepath` directory (default = 0)
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)
- `EqHeatmapPlotter.m`: Produces plots of equilibrium heatmaps with respect to `delta`, `epsilon`, and three values of a specified parameter. Model outputs used by this function were generated using the `TauTransportPDE_Wrapper_EqHeatmap.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - **`paramname`**: Character array specifying the parameter of interest, which can take values of 'beta', 'gamma', and 'frac'. 
        - `savenclose`: Binary flag that, if set to 1, will save all figures as 300 DPI .png files in a folder called "OutputFigures" within the `basepath` directory (default = 0)
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)
- `redblue.m` (**Not generated by author of this repository**): Developed by Adam Auton in 2009, this function was obtained off of the [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap) and generates a red-blue colormap, which is used in `EqHeatmapPlotter.m`. See source documentation for more details.
- `SnapshotPlotter.m`: Produces plots of soluble and insoluble tau concentration profiles at specified time points. Model outputs used by this function can be generated using the `TauTransportPDE_Wrapper_Single.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - **`matstr`**: Character array specifying the name of the output file (without the '.mat' appendix) within the "SampleFiles" folder to load
        - **`customt`**: Vector of valid indices within `trange` at which to display outputs. *Not recommended for more than five time points*
        - `savenclose`: Binary flag that, if set to 1, will save all figures as 300 DPI .png files in a folder called "OutputFigures" within the `basepath` directory (default = 0)
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)
- `SummaryPlotter.m`: Produces plots of somatodendritic tau deposition over time for three combinations of `delta` and `epsilon` values. Model outputs used by this function were generated using the `TauTransportPDE_Wrapper_Single.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - `savenclose`: Binary flag that, if set to 1, will save all figures as 300 DPI .png files in a folder called "OutputFigures" within the `basepath` directory (default = 0)
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)
- `VideoMaker.m`: Produces videos of model outputs over time. Model outputs used by this function can be generated using the `TauTransportPDE_Wrapper_Single.m` wrapper script structure, with different parameter inputs.
    - ***Inputs***:
        - **`matstr`**: Character array specifying the name of the output file (without the '.mat' appendix) within the "SampleFiles" folder to load
        - `basepath`: Character array specifying the enclosing folder of the "SampleFiles" directory, where the required dependencies are stored (default = cd)
    - ***Outputs*** (none)