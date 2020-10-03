# Imoto_Cancers_2020

This repository contains modeling code for the following paper:<br>
- Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data — Application to the ErbB Receptor Signaling Pathway. *Cancers* (2020)

The paper can be accessed at the [Cancers website]().

## Description

- training/
    - Training model parameters using [ParamEstim](https://github.com/himoto/ParamEstim) (v0.2.0).
    - Requirements
        > - [Sundials](https://github.com/SciML/Sundials.jl)
        > - [SteadyStateDiffEq](https://github.com/SciML/SteadyStateDiffEq.jl)
        > - [StatsBase](https://github.com/JuliaStats/StatsBase.jl)
        > - [CSV](https://github.com/JuliaData/CSV.jl)
        > - [DataFrames](https://github.com/JuliaData/DataFrames.jl)
        > - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
        > - [Seaborn](https://github.com/JuliaPy/Seaborn.jl)

- python/
    - Validation and prediction using [BioMASS](https://github.com/okadalabipr/biomass) (v0.1.0).
    - Requirements
        > - numpy
        > - scipy
        > - matplotlib
        > - seaborn

- gene_expression/
    - CCLE RNA-seq gene expression data used for the model individualization

## Usage
1. Parameter estimation
    ```bash
    $ cd trainig
    $ sh optimize_parallel.sh
    ```

1. Visualization of simulation results
    ```bash
    $ cd python
    ```
    ```python
    import SKBR3
    from biomass import run_simulation
    run_simulation(SKBR3, viz_type='average', show_all=False, stdev=True)
    ```

1. Sensitivity analysis
    ```python
    from biomass import run_analysis
    run_analysis(SKBR3, target='initial_condition', metric='integral', style='heatmap')
    ```

## License
[MIT](LICENSE)