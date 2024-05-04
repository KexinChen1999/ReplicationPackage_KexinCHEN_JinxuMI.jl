```@meta
CurrentModule = ReplicationPackage_KexinCHEN_JinxuMI
```

# ReplicationPackage_KexinCHEN_JinxuMI

Documentation for [ReplicationPackage_KexinCHEN_JinxuMI](https://github.com/KexinChen1999/ReplicationPackage_KexinCHEN_JinxuMI.jl).

<p align="justify"> We replicate the simulations of the three models of land distribution after reform in the original paper. The original models are calibrated and simulated in Matlab, while we used Julia (http://julialang.org/), an open source, free of charge computational lanaguage, to replicate the paper. Our work provides a new way to use the model in the original paper without buying a liscense from MathWorks. The package sheds new light on using Julia to estimate the effects of land reforms. </p>




# Summary of the Paper

<p align="justify"> The paper's primary purpose is to assess the effects of a significant land reform on agricultural productivity in the Philippines. The 1988 Comprehensive Agrarian Reform Program put a ceiling on the quantity of private land holdings, leading to a drop in farm size and productivity. Using the Decennial Agricultural Census Data and Philippines Cash Cropping Project survey data, the authors can construct balanced panel data that tracks the farm-level changes before and after the reform, controlling for farmer ability and location. </p>

<p align="justify"> The Philippines' 1988 Comprehensive Agrarian Reform Law put a 5-hectare ceiling on individual land ownership. The redistribution at fair market value was implemented by the end of 1990. The data show that land redistribution led to a -11.6% decrease in land productivity from 1990 to 1993. Even though the ceiling was not perfectly implemented, this study still discovered a significant change in farm size. The share of small farms almost doubled while the share of more giant farms shrank. This phenomenon was verified by farm-level survey data as well.  </p>

<p align="justify"> The authors also build a model to estimate the effect of reform distortion quantitatively. In the production sector, a farm produces with farmer productivity, economy-wide productivity, land input, and labor input. The farm can choose from producing cash crops or food crops. Therefore, production would be influenced by crop-specific productivity parameters and fixed costs. After profit maximization, the authors have the input demand functions and profits linear in productivity and market distortion parameters. Farmers will choose to become workers, food crop farm operators, or cash crop farm operators based on their individual productivity and occupation-choice thresholds.</p>

<p align="justify"> Then, the authors simulated the land holding reform, embodied in a land ceiling, degree of enforcement, and reform beneficiaries who receive an amount of land. Assuming every recipient receives the same amount of land, an equilibrium after the government-mandated land redistribution characterized could be described by two equations. The authors also estimate a market-based redistribution where the excessive land above the ceiling is distributed via a rental market. The results show that the market-based counterfactual land reform influenced both the land demand for hired labor and the occupational choice of farmers. </p>

<p align="justify"> After parameterization, the authors calibrate the model and match the parameters with data. The joint population distribution of farm-level distortion and farmer ability is sampled from a bivariate log-normal distribution. Other parameters are obtained by jointly matching population moments with observed moments. After calibration, the model matches quite well with the observed data, especially among the smaller farms, which are the primary targets of the reform. Therefore, the authors can estimate the reform's effect in the model after choosing the four targeted parameters of land ceiling, enforcement intensity, and fraction of landless or small-holders. </p>

<p align="justify"> This work clearly shows that government-mandated land reform leads to a drop in farm size, agricultural labor productivity, and share of landless farmers, and the drop is increasing as enforcement gets more strict. The mechanism of government-mandated reform is driven by the redistribution based on land holding, not productivity. The increase in landed farmers working on food crop farms reduces productivity, labor shortage, and higher wages, which will spill over to the cash crop sector, leading to lower-ability cash croppers switching to food crops. The market-based reform mechanism reduces the holding of all constrained farmers, reducing land rents and wages. Every unconstrained farmer can raise their demand for labor and land. Low factor prices also increase cash crops' relative profits, encouraging more food crop farmers to switch to cash crops. Meanwhile, long-term productivity growth can be explained by other exogenous factors. The adverse effects are also shown by reduced form regression analysis. </p>

<p align="justify"> The significant contribution of the paper is that it focuses on a specific farm-size land reform distortion using micro-level data; the occupational and technological choices between food and cash crops caused by the distortion were another primary focus of the paper. It sheds new light on the existing literature about land reform and agricultural productivity, especially the role of ownership reform. </p>




# Installation

using latest master directly from github:

```julia
Pkg.add(url="https://github.com/KexinChen1999/ReplicationPackage_KexinCHEN_JinxuMI.jl")
```

if you are just interested in using the Calibration of Benchmark Economy (BE) part:
```julia
Pkg.add(url="https://github.com/KexinChen1999/BE.jl")
```

if you are just interested in using Goverment-mandated Land Reform (LR_main) part:
```julia
Pkg.add(url="https://github.com/KexinChen1999/LR_main.jl")
```

if you are just interested in using the Market-based Land Reform (LR_market) part:
```julia
Pkg.add(url="https://github.com/KexinChen1999/LR_market.jl")
```




# Usage

To show how the `ReplicationPackage_KexinCHEN_JinxuMI` package can be used, several Julia packages are required to be installed and used, including Pkg, Printf, MAT, Distributions, LinearAlgebra, Statistics, DelimitedFiles, Optim, NLsolve, Plots, Serialization, CSV, and DataFrames. Also, we assume you have already installed `ReplicationPackage_KexinCHEN_JinxuMI` as described above.

First, we focus on the the Calibration of Benchmark Economy (BE) part and define the `BE_eval` function:

```julia
function BE_eval(x, A)
    global KAPPAc, KAPPAf, LN, g_vec, s_vec, phi_vec, Pc, Pf, GAMMA, ALPHA, N, Nw, hired_lab_sh, cash_oper_sh, cdf_g

    Cf = x[1]
    Cc = x[2]

    # Compute occupational choice cutoffs
    INl = findfirst(cdf_g .> hired_lab_sh)
    g_lbar = g_vec[INl]
    sh_l = cdf_g[INl]
    g_lbar_Indic = findlast(g_vec .< g_lbar)

    INu = findfirst(cdf_g .> (hired_lab_sh .+ (1 .- hired_lab_sh) .* (1 .- cash_oper_sh)))
    g_ubar = g_vec[INu]
    sh_u = cdf_g[INu]
    g_ubar_Indic = findfirst(g_vec .> g_ubar)

    # Factor prices - from model equations
    w = (Cc .- Cf) ./ ((g_ubar ./ g_lbar) .* ((Pc ./ Pf).^(1 ./ (1 .- GAMMA)) .* (KAPPAc ./ KAPPAf) .- 1)) .- Cf
    q = ALPHA .* ((g_lbar .* (GAMMA.^(GAMMA ./ (1 .- GAMMA))) .* (1 .- GAMMA) .* (((1 .- ALPHA) ./ w).^(GAMMA .* (1 .- ALPHA) ./ (1 .- GAMMA))) .* (Pf).^(1 ./ (1 .- GAMMA)) .* (A .* KAPPAf)) ./ (w .+ Cf)).^((1 .- GAMMA) ./ (ALPHA .* GAMMA))
    qw_ratio = q ./ w

    # Compute occupational choice vectors based on cutoffs
    of_vec = zeros(Int, N)
    of_vec[g_lbar_Indic+1:g_ubar_Indic] .= 1

    oc_vec = zeros(Int, N)
    oc_vec[g_ubar_Indic+1:end] .= 1

    # Solve problem under each technology for every individual
    lf_vec = ((ALPHA ./ q).^((1 .- (1 .- ALPHA) .* GAMMA) ./ (1 .- GAMMA))) .* (((1 .- ALPHA) ./ w).^(GAMMA .* (1 .- ALPHA) ./ (1 .- GAMMA))) .* ((GAMMA .* Pf).^(1 ./ (1 .- GAMMA))) .* (A .* KAPPAf) .* g_vec
    lc_vec = ((ALPHA ./ q).^((1 .- (1 .- ALPHA) .* GAMMA) ./ (1 .- GAMMA))) .* (((1 .- ALPHA) ./ w).^(GAMMA .* (1 .- ALPHA) ./ (1 .- GAMMA))) .* ((GAMMA .* Pc).^(1 ./ (1 .- GAMMA))) .* (A .* KAPPAc) .* g_vec
    nl_ratio = ((1 - ALPHA) ./ ALPHA) .* qw_ratio
    nf_vec = nl_ratio .* lf_vec
    nc_vec = nl_ratio .* lc_vec
    yf_vec = (A * KAPPAf .* s_vec).^(1 - GAMMA) .* (lf_vec.^ALPHA .* nf_vec.^(1 .- ALPHA)).^GAMMA
    yc_vec = (A * KAPPAc .* s_vec).^(1 - GAMMA) .* (lc_vec.^ALPHA .* nc_vec.^(1 .- ALPHA)).^GAMMA

    # Compute aggregates using occupational choice vectors
    Nw = sum(1 .- oc_vec .- of_vec) ./ N  # share of hired labor
    Nf = sum(of_vec) ./ N  # share of food farm operators
    Nc = sum(oc_vec) ./ N  # share of cash farm operators
    LAB_f = (sum(of_vec .* nf_vec) ./ N) .+ Nf  # total share of labor in food farms (hired+operators)
    LAB_c = (sum(oc_vec .* nc_vec) ./ N) .+ Nc  # total share of labor in cash farms (hired+operators)

    # Clear land and labor markets
    f1 = (sum(of_vec .* lf_vec) ./ N) + (sum(oc_vec .* lc_vec) ./ N) .- LN  # land market clearing condition
    f2 = (sum(of_vec .* nf_vec) ./ N) + (sum(oc_vec .* nc_vec) ./ N) .- Nw  # labor market clearing condition

    return [f1, f2]
end
```








# Providing our results




# Comparison with the original results




# Conclusion




# Acknowledgements
<p align="justify"> We extend our heartfelt appreciation to Professor Florian Oswald for his exceptional mentorship and dedication to our academic growth. His teaching in Introduction to Programming last year and Computational Economics this year has been nothing short of transformative. Professor Oswald's clarity of instruction, insightful guidance, and unwavering support have enabled us to grasp complex concepts and excel in our studies. We are immensely grateful for the knowledge and skills we have gained under his tutelage. </p>




