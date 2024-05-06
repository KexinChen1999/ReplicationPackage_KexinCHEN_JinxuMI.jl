```@meta
CurrentModule = ReplicationPackage_KexinCHEN_JinxuMI
```

# ReplicationPackage_KexinCHEN_JinxuMI

Documentation for [ReplicationPackage_KexinCHEN_JinxuMI](https://github.com/KexinChen1999/ReplicationPackage_KexinCHEN_JinxuMI.jl).

We replicate the simulations of the three models of land distribution after reform in the original paper[(Adamopoulos & Restuccia, 2020)](https://www.aeaweb.org/articles?id=10.1257/mac.20150222). 
<p align="justify"> The original models are calibrated and simulated in Matlab, while we used Julia (http://julialang.org/), an open source, free of charge computational lanaguage, to replicate the paper. Our work provides a new way to use the model in the original paper without buying a liscense from MathWorks. Furthermore, our results can also be used to verify the original article. </p>




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




# Datasets

You can rebuild all datasets except the raw data file with the  `BE.jl` or the first part of the `ReplicationPackage_KexinCHEN_JinxuMI`:

Dataset Name |  Description
-------------| -------------
AbDist_matrix.mat	           | 	 Raw Data File
BE_var_julia.mat	       | 	 Output - The first set of variables
BE_parameters_julia.mat	               | 	Output - The second set of variables
BE_values_julia.mat	       | 	 Output - The third set of variables





# Usage

To show how the `ReplicationPackage_KexinCHEN_JinxuMI` package can be used, several Julia packages are required to be installed and used, including `Pkg`, `Printf`, `MAT`, `Distributions`, `LinearAlgebra`, `Statistics`, `DelimitedFiles`, `Optim`, `NLsolve`, `Plots`, `Serialization`, `CSV`, and `DataFrames`. Also, we assume you have already installed `ReplicationPackage_KexinCHEN_JinxuMI` as described above.

## Attention⚠️
You must change the current working directory using the cd() function where you put the data.
```julia
cd("path/to/directory")
```

## Calibration of Benchmark Economy
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

The `BE_eval function` defines the market-clearing conditions of our optimization problem. It acts as the constraint in our final step of optimization.

<p align="justify"> The function first determines cutoff points (g_lbar and g_ubar) for choosing between food farming, cash crop farming, and possibly other occupations. These are calculated based on cumulative distribution functions and certain market share thresholds. Then it performs factor price calculations. It computes the wages (w) and another price factor (q), based on the model equations which involve ratios of prices, costs, and technological parameters. Then, it derives the occupational Vectors. It creates vectors that define whether individuals choose food farming, cash crop farming, or neither based on the earlier calculated g-cutoffs. After that, it calculates individual productivity and income by deriving labor and capital inputs (lf_vec, lc_vec, nf_vec, nc_vec), outputs (yf_vec, yc_vec), and profits (PIf_vec, PIc_vec) for each individual, adjusted by their ability, soil quality, and chosen occupation. Finally, it derives the aggregates and market clearing conditions. It computes total outputs, labor shares, and productivity ratios for the entire economy. It ensures that the land and labor markets are clear, i.e., total demand for land and labor matches supply. The returned values are a vector f containing the residuals of the market clearing conditions for land and labor. These residuals are used to assess how well the model's assumptions and parameters fit the actual economic conditions. </p>

Initial Guess for the variables:
```julia
guess = [-2.4, -2.0]
```

Parameters:
```julia
GAMMA       = 0.7;
ALPHA       = 0.3;
A           = 1; 
KAPPAf      = 1; 
Pf          = 1; 
Pc          = 1; 
LN          = AFS*(1-hired_lab_sh); 
KAPPAc      = 1.25
```

Additional parameters:
```julia
params = [A] 
```

<p align="justify"> Then we use nlsolve function to find the numerical solutions of our non-linear system of equations. Given the market-clearing conditions calculated from BE_eval() and specified set of varaibles we are interested in,  the algorithms uses numerical methods to approximate the roots of multiple equations simultaneously. We can derive the optimal distribution of lands and labor from the optimization results. </p>

```julia
result = nlsolve((res, x) -> res .= BE_eval(x, params), guess, show_trace=true, xtol=1e-16)
```

Extracting the solution:
```julia
x = result.zero
Cf = x[1]
Cc = x[2]
```

Check for convergence:
```julia
converged = result.f_converged
println("Converged: ", converged)
```

Compute occupational choice cutoffs:
```julia
#Choose g_lbar to match a share of hired labor in total labor
INl = findfirst(cdf_g .> hired_lab_sh)
g_lbar = g_vec[INl]
sh_l = cdf_g[INl]
g_lbar_Indic = findfirst(x -> x > g_lbar, g_vec) - 1

#Choose g_ubar to match a share of food crop operators in total operators
INu = findfirst(cdf_g .> (hired_lab_sh + (1 - hired_lab_sh) * (1 - cash_oper_sh)))
g_ubar = g_vec[INu]
sh_u = cdf_g[INu]
g_ubar_Indic = findfirst(x -> x > g_ubar, g_vec) - 1
```

Factor prices - from model equations:
```julia
w = (Cc - Cf) / ((g_ubar / g_lbar) * ((Pc / Pf)^(1 / (1 - GAMMA)) * (KAPPAc / KAPPAf) - 1)) - Cf
q = ALPHA * (((g_lbar * (GAMMA^(GAMMA / (1 - GAMMA))) * (1 - GAMMA) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * (Pf)^(1 / (1 - GAMMA)) * (A * KAPPAf)) / (w + Cf))^((1 - GAMMA) / (ALPHA * GAMMA)))
qw_ratio = q / w
```

Initialize of_vec with ones and then set specific elements to zero based on conditions:
```julia
of_vec = ones(Int, N)
of_vec[1:g_lbar_Indic-1] .= 0
of_vec[g_ubar_Indic+1:N] .= 0
```

Initialize oc_vec with ones and then set specific elements to zero based on conditions:
```julia
oc_vec = ones(Int, N)
oc_vec[1:g_ubar_Indic] .= 0
```

Solve problem under each technology for every individual:
```julia
lf_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pf)^(1 / (1 - GAMMA))) * (A * KAPPAf) .* g_vec
lc_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pc)^(1 / (1 - GAMMA))) * (A * KAPPAc) .* g_vec
nl_ratio = ((1 - ALPHA) / ALPHA) * qw_ratio
nf_vec = nl_ratio .* lf_vec
nc_vec = nl_ratio .* lc_vec
yf_vec = (A * KAPPAf .* s_vec) .^ (1 - GAMMA) .* (lf_vec .^ ALPHA .* nf_vec .^ (1 - ALPHA)) .^ GAMMA
yc_vec = (A * KAPPAc .* s_vec) .^ (1 - GAMMA) .* (lc_vec .^ ALPHA .* nc_vec .^ (1 - ALPHA)) .^ GAMMA
PIf_vec = (1 - GAMMA) * Pf * yf_vec .* (phi_vec .^ (1 - GAMMA)) - Cf .* ones(N)
PIc_vec = (1 - GAMMA) * Pc * yc_vec .* (phi_vec .^ (1 - GAMMA)) - Cc .* ones(N)
```


Compute aggregates using occupational choice vectors:
```julia
# Labor Shares - based on occupational choices
Nw = sum(1 .- oc_vec .- of_vec) / N  # Share of hired labor
Nf = sum(of_vec) / N  # Share of food farm operators
Nc = sum(oc_vec) / N  # Share of cash farm operators
LAB_f = (sum(of_vec .* nf_vec) / N) + Nf  # Total share of labor in food farms (hired + operators)
LAB_c = (sum(oc_vec .* nc_vec) / N) + Nc  # Total share of labor in cash farms (hired + operators)

# Crop Outputs - based on occupational choices
Yc = sum(oc_vec .* yc_vec) / N  # Total output of cash crop farms
Yf = sum(of_vec .* yf_vec) / N  # Total output of food crop farms
```

Ratio of cash-to-food labor productivities (INCLUDES OPERATORS):
```julia
YNR_cf_model = (Yc / LAB_c) / (Yf / LAB_f)
```

Clear land and labor markets:
```julia
f1 = (sum(of_vec .* lf_vec) / N) + (sum(oc_vec .* lc_vec) / N) - LN  # Land market clearing condition
f2 = (sum(of_vec .* nf_vec) / N) + (sum(oc_vec .* nc_vec) / N) - Nw  # Labor market clearing condition
f = [f1, f2]  # Construct the vector f from f1 and f2
```

Production choices for all individuals (including non-operators)
```julia
# Farm size vector (land input)
l_vec = of_vec .* lf_vec + oc_vec .* lc_vec

# Hired labor vector
n_vec = of_vec .* nf_vec + oc_vec .* nc_vec

# Output vector
y_vec = of_vec .* yf_vec + oc_vec .* yc_vec
```


Calculate distributions and statistics of interest
```julia
# Truncated distribution - conditional on operating (to calculate averages across active units)
Indic = oc_vec .+ of_vec
COND_distr = Indic ./ sum(Indic)

# Average Farm Size (AFS) - across operators
l_value = sum(COND_distr .* lf_vec .* of_vec) + sum(COND_distr .* lc_vec .* oc_vec)
AFS_BE = l_value

# Landless share
landless_BE = Nw

# parameters used to generate the size distribution (bins)
bin1 = 1;    
bin2 = 2;
bin3 = 5; 
bin4 = 7;
bin5 = 10;
bin6 = 15;


# Size distribution (culmulative percentage of land/farms) - find farms (operators) in each of the above bins
lhat1in = findfirst(x -> x > bin1, l_vec)
lhat2in = findfirst(x -> x > bin2, l_vec)
lhat3in = findfirst(x -> x > bin3, l_vec)
lhat4in = findfirst(x -> x > bin4, l_vec)
lhat5in = findfirst(x -> x > bin5, l_vec)
lhat6in = findfirst(x -> x > bin6, l_vec)

# Check if the last bin is empty
if isnothing(lhat6in)
    
    # Farm-Size Distribution over rest of bins
    farm1 = sum(COND_distr[1:lhat1in-1])
    farm2 = sum(COND_distr[lhat1in:lhat2in-1])
    farm3 = sum(COND_distr[lhat2in:lhat3in-1])
    farm4 = sum(COND_distr[lhat3in:lhat4in-1])
    farm5 = sum(COND_distr[lhat4in:lhat5in])
    farm6 = sum(COND_distr[lhat5in:N])
    farm7 = 0
    
    # Distribution of Land over rest of bins
    land1 = sum(l_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1]) / l_value
    land2 = sum(l_vec[lhat1in:lhat2in-1] .* COND_distr[lhat1in:lhat2in-1]) / l_value
    land3 = sum(l_vec[lhat2in:lhat3in-1] .* COND_distr[lhat2in:lhat3in-1]) / l_value
    land4 = sum(l_vec[lhat3in:lhat4in-1] .* COND_distr[lhat3in:lhat4in-1]) / l_value
    land5 = sum(l_vec[lhat4in:lhat5in] .* COND_distr[lhat4in:lhat5in]) / l_value
    land6 = sum(l_vec[lhat5in:N] .* COND_distr[lhat5in:N]) / l_value
    land7 = 0
    
    # Distribution of Output over rest of bins
    output1 = sum(y_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1])
    # Similar for output2 to output6, as with MATLAB code
    output7 = 0
    
    # Distribution of Hired Labor over rest of bins
    hirelab1 = sum(n_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1])
    # Similar for hirelab2 to hirelab6, as with MATLAB code
    hirelab7 = 0

else

    # Farm-Size distribution of Farms over all bins
    farm1 = sum(COND_distr[1:lhat1in-1])
    farm2 = sum(COND_distr[lhat1in:lhat2in-1])
    farm3 = sum(COND_distr[lhat2in:lhat3in-1])
    farm4 = sum(COND_distr[lhat3in:lhat4in-1])
    farm5 = sum(COND_distr[lhat4in:lhat5in-1])
    farm6 = sum(COND_distr[lhat5in:lhat6in])
    farm7 = sum(COND_distr[lhat6in:end])

    # Distribution of Land over all bins
    land1 = sum(l_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1]) / l_value
    land2 = sum(l_vec[lhat1in:lhat2in-1] .* COND_distr[lhat1in:lhat2in-1]) / l_value
    land3 = sum(l_vec[lhat2in:lhat3in-1] .* COND_distr[lhat2in:lhat3in-1]) / l_value
    land4 = sum(l_vec[lhat3in:lhat4in-1] .* COND_distr[lhat3in:lhat4in-1]) / l_value
    land5 = sum(l_vec[lhat4in:lhat5in-1] .* COND_distr[lhat4in:lhat5in-1]) / l_value
    land6 = sum(l_vec[lhat5in:lhat6in] .* COND_distr[lhat5in:lhat6in]) / l_value
    land7 = sum(l_vec[lhat6in:end] .* COND_distr[lhat6in:end]) / l_value

    # Distribution of Output over all bins
    output1 = sum(y_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1])
    output2 = sum(y_vec[lhat1in:lhat2in-1] .* COND_distr[lhat1in:lhat2in-1])
    output3 = sum(y_vec[lhat2in:lhat3in-1] .* COND_distr[lhat2in:lhat3in-1])
    output4 = sum(y_vec[lhat3in:lhat4in-1] .* COND_distr[lhat3in:lhat4in-1])
    output5 = sum(y_vec[lhat4in:lhat5in-1] .* COND_distr[lhat4in:lhat5in-1])
    output6 = sum(y_vec[lhat5in:lhat6in] .* COND_distr[lhat5in:lhat6in])
    output7 = sum(y_vec[lhat6in:end] .* COND_distr[lhat6in:end])

    # Distribution of Hired Labor over all bins
    hirelab1 = sum(n_vec[1:lhat1in-1] .* COND_distr[1:lhat1in-1])
    hirelab2 = sum(n_vec[lhat1in:lhat2in-1] .* COND_distr[lhat1in:lhat2in-1])
    hirelab3 = sum(n_vec[lhat2in:lhat3in-1] .* COND_distr[lhat2in:lhat3in-1])
    hirelab4 = sum(n_vec[lhat3in:lhat4in-1] .* COND_distr[lhat3in:lhat4in-1])
    hirelab5 = sum(n_vec[lhat4in:lhat5in-1] .* COND_distr[lhat4in:lhat5in-1])
    hirelab6 = sum(n_vec[lhat5in:lhat6in] .* COND_distr[lhat5in:lhat6in])
    hirelab7 = sum(n_vec[lhat6in:end] .* COND_distr[lhat6in:end])

end
```


Conditional (on operating) distributions over specified bins
```julia
# Conditional (on operating) distribution of LAND INPUT across specified bins - note this is the same as land_pdf_model
l_model = [land1, land2, land3, land4, land5, land6, land7]

# Conditional (on operating) distribution of FARMS (also OPERATORS) across specified bins - note this is the same as farm_pdf_model
f_model = [farm1, farm2, farm3, farm4, farm5, farm6, farm7]

# Conditional (on operating) distribution of OUTPUT across specified bins
y_model = [output1, output2, output3, output4, output5, output6, output7]

# Conditional (on operating) distribution of HIRED LABOR across specified bins
h_model = [hirelab1, hirelab2, hirelab3, hirelab4, hirelab5, hirelab6, hirelab7]

# Conditional (on operating) distribution of HIRED LABOR+OPERATORS across specified bins
n_model = h_model .+ f_model
```

Find active farm operators
```julia
indACTIVE = findall(x -> x == 1, Indic)
```

Distribution of TFP
```julia
# TFP of Active Units
TFPf_vec = (A*KAPPAf.*s_vec).^(1-GAMMA)
TFPc_vec = (A*KAPPAc.*s_vec).^(1-GAMMA)

# Initialize of_vec with ones and then set specific elements to zero based on conditions
#of_vec = ones(Int, N)
#of_vec[1:g_lbar_Indic-1] .= 0
#of_vec[g_ubar_Indic+1:N] .= 0

# Initialize oc_vec with ones and then set specific elements to zero based on conditions
#oc_vec = ones(Int, N)
#oc_vec[1:g_ubar_Indic] .= 0

TFP_vec = of_vec .* TFPf_vec .+ oc_vec .* TFPc_vec
TFP_vec_ac = TFP_vec[indACTIVE]

# Log-TFP of Active Units
lTFP_vec_ac = log.(TFP_vec_ac)  # Use broadcasting with 'log.'

# Compute the STD of log TFP (including the KAPPAs)
STDlTFPac = std(lTFP_vec_ac)

# Distribution of TFPR
# TFPRs of all individuals
TFPR_vec = phi_vec .^ (-(1 - GAMMA))

# TFPRs of active farmers
TFPR_vec_ac = TFPR_vec[indACTIVE]

# log of active TFPRs
lTFPR_vec_ac = log.(TFPR_vec_ac)  # Use broadcasting with 'log.'

# Standard Deviation of active TFPRs
println("\nSTD OF log-TFPR")
STD_lTFPRac = std(lTFPR_vec_ac)
VAR_lTFPR_model = STD_lTFPRac^2


#Correlation of log(TFP)-log(TFPR) across active farms
println("\nCORR of log-TFP AND log-TFPR")
CORR_lTFP_lTFPR_model = cor(lTFP_vec_ac, lTFPR_vec_ac)
println(CORR_lTFP_lTFPR_model)
```

OUTPUT / OUTPUT PER WORKER 🌟
```julia
# Aggregate Output per capita (since persons = workers, and only one sector this is also aggregate output per worker - and includes operators + hired workers)
VApw_BE = (sum(oc_vec .* yc_vec) + sum(of_vec .* yf_vec)) / N

# Aggregate Output Per Hired Worker
VAphw = VApw_BE / Nw

# Output per worker (hired + operators) in cash crops
VApw_cash = Yc / LAB_c

# Output per worker (hired + operators) in food crops
VApw_food = Yf / LAB_f




# Save the first set of variables to BE_var.mat
vars_BE_var = Dict(
    "g_vec" => g_vec,
    "s_vec" => s_vec,
    "phi_vec" => phi_vec,
    "oc_vec" => oc_vec,
    "of_vec" => of_vec,
    "COND_distr" => COND_distr,
    "lf_vec" => lf_vec,
    "lc_vec" => lc_vec,
    "l_value" => l_value,
    "l_vec" => l_vec,
    "Indic" => Indic,
    "cdf_g" => cdf_g
)
matwrite("BE_var_julia.mat", vars_BE_var)

# Save the second set of variables to BE_parameters.mat
# Save the first set of variables to BE_var.mat
vars_BE_var = Dict(
    "g_vec" => g_vec,
    "s_vec" => s_vec,
    "phi_vec" => phi_vec,
    "oc_vec" => oc_vec,
    "of_vec" => of_vec,
    "COND_distr" => COND_distr,
    "lf_vec" => lf_vec,
    "lc_vec" => lc_vec,
    "l_value" => l_value,
    "l_vec" => l_vec,
    "Indic" => Indic,
    "cdf_g" => cdf_g
)
matwrite("BE_var_julia.mat", vars_BE_var)

# Save the second set of variables to BE_parameters.mat
vars_BE_parameters = Dict(
    "Cf" => Cf,
    "Cc" => Cc,
    "ALPHA" => ALPHA,
    "GAMMA" => GAMMA,
    "A" => A,
    "KAPPAf" => KAPPAf,
    "KAPPAc" => KAPPAc,
    "Pc" => Pc,
    "Pf" => Pf,
    "N" => N,
    "LN" => LN,
    "hired_lab_sh" => hired_lab_sh
)
matwrite("BE_parameters_julia.mat", vars_BE_parameters)

# Save the third set of variables to BE_values.mat
vars_BE_values = Dict(
    "VApw_BE" => VApw_BE,
    "AFS_BE" => AFS_BE,
    "landless_BE" => landless_BE
)
matwrite("BE_values_julia.mat", vars_BE_values)
```


## Government-mandated Land Reform
Second, we focus on the Government-mandated Land Reform (LR_main) part, using the data `BE_values_julia.mat` calculated in the previous sector, and define the `LR_main_eval` function:

```julia
function LR_main_eval(x, A)
    global KAPPAc, KAPPAf, LN, l_vec, Indic, g_vec, s_vec, phi_vec, Pc, Pf, GAMMA, ALPHA, Cf, Cc, N
    w = x[1]

     # Initialize vectors before in-place updating
     nf_vec = zeros(Float64, N)
     nc_vec = zeros(Float64, N)
     yf_vec = zeros(Float64, N)
     yc_vec = zeros(Float64, N)
     PIf_vec = zeros(Float64, N)
     PIc_vec = zeros(Float64, N)
     oc_vec = zeros(Int, N)
     of_vec = zeros(Int, N)

    # Solve Farmer's Problem
    nf_vec .= ((1 .- ALPHA) .* GAMMA .* (Pf ./ w) .* ((A .* KAPPAf .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA .* (1 .- ALPHA)))
    nc_vec .= ((1 .- ALPHA) .* GAMMA .* (Pc ./ w) .* ((A .* KAPPAc .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA .* (1 .- ALPHA)))
    yf_vec .= (A .* KAPPAf .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nf_vec .^ (1 .- ALPHA)) .^ GAMMA
    yc_vec .= (A .* KAPPAc .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nc_vec .^ (1 .- ALPHA)) .^ GAMMA
    PIf_vec .= (1 .- GAMMA .* (1 .- ALPHA)) .* Pf .* yf_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pf .* Cf .* ones(N)
    PIc_vec .= (1 .- GAMMA .* (1 .- ALPHA)) .* Pc .* yc_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pc .* Cc .* ones(N)

    # Determine technology choice vector
    oc_vec .= zeros(Int, N)
    max_value_Pw = max.(PIf_vec, w)
    I = findall(PIc_vec .> max_value_Pw)
    oc_vec[I] .= 1


    # Reform implied occupational choice vector
    of_vec .= zeros(Int, N)
    for ij = 1:minimum(I) .- 1
        of_vec[ij] = 1
    end

    Indic = float(Indic)
    of_vec = of_vec .* Indic

    # Remaining hired workers (=landless)
    Nww = sum((1 .- oc_vec .- of_vec)) ./ N
    f = (sum(of_vec .* nf_vec) ./ N) .+ (sum(oc_vec .* nc_vec) ./ N) .- Nww

    return f
end
```

The function `LR_main_eval` does the same thing as `BE_eval` but under the condition of government-mandated land reform.

Initial Guess:
```julia
guess  = [0.7]
```

Additional parameters:
```julia
params = [A] 
```


<p align="justify"> Then we use nlsolve function to find the numerical solutions of our non-linear system of equations. Given the market-clearing conditions calculated from LR_main_eval() and specified set of varaibles we are interested in, the algorithms uses numerical methods to approximate the roots of multiple equations simultaneously. We can derive the optimal distribution of lands and labor from the optimization results. </p>


```julia
result = nlsolve((res, x) -> res .= LR_main_eval(x, params), guess, show_trace=true, xtol=1e-16)
```



Solve farmer's problem:
```julia
nf_vec = ((1 .- ALPHA) .* GAMMA .* (Pf ./ w) .* ((A .* KAPPAf .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA * (1 .- ALPHA)))
nc_vec = ((1 .- ALPHA) .* GAMMA .* (Pc ./ w) .* ((A .* KAPPAc .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA * (1 .- ALPHA)))
yf_vec = (A .* KAPPAf .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nf_vec .^ (1 .- ALPHA)) .^ GAMMA
yc_vec = (A .* KAPPAc .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nc_vec .^ (1 .- ALPHA)) .^ GAMMA
PIf_vec = (1 .- GAMMA .* (1 .- ALPHA)) .* Pf .* yf_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pf .* Cf .* ones(length(N))
PIc_vec = (1 .- GAMMA .* (1 .- ALPHA)) .* Pc .* yc_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pc .* Cc .* ones(length(N))
```



## Market-based Land Reform
Third, we focus on the Market-based Land Reform (LR_market) part, using the data `BE_values_julia.mat` calculated in the previous sector, and define the `LR_market_eval` function:

```julia
function LR_market_eval(x, A)
    global KAPPAc, KAPPAf, LN, g_vec, s_vec, phi_vec, Pc, Pf, GAMMA, ALPHA, Cf, Cc, N, l_max, THETA

    w = x[1]
    q = x[2]

    # Factor price ratio
    qw_ratio = q ./ w

    # Solve unconstrained problem under each technology for all individuals
    lf_vec .= ((ALPHA ./ q) .^ ((1 .- (1 .- ALPHA) .* GAMMA) ./ (1 .- GAMMA)) .* ((1 .- ALPHA) ./ w) .^ (GAMMA .* (1 .- ALPHA) ./ (1 .- GAMMA)) .* (GAMMA .* Pf) .^ (1 ./ (1 .- GAMMA)) .* (A .* KAPPAf) .* g_vec)
    lc_vec .= ((ALPHA ./ q) .^ ((1 .- (1 .- ALPHA) .* GAMMA) ./ (1 .- GAMMA)) .* ((1 .- ALPHA) ./ w) .^ (GAMMA .* (1 .- ALPHA) ./ (1 .- GAMMA)) .* (GAMMA .* Pc) .^ (1 ./ (1 .- GAMMA)) .* (A .* KAPPAc) .* g_vec)
    nl_ratio = ((1 .- ALPHA) ./ ALPHA) .* qw_ratio
    nf_vec = nl_ratio .* lf_vec
    nc_vec = nl_ratio .* lc_vec
    yf_vec = (A .* KAPPAf .* s_vec) .^ (1 .- GAMMA) .* (lf_vec .^ ALPHA .* nf_vec .^ (1 - ALPHA)) .^ GAMMA
    yc_vec = (A .* KAPPAc .* s_vec) .^ (1 .- GAMMA) .* (lc_vec .^ ALPHA .* nc_vec .^ (1 - ALPHA)) .^ GAMMA
    PIf_vec = (1 .- GAMMA) .* Pf .* yf_vec .* (phi_vec .^ (1 .- GAMMA)) .- Cf .* ones(N)
    PIc_vec = (1 .- GAMMA) .* Pc .* yc_vec .* (phi_vec .^ (1 .- GAMMA)) .- Cc .* ones(N)

    # Find potentially constrained farmers
    cf_vec = zeros(Int, N)
    cc_vec = zeros(Int, N)
    If = findall(lf_vec .>= l_max)
    cf_vec[If] .= 1
    Ic = findall(lc_vec .>= l_max)
    cc_vec[Ic] .= 1

    # Land input demand accounting for constraint and its enforcement
    lfhat_vec = (1 .- cf_vec) .* lf_vec .+ cf_vec .* (THETA .* lf_vec .+ (1 .- THETA) .* l_max)
    lchat_vec = (1 .- cc_vec) .* lc_vec .+ cc_vec .* (THETA .* lc_vec .+ (1 .- THETA) .* l_max)

    # Auxiliary demand and profit functions that account for constrained farmers
    nf_max = (((1 .- ALPHA) .* GAMMA .* Pf .* (A .* KAPPAf .* g_vec) .^ (1 .- GAMMA) .* l_max .^ (GAMMA .* ALPHA)) ./ w) .^ (1 ./ (1 .- GAMMA .* (1 .- ALPHA)))
    nc_max = (((1 .- ALPHA) .* GAMMA .* Pc .* (A .* KAPPAc .* g_vec) .^ (1 .- GAMMA) .* l_max .^ (GAMMA .* ALPHA)) ./ w) .^ (1 ./ (1 .- GAMMA .* (1 .- ALPHA)))
    nfhat_vec = (1 .- cf_vec) .* nf_vec .+ cf_vec .* (THETA .* nf_vec .+ (1 .- THETA) .* nf_max)
    nchat_vec = (1 .- cc_vec) .* nc_vec .+ cc_vec .* (THETA .* nc_vec .+ (1 .- THETA) .* nc_max)
    yf_max = (A .* KAPPAf .* s_vec) .^ (1 .- GAMMA) .* (l_max .^ ALPHA .* nf_max .^ (1 .- ALPHA)) .^ GAMMA
    yc_max = (A .* KAPPAc .* s_vec) .^ (1 .- GAMMA) .* (l_max .^ ALPHA .* nc_max .^ (1 .- ALPHA)) .^ GAMMA
    yfhat_vec = (1 .- cf_vec) .* yf_vec .+ cf_vec .* (THETA .* yf_vec .+ (1 .- THETA) .* yf_max)
    ychat_vec = (1 .- cc_vec) .* yc_vec .+ cc_vec .* (THETA .* yc_vec .+ (1 .- THETA) .* yc_max)
    PIf_max = (1 .- GAMMA) .* Pf .* yf_max .* (phi_vec .^ (1 .- GAMMA)) .- Cf .* ones(N)
    PIc_max = (1 .- GAMMA) .* Pc .* yc_max .* (phi_vec .^ (1 .- GAMMA)) .- Cc .* ones(N)
    PIfhat_vec = (1 .- cf_vec) .* PIf_vec .+ cf_vec .* (THETA .* PIf_vec .+ (1 .- THETA) .* PIf_max)
    PIchat_vec = (1 .- cc_vec) .* PIc_vec .+ cc_vec .* (THETA .* PIc_vec .+ (1 .- THETA) .* PIc_max)

    # Solve for associated occupational choices
    ofhat_vec = zeros(Int, N)
    ochat_vec = zeros(Int, N)
    max_value_PIchat = max.(PIchat_vec, w)
    Ifhat = findall(PIfhat_vec .>= max_value_PIchat)
    ofhat_vec[Ifhat] .= 1
    max_value_PIfhat = max.(PIfhat_vec, w)
    Ichat = findall(PIchat_vec .> max_value_PIfhat)
    ochat_vec[Ichat] .= 1

    # Implied hired workers
    Nw = sum((1 .- ochat_vec .- ofhat_vec)) ./ N

    # Check whether labor and land markets clear
    f1 = (sum(ofhat_vec .* lfhat_vec) ./ N) .+ (sum(ochat_vec .* lchat_vec) ./ N) .- LN
    f2 = (sum(ofhat_vec .* nfhat_vec) ./ N) .+ (sum(ochat_vec .* nchat_vec) ./ N) .- Nw
    f = [f1, f2]

    return f
end
```


The function `LR_main_eval` does the same thing as `BE_eval` but under the condition of market-based reform in distributing the land above the government-mandated ceiling.

Initial Guess:
```julia
guess  = [0.88 0.09]
```

Additional parameters:
```julia
params = [A] 
```

<p align="justify"> Then we use nlsolve function to find the numerical solutions of our non-linear system of equations. Given the market-clearing conditions calculated from LR_market_eval() and specified set of varaibles we are interested in, the algorithms uses numerical methods to approximate the roots of multiple equations simultaneously. We can derive the optimal distribution of lands and labor from the optimization results. </p>

```julia
result = nlsolve(x -> LR_market_eval(x, A), guess) 
```

Extracting the solution:
```julia
w = result.zero[1]
q = result.zero[2]
```

Factor price ratio:
```julia
qw_ratio = q / w
```

Solve unconstrained problem under each technology for all individuals:
```julia
lf_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pf)^(1 / (1 - GAMMA))) * ((A * KAPPAf) .* g_vec)
lc_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pc)^(1 / (1 - GAMMA))) * ((A * KAPPAc) .* g_vec)
nl_ratio = ((1 - ALPHA) / ALPHA) * qw_ratio
nf_vec = nl_ratio .* lf_vec
nc_vec = nl_ratio .* lc_vec
yf_vec = (A * KAPPAf .* s_vec).^(1 - GAMMA) .* (lf_vec.^ALPHA .* nf_vec.^(1 - ALPHA)).^GAMMA
yc_vec = (A * KAPPAc .* s_vec).^(1 - GAMMA) .* (lc_vec.^ALPHA .* nc_vec.^(1 - ALPHA)).^GAMMA
PIf_vec = (1 - GAMMA) * Pf * yf_vec .* (phi_vec.^(1 - GAMMA)) - Cf * ones(N)
PIc_vec = (1 - GAMMA) * Pc * yc_vec .* (phi_vec.^(1 - GAMMA)) - Cc * ones(N)
```

Initialize constraint vectors for farmers:
```julia
cf_vec = zeros(N)
cc_vec = zeros(N)
```

Find potentially constrained farmers:
```julia
If = findall(lf_vec .>= l_max)
cf_vec[If] .= 1
Ic = findall(lc_vec .>= l_max)
cc_vec[Ic] .= 1
```

Land input demand accounting for constraint and its enforcement:
```julia
lfhat_vec = (1 .- cf_vec) .* lf_vec .+ cf_vec .* (THETA .* lf_vec .+ (1 .- THETA) .* l_max)
lchat_vec = (1 .- cc_vec) .* lc_vec .+ cc_vec .* (THETA .* lc_vec .+ (1 .- THETA) .* l_max)
```

Auxiliary demand and profit functions that account for constrained farmers:
```julia
nf_max = (((1 - ALPHA) * GAMMA * Pf * (A * KAPPAf .* g_vec) .^ (1 - GAMMA) .* l_max .^ (GAMMA * ALPHA)) / w) .^ (1 / (1 - GAMMA * (1 - ALPHA)))
nc_max = (((1 - ALPHA) * GAMMA * Pc * (A * KAPPAc .* g_vec) .^ (1 - GAMMA) .* l_max .^ (GAMMA * ALPHA)) / w) .^ (1 / (1 - GAMMA * (1 - ALPHA)))
nfhat_vec = (1 .- cf_vec) .* nf_vec + cf_vec .* (THETA .* nf_vec + (1 - THETA) .* nf_max)
nchat_vec = (1 .- cc_vec) .* nc_vec + cc_vec .* (THETA .* nc_vec + (1 - THETA) .* nc_max)
yf_max = (A * KAPPAf .* s_vec) .^ (1 - GAMMA) .* (l_max .^ ALPHA .* nf_max .^ (1 - ALPHA)) .^ GAMMA
yc_max = (A * KAPPAc .* s_vec) .^ (1 - GAMMA) .* (l_max .^ ALPHA .* nc_max .^ (1 - ALPHA)) .^ GAMMA
yfhat_vec = (1 .- cf_vec) .* yf_vec + cf_vec .* (THETA .* yf_vec + (1 - THETA) .* yf_max)
ychat_vec = (1 .- cc_vec) .* yc_vec + cc_vec .* (THETA .* yc_vec + (1 - THETA) .* yc_max)
PIf_max = (1 - GAMMA) * Pf * yf_max .* (phi_vec .^ (1 - GAMMA)) - Cf * ones(N)
PIc_max = (1 - GAMMA) * Pc * yc_max .* (phi_vec .^ (1 - GAMMA)) - Cc * ones(N)
PIfhat_vec = (1 .- cf_vec) .* PIf_vec + cf_vec .* (THETA .* PIf_vec + (1 - THETA) .* PIf_max)
PIchat_vec = (1 .- cc_vec) .* PIc_vec + cc_vec .* (THETA .* PIc_vec + (1 - THETA) .* PIc_max)
```

Solve for associated occupational choices:
```julia
ofhat_vec = zeros(N)
ochat_vec = zeros(N)
Ifhat = findall(PIfhat_vec .>= max.(PIchat_vec, w))
ofhat_vec[Ifhat] .= 1
Ichat = findall(PIchat_vec .> max.(PIfhat_vec, w))
ochat_vec[Ichat] .= 1
```


Implied hired workers:
```julia
Nw = sum((1 .- ochat_vec .- ofhat_vec)) ./ N 
```

Check whether labor and land markets clear:
```julia
f1 = (sum(ofhat_vec .* lfhat_vec) / N) + (sum(ochat_vec .* lchat_vec) / N) - LN 
f2 = (sum(ofhat_vec .* nfhat_vec) / N) + (sum(ochat_vec .* nchat_vec) / N) - Nw 
f = [f1, f2]
```





# Comparison with the original results
In this section, we verify the replication results by comparing our figures with those in the original paper.

## BE Farm Size Distribution
![BE Farm Size Distribution](BE_Farm%20Size%20Distribution%20Across%20Specified%20Bins.png)

## BE Hired Labor per Hectare Distribution Across Specified Bins.png
![BE Hired Labor per Hectare Distribution](BE_Hired%20Labor%20per%20Hectare%20Distribution%20Across%20Specified%20Bins.png)

## BE Land Size Distribution
![BE Land Size Distribution](BE_Land%20Size%20Distribution%20Across%20Specified%20Bins.png)

## BE Value Added per Worker Distribution
![BE Value Added per Worker Distribution](BE_Value%20Added%20per%20Worker%20Distribution%20Across%20Specified%20Bins.png)

## LR Main Farm Size Distribution 
![LR Main Farm Size Distribution](LR_main_Farm%20Size%20Distribution%20Across%20Specified%20Bins.png)

## LR Main Hired Labor per Hectare Distribution
![LR Main Hired Labor per Hectare Distribution](LR_main_Hired%20Labor%20per%20Hectare%20Distribution%20Across%20Specified%20Bins.png)

## LR Main Land Size Distribution
![LR Main Land Size Distribution](LR_main_Land%20Size%20Distribution%20Across%20Specified%20Bins.png)

## LR Main Value Added per Worker Distribution
![LR Main Value Added per Worker Distribution](LR_main_Value%20Added%20per%20Worker%20Distribution%20Across%20Specified%20Bins.png)

## LR Market Farm Size Distribution
![LR Market Farm Size Distribution](LR_market_Farm%20Size%20Distribution%20Across%20Specified%20Bins.png)

## LR Market Hired Labor per Hectare Distribution
![LR Market Hired Labor per Hectare Distribution](LR_market_Hired%20Labor%20per%20Hectare%20Distribution%20Across%20Specified%20Bins.png)

## LR Market Land Size Distribution
![LR Market Land Size Distribution](LR_market_Land%20Size%20Distribution%20Across%20Specified%20Bins.png)

## LR Market Value Added per Worker Distribution
![LR Market Value Added per Worker Distribution](LR_market_Value%20Added%20per%20Worker%20Distribution%20Across%20Specified%20Bins.png)








# Conclusion




# Acknowledgements
<p align="justify"> We extend our heartfelt appreciation to Professor Florian Oswald for his exceptional mentorship and dedication to our academic growth. His teaching in Introduction to Programming last year and Computational Economics this year has been nothing short of transformative. Professor Oswald's clarity of instruction, insightful guidance, and unwavering support have enabled us to grasp complex concepts and excel in our studies. We are immensely grateful for the knowledge and skills we have gained under his tutelage. </p>




