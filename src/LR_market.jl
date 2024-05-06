# ReplicationPackage_KexinCHEN_JinxuMI.jl
# Replicating the article "Land Reform and Productivity: 
# A Quantitative Analysis with Micro Data" (Adamopoulos & Restuccia, 2020)

# LR_market.jl
# Market-based land reform.

# Import Julia packages
using Pkg
Pkg.add("Printf")
Pkg.add("MAT")
Pkg.add("Distributions")
Pkg.add("LinearAlgebra")
Pkg.add("Statistics")
Pkg.add("DelimitedFiles")
Pkg.add("Optim")
Pkg.add("NLsolve")
Pkg.add("Plots")
Pkg.add("Serialization")
Pkg.add("CSV")
Pkg.add("DataFrames")

using Printf
using MAT # For handling .mat files
using Distributions
using LinearAlgebra
using Statistics
using DelimitedFiles
using Optim
using NLsolve
using Plots
using Serialization
using CSV
using DataFrames

# Load parameters and values from benchmark economy
file_BE_var = matopen("BE_var_julia.mat")
BE_var = read(file_BE_var)
close(file_BE_var)

g_vec = BE_var["g_vec"]
s_vec = BE_var["s_vec"]
phi_vec = BE_var["phi_vec"]
oc_vec = BE_var["oc_vec"]
of_vec = BE_var["of_vec"]
COND_distr = BE_var["COND_distr"]
lf_vec = BE_var["lf_vec"]
lc_vec = BE_var["lc_vec"]
l_value = BE_var["l_value"]
l_vec = BE_var["l_vec"]
Indic = BE_var["Indic"]
cdf_g = BE_var["cdf_g"]


file_BE_parameters = matopen("BE_parameters_julia.mat")
BE_parameters = read(file_BE_parameters)
close(file_BE_parameters)

Cf = BE_parameters["Cf"]
Cc = BE_parameters["Cc"]
ALPHA = BE_parameters["ALPHA"]
GAMMA = BE_parameters["GAMMA"]
A = BE_parameters["A"]
KAPPAf = BE_parameters["KAPPAf"]
KAPPAc = BE_parameters["KAPPAc"]
Pc = BE_parameters["Pc"]
Pf = BE_parameters["Pf"]
N = BE_parameters["N"]
LN = BE_parameters["LN"]
hired_lab_sh = BE_parameters["hired_lab_sh"]


file_BE_values = matopen("BE_values_julia.mat")
BE_values = read(file_BE_values)
close(file_BE_values)

VApw_BE = BE_values["VApw_BE"]
AFS_BE = BE_values["AFS_BE"]
landless_BE = BE_values["landless_BE"]


%#DDITIONAL FACTORS
#A=1.1;
#KAPPAf = KAPPAf*1.65;
#KAPPAc = KAPPAc*1.65;

#Parameters of the market-based land reform
#--------------------------------------------------------------------------
#ceiling
l_max = 5
#degree of enforcement: a fraction THETA of potentially constrained farmers gets to keep their land
THETA = 0


# Function LR_market_eval
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

guess  = [0.88 0.09]

# Additional parameters, assuming 'A' is needed as demonstrated in the previous explanation
params = [A] 


# Solving the system
result = nlsolve(x -> LR_market_eval(x, A), guess) 

w = result.zero[1]
q = result.zero[2]

# Factor price ratio
qw_ratio = q / w

# Solve unconstrained problem under each technology for all individuals
lf_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pf)^(1 / (1 - GAMMA))) * ((A * KAPPAf) .* g_vec)
lc_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pc)^(1 / (1 - GAMMA))) * ((A * KAPPAc) .* g_vec)
nl_ratio = ((1 - ALPHA) / ALPHA) * qw_ratio
nf_vec = nl_ratio .* lf_vec
nc_vec = nl_ratio .* lc_vec
yf_vec = (A * KAPPAf .* s_vec).^(1 - GAMMA) .* (lf_vec.^ALPHA .* nf_vec.^(1 - ALPHA)).^GAMMA
yc_vec = (A * KAPPAc .* s_vec).^(1 - GAMMA) .* (lc_vec.^ALPHA .* nc_vec.^(1 - ALPHA)).^GAMMA
PIf_vec = (1 - GAMMA) * Pf * yf_vec .* (phi_vec.^(1 - GAMMA)) - Cf * ones(N)
PIc_vec = (1 - GAMMA) * Pc * yc_vec .* (phi_vec.^(1 - GAMMA)) - Cc * ones(N)

# Initialize constraint vectors for farmers
cf_vec = zeros(N)
cc_vec = zeros(N)

# Find potentially constrained farmers
If = findall(lf_vec .>= l_max)
cf_vec[If] .= 1
Ic = findall(lc_vec .>= l_max)
cc_vec[Ic] .= 1

# Land input demand accounting for constraint and its enforcement
lfhat_vec = (1 .- cf_vec) .* lf_vec .+ cf_vec .* (THETA .* lf_vec .+ (1 .- THETA) .* l_max)
lchat_vec = (1 .- cc_vec) .* lc_vec .+ cc_vec .* (THETA .* lc_vec .+ (1 .- THETA) .* l_max)

# Auxiliary demand and profit functions that account for constrained farmers
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

# Solve for associated occupational choices
ofhat_vec = zeros(N)
ochat_vec = zeros(N)
Ifhat = findall(PIfhat_vec .>= max.(PIchat_vec, w))
ofhat_vec[Ifhat] .= 1
Ichat = findall(PIchat_vec .> max.(PIfhat_vec, w))
ochat_vec[Ichat] .= 1


# Implied hired workers
Nw = sum((1 .- ochat_vec .- ofhat_vec)) ./ N 


# Check whether labor and land markets clear
f1 = (sum(ofhat_vec .* lfhat_vec) / N) + (sum(ochat_vec .* lchat_vec) / N) - LN 
f2 = (sum(ofhat_vec .* nfhat_vec) / N) + (sum(ochat_vec .* nchat_vec) / N) - Nw 
f = [f1, f2]


#--------------------------------------------------------------------------
#Calculate distributions and stats of interest
#--------------------------------------------------------------------------

# Share of food farm operators
Nf = sum(ofhat_vec) / N

# Share of cach farm operators
Nc = sum(ochat_vec) / N

# Land input for each individual
l_vec = ofhat_vec .* lfhat_vec + ochat_vec .* lchat_vec

# Output input for each individual
y_vec = ofhat_vec .* yfhat_vec + ochat_vec .* ychat_vec

# Labor input for each individual
n_vec = ofhat_vec .* nfhat_vec + ochat_vec .* nchat_vec

# Conditional (on operating) distribution
Indic = ochat_vec + ofhat_vec
COND_distr = Indic / sum(Indic)  

# Average Farm Size
l_value = sum(COND_distr .* lfhat_vec .* ofhat_vec) + sum(COND_distr .* lchat_vec .* ochat_vec)
AFS_post = l_value

# Bins for distribution
bin1 = 1
bin2 = 2
bin3 = 5
bin4 = 7
bin5 = 10
bin6 = 15

# Cutoffs for bins
lhat1in = findall(l_vec .<= bin1)
lhat2in = findall((l_vec .<= bin2) .& (l_vec .> bin1))
lhat3in = findall((l_vec .<= bin3) .& (l_vec .> bin2))
lhat4in = findall((l_vec .<= bin4) .& (l_vec .> bin3))
lhat5in = findall((l_vec .<= bin5) .& (l_vec .> bin4))
lhat6in = findall((l_vec .<= bin6) .& (l_vec .> bin5))
lhat7in = findall(l_vec .> bin6)

#Compute distributions: farms, land, output, hired labor
if isempty(lhat1in)==1
    farm1   = 0
    land1   = 0
    output1 = 0
    hirelab1=0
else
    farm1   = sum(COND_distr[lhat1in])
    land1   = sum(l_vec[lhat1in].*COND_distr[lhat1in])/l_value
    output1 = sum(y_vec[lhat1in].*COND_distr[lhat1in])
    hirelab1 = sum(n_vec[lhat1in].*COND_distr[lhat1in])
end

if isempty(lhat2in)
    farm2 = 0
    land2 = 0
    output2 = 0
    hirelab2 = 0
else
    farm2 = sum(COND_distr[lhat2in])
    land2 = sum(l_vec[lhat2in] .* COND_distr[lhat2in]) / l_value
    output2 = sum(y_vec[lhat2in] .* COND_distr[lhat2in])
    hirelab2 = sum(n_vec[lhat2in] .* COND_distr[lhat2in])
end

if isempty(lhat3in)
    farm3 = 0
    land3 = 0
    output3 = 0
    hirelab3 = 0
else
    farm3 = sum(COND_distr[lhat3in])
    land3 = sum(l_vec[lhat3in] .* COND_distr[lhat3in]) / l_value
    output3 = sum(y_vec[lhat3in] .* COND_distr[lhat3in])
    hirelab3 = sum(n_vec[lhat3in] .* COND_distr[lhat3in])
end

if isempty(lhat4in)
    farm4 = 0
    land4 = 0
    output4 = 0
    hirelab4 = 0
else
    farm4 = sum(COND_distr[lhat4in])
    land4 = sum(l_vec[lhat4in] .* COND_distr[lhat4in]) / l_value
    output4 = sum(y_vec[lhat4in] .* COND_distr[lhat4in])
    hirelab4 = sum(n_vec[lhat4in] .* COND_distr[lhat4in])
end

if isempty(lhat5in)
    farm5 = 0
    land5 = 0
    output5 = 0
    hirelab5 = 0
else
    farm5 = sum(COND_distr[lhat5in])
    land5 = sum(l_vec[lhat5in] .* COND_distr[lhat5in]) / l_value
    output5 = sum(y_vec[lhat5in] .* COND_distr[lhat5in])
    hirelab5 = sum(n_vec[lhat5in] .* COND_distr[lhat5in])
end

if isempty(lhat6in)
    farm6 = 0
    land6 = 0
    output6 = 0
    hirelab6 = 0
else
    farm6 = sum(COND_distr[lhat6in])
    land6 = sum(l_vec[lhat6in] .* COND_distr[lhat6in]) / l_value
    output6 = sum(y_vec[lhat6in] .* COND_distr[lhat6in])
    hirelab6 = sum(n_vec[lhat6in] .* COND_distr[lhat6in])
end

if isempty(lhat7in)
    farm7 = 0
    land7 = 0
    output7 = 0
    hirelab7 = 0
else
    farm7 = sum(COND_distr[lhat7in])
    land7 = sum(l_vec[lhat7in] .* COND_distr[lhat7in]) / l_value
    output7 = sum(y_vec[lhat7in] .* COND_distr[lhat7in])
    hirelab7 = sum(n_vec[lhat7in] .* COND_distr[lhat7in])
end


# Distributions of output, hired labor, farms, farmers in the model for specified bins
l_model = [land1, land2, land3, land4, land5, land6, land7]
y_model = [output1, output2, output3, output4, output5, output6, output7]
h_model = [hirelab1, hirelab2, hirelab3, hirelab4, hirelab5, hirelab6, hirelab7]
f_model = [farm1, farm2, farm3, farm4, farm5, farm6, farm7]
n_model = h_model + f_model


# Benchmark prices for crops
Pf_star = 1
Pc_star = 1

# Aggregate Output per capita
VApw = (Pc_star * sum(ochat_vec .* ychat_vec) + Pf_star * sum(ofhat_vec .* yfhat_vec)) / N
VApw_post = VApw

# Landless
landless_post = Nw

# Aggregate hired labor to land ratio
AggHLph = sum(h_model) / sum(l_model)

# Value added per worker (hired + operators) across specified bins in model
VApw_model = (y_model ./ n_model) / VApw

# Hired labor per hectare across specified bins in model
HLph_model = (h_model ./ l_model) / AggHLph


#DISTRIBUTIONS IN THE DATA ACROSS SPECIFIED BINS (2003-04 from panel sample)
#Range in Ha      0--1         1--2            2--5            5--7            7--10           10-15        15+      
#                 bin1         bin2            bin3            bin4             bin5            bin6        bin7


# Farm size distribution in 2003-04
farm_pdf_03panel = [0.347, 0.275, 0.222, 0.048, 0.048, 0.036, 0.024]

# Distribution of land in 2003-04
land_pdf_03panel = [0.067, 0.137, 0.229, 0.093, 0.137, 0.157, 0.180]

# Distribution of value added per worker in 2003-04
VApw_2003data = [0.49359913, 0.844372155, 0.712228608, 1.451573247, 1.284452801, 0.833275073, 1.331615042]

# Distribution of hired labor in 2003-04
HLph_data_03 = [1.068650632, 1.135513794, 0.943871958, 0.861784207, 0.513856394, 1.161874909, 1.307819784]


#Figures
#--------------------------------------------------------------------------
Ha = [1, 2, 3, 4, 5, 6, 7]
xticklabels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
begin
    farm_pdf_03panel_matrix = reshape(farm_pdf_03panel, :, 1)
    data = hcat(reshape(l_model, :, 1), farm_pdf_03panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p1 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 1.0),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Farms")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end

    savefig(p1, "LR_market_Farm Size Distribution Across Specified Bins.png")
end

begin
    land_pdf_03panel_matrix = reshape(land_pdf_03panel, :, 1)
    data = hcat(reshape(l_model, :, 1), land_pdf_03panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p2 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 1.0),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Land")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p2, "LR_market_Land Size Distribution Across Specified Bins.png")
end

begin
    VApw_model_col = reshape(VApw_model, :, 1)
    VApw_data_col = reshape(VApw_2003data, :, 1)
    data = hcat(VApw_model_col, VApw_data_col)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p3 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 2),
              xlabel="Farm Size Class in Ha", ylabel="Value added per worker")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p3, "LR_market_Value Added per Worker Distribution Across Specified Bins.png")

end

begin
    HLph_model_col = reshape(HLph_model, :, 1)
    HLph_data_col = reshape(HLph_data_03, :, 1)
    data = hcat(HLph_model_col, HLph_data_col)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p4 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 1.5),
              xlabel="Farm Size Class in Ha", ylabel="Value added per worker")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p4, "LR_market_Hired Labor per Hectare Distribution Across Specified Bins.png")

end


# Aggregate Output Per Hired Worker
VAphw = VApw / Nw

# All labor (hired+operators) by crop
LAB_f = (sum(ofhat_vec .* nfhat_vec) / N) + Nf
LAB_c = (sum(ochat_vec .* nchat_vec) / N) + Nc

# Crop Outputs
Yc = sum(ochat_vec .* ychat_vec) / N
Yf = sum(ofhat_vec .* yfhat_vec) / N

# Labor productivity by crop
VApw_cash = Yc / LAB_c
VApw_food = Yf / LAB_f

# Total effects of land reform (Table 7, second column results)
println("\n% CHANGE IN AVERAGE FARM SIZE FOLLOWING THE REFORM")
AFS_EFF = 100 * (AFS_post - AFS_BE) / AFS_BE

println("\n% CHANGE IN VALUE ADDED PER WORKER FOLLOWING THE REFORM")
VApw_EFF = 100 * (VApw_post - VApw_BE) / VApw_BE

println("\nCHANGE IN % OF LANDLESS FOLLOWING THE REFORM")
landless_CHANGE = 100 * (landless_post - landless_BE)

