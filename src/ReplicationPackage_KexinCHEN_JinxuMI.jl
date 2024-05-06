module ReplicationPackage_KexinCHEN_JinxuMI

# Write your package code here.

# ReplicationPackage_KexinCHEN_JinxuMI.jl
# Replicating the article "Land Reform and Productivity: 
# A Quantitative Analysis with Micro Data" (Adamopoulos & Restuccia, 2020)

# BE.jl 
# Simulation/Calibration of Benchmark Economy (BE).

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

# Define a structure with parameters and variables
struct BenchmarkEconomy
    AFS::Float64
    hired_lab_sh::Float64
    cash_oper_sh::Float64
    farm_pdf_85panel::Vector{Float64}
    YNR_cf_target::Float64
    VAR_lTFPR_data::Float64
    CORR_lTFP_lTFPR_data::Float64
    land_pdf_85panel::Vector{Float64}
    VAR_lTFP_data::Float64
    GAMMA::Float64
    ALPHA::Float64
    A::Float64
    KAPPAf::Float64
    Pf::Float64
    Pc::Float64
    LN::Float64
    KAPPAc::Float64
end

# Initialize the BenchmarkEconomy with default values
function init_benchmark_economy()

    AFS = 3.7
    hired_lab_sh = 0.611
    LN = AFS * (1 - hired_lab_sh) 

    BenchmarkEconomy(
        AFS,
        hired_lab_sh,
        0.617, # cash_oper_sh
        [0.168, 0.269, 0.347, 0.090, 0.078, 0.018, 0.030], # farm_pdf_85panel
        2.95, # YNR_cf_target
        1.01^2, # VAR_lTFPR_data
        0.97, # CORR_lTFP_lTFPR_data
        [0.036, 0.114, 0.304, 0.143, 0.184, 0.061, 0.158], # land_pdf_85panel
        1.15^2, # VAR_lTFP_data
        0.7, # GAMMA
        0.3, # ALPHA
        1, # A
        1, # KAPPAf
        1, # Pf
        1, # Pc
        LN, # Calculated LN
        1.25 # KAPPAc
    )
end




#--------------------------------------------------------------------------
#Targets in Calibration
#--------------------------------------------------------------------------

# Define the calibration targets and parameters
AFS = 3.7 # Average Farm Size in 1984-85 (Total Land Area/Total Number of Holdings) - from panel sample
hired_lab_sh = 0.611 # Share of hired labor in total labor
cash_oper_sh = 0.617 # Share of cash crop operators in total operators

# Farm size Distribution in 1984-85 (from panel sample)
farm_pdf_85panel = [0.168, 0.269, 0.347, 0.090, 0.078, 0.018, 0.030]

# Ratio of average productivities between cash crop farms and food crop farms
YNR_cf_target = 2.95

# Variance of land quality adjusted log-distortions (TFPR) for 1985 (panel data)
VAR_lTFPR_data = 1.01^2

# Correlation of land quality adjusted log-TFP and log-TFPR for 1985
CORR_lTFP_lTFPR_data = 0.97




#--------------------------------------------------------------------------
#Additional Data Moments (not targeted)
#--------------------------------------------------------------------------

# Additional Data Moments (not targeted)
land_pdf_85panel = [0.036, 0.114, 0.304, 0.143, 0.184, 0.061, 0.158]

# Variance of land quality adjusted log-TFP from data (among all farms in data)
VAR_lTFP_data = 1.15^2 # Panel for 1985




#--------------------------------------------------------------------------
#Parameters
#--------------------------------------------------------------------------

# Parameters
GAMMA = 0.7
ALPHA = 0.3
A = 1 # Normalization
KAPPAf = 1 # Normalization
Pf = 1 # Normalization
Pc = 1 # Normalization
LN = AFS * (1 - hired_lab_sh) # Choose land per capita to match AFS
KAPPAc = 1.25 # Chosen to match closely the cash-food productivity ratio

# Output some of the parameters to verify
println("AFS: ", AFS)
println("hired_lab_sh: ", hired_lab_sh)
println("cash_oper_sh: ", cash_oper_sh)
println("farm_pdf_85panel: ", farm_pdf_85panel)
println("YNR_cf_target: ", YNR_cf_target)
println("VAR_lTFPR_data: ", VAR_lTFPR_data)
println("CORR_lTFP_lTFPR_data: ", CORR_lTFP_lTFPR_data)
println("land_pdf_85panel: ", land_pdf_85panel)
println("VAR_lTFP_data: ", VAR_lTFP_data)
println("GAMMA: ", GAMMA)
println("ALPHA: ", ALPHA)
println("LN (Land per capita): ", LN)
println("KAPPAc: ", KAPPAc)




#--------------------------------------------------------------------------
#Matrix of correlated data on abilities (s) and distortions (phi) - generated from joint log-normal distribution
#--------------------------------------------------------------------------

# Number of data points (individuals in the economy)
# N = 100000

# Normalize vector of means for ability and distortions distribution
# MU_vec = [0.0, 0.0]

# Specify the variance-covariance matrix of joint ability-distortions distribution
# SIGMA_mat = [13.421 -11.586; -11.586 11.771]

# Generate N x 2 matrix of N data points for each of the two variables s, phi from multivariate normal distribution with above mean and var-cov matrix
# dist = MvNormal(MU_vec, SIGMA_mat)
# Ngen = rand(dist, N)'

# Load the .mat file
matfile_AbDist = matopen("AbDist_matrix.mat")

Ngen = read(matfile_AbDist, "Ngen")

# Exponential transformation
LNgen = exp.(Ngen)

# Close the file
close(matfile_AbDist)

# Extract the first column of the generated data as the s vector and the second column as the phi vector (unsorted).
unsorted_s_vec = LNgen[:, 1]
unsorted_phi_vec = LNgen[:, 2]

# Combine to generate the unsorted g vector
unsorted_g_vec = unsorted_s_vec .* unsorted_phi_vec

# Get the indices that would sort the array
idx = sortperm(unsorted_g_vec)

# Sort the g vector in ascending order and get the sorted vector
g_vec = sort(unsorted_g_vec)

# Sort s and phi vectors in ascending order of g
s_vec = unsorted_s_vec[idx]
phi_vec = unsorted_phi_vec[idx]

# Number of data points (individuals in the economy)
N = 100000

# Vector of 1 for each individual
pdf_g = ones(N)

# The weight of each observation is 1/N
sh_g = pdf_g / N

# The sum of shares should be 1
cdf_g = cumsum(sh_g)




#--------------------------------------------------------------------------
#Solve equilibrium to determine fixed costs Cf, Cc that match (i) share of
#operators; (ii) share of operators in cash crops
#--------------------------------------------------------------------------

# Function BE_eval
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

# Your initial guess for the variables
guess = [-2.4, -2.0]


#--------------------------------------------------------------------------
#Targets in Calibration
#--------------------------------------------------------------------------

#Average Farm Size in 1984-85 (Total Land Area/Total Number of Holdings) - from panel sample
AFS             = 3.7;

#Share of hired labor in total labor
hired_lab_sh    = 0.611;

#Share of cash crop operators in total operators
cash_oper_sh    = 0.617;

#Farm size Distribution in 1984-85 (from panel sample) 
#Range in Ha      0--1         1--2            2--5            5--7            7--10           10-15        15+      
#                 bin1         bin2            bin3            bin4             bin5            bin6        bin7
#Farms (%)       0.168       0.269           0.347           0.090             0.078           0.048      0.030 
farm_pdf_85panel   = [0.168 0.269 0.347 0.090 0.078 0.018 0.030];

#Ratio of average productivities between cash crop farms and food crop farms
YNR_cf_target   = 2.95;

#Variance of land quality adjusted log-distortions (TFPR) for 1985 (panel data)
VAR_lTFPR_data  = 1.01^2;

#Correlation of land quality adjusted log-TFP and log-TFPR for 1985
CORR_lTFP_lTFPR_data = 0.97;


#--------------------------------------------------------------------------
#Additional Data Moments (not targeted)
#--------------------------------------------------------------------------

#Land Distribution in 1984-85 (from panel sample) 
#Range in Ha      0--1         1--2            2--5            5--7            7--10           10-15        15+      
#                 bin1         bin2            bin3            bin4             bin5            bin6        bin7
#LAND  (%)       0.036        0.114           0.304           0.143             0.184           0.220    0.158
land_pdf_85panel   = [0.036 0.114 0.304 0.143 0.184 0.061 0.158];

#Variance of land quality adjusted log-TFP from data (among all farms in data)
VAR_lTFP_data = 1.15^2; #panel for 1985


#--------------------------------------------------------------------------
#Parameters
#--------------------------------------------------------------------------
GAMMA       = 0.7;
ALPHA       = 0.3;
A           = 1; #normalization
KAPPAf      = 1; #normalization
Pf          = 1; #normalization
Pc          = 1; #normalization 
LN          = AFS*(1-hired_lab_sh); #choose land per capita to match AFS 
KAPPAc      = 1.25; #chosen to match closely the cash-food productivity ratio


# Additional parameters, assuming 'A' is needed as demonstrated in the previous explanation
params = [A] 

# Solve the system, using an anonymous function to pass extra parameters to BE_eval
result = nlsolve((res, x) -> res .= BE_eval(x, params), guess, show_trace=true, xtol=1e-16)

# Extracting the solution
x = result.zero
Cf = x[1]
Cc = x[2]

# Check for convergence
converged = result.f_converged
println("Converged: ", converged)




# Compute occupational choice cutoffs
# choose g_lbar to match a share of hired labor in total labor
INl = findfirst(cdf_g .> hired_lab_sh)
g_lbar = g_vec[INl]
sh_l = cdf_g[INl]
g_lbar_Indic = findfirst(x -> x > g_lbar, g_vec) - 1

# choose g_ubar to match a share of food crop operators in total operators
INu = findfirst(cdf_g .> (hired_lab_sh + (1 - hired_lab_sh) * (1 - cash_oper_sh)))
g_ubar = g_vec[INu]
sh_u = cdf_g[INu]
g_ubar_Indic = findfirst(x -> x > g_ubar, g_vec) - 1

# Factor prices - from model equations
w = (Cc - Cf) / ((g_ubar / g_lbar) * ((Pc / Pf)^(1 / (1 - GAMMA)) * (KAPPAc / KAPPAf) - 1)) - Cf
q = ALPHA * (((g_lbar * (GAMMA^(GAMMA / (1 - GAMMA))) * (1 - GAMMA) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * (Pf)^(1 / (1 - GAMMA)) * (A * KAPPAf)) / (w + Cf))^((1 - GAMMA) / (ALPHA * GAMMA)))
qw_ratio = q / w

# Initialize of_vec with ones and then set specific elements to zero based on conditions
of_vec = ones(Int, N)
of_vec[1:g_lbar_Indic-1] .= 0
of_vec[g_ubar_Indic+1:N] .= 0

# Initialize oc_vec with ones and then set specific elements to zero based on conditions
oc_vec = ones(Int, N)
oc_vec[1:g_ubar_Indic] .= 0


# Solve problem under each technology for every individual
lf_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pf)^(1 / (1 - GAMMA))) * (A * KAPPAf) .* g_vec
lc_vec = ((ALPHA / q)^((1 - (1 - ALPHA) * GAMMA) / (1 - GAMMA))) * (((1 - ALPHA) / w)^(GAMMA * (1 - ALPHA) / (1 - GAMMA))) * ((GAMMA * Pc)^(1 / (1 - GAMMA))) * (A * KAPPAc) .* g_vec
nl_ratio = ((1 - ALPHA) / ALPHA) * qw_ratio
nf_vec = nl_ratio .* lf_vec
nc_vec = nl_ratio .* lc_vec
yf_vec = (A * KAPPAf .* s_vec) .^ (1 - GAMMA) .* (lf_vec .^ ALPHA .* nf_vec .^ (1 - ALPHA)) .^ GAMMA
yc_vec = (A * KAPPAc .* s_vec) .^ (1 - GAMMA) .* (lc_vec .^ ALPHA .* nc_vec .^ (1 - ALPHA)) .^ GAMMA
PIf_vec = (1 - GAMMA) * Pf * yf_vec .* (phi_vec .^ (1 - GAMMA)) - Cf .* ones(N)
PIc_vec = (1 - GAMMA) * Pc * yc_vec .* (phi_vec .^ (1 - GAMMA)) - Cc .* ones(N)

# Compute aggregates using occupational choice vectors
# Labor Shares - based on occupational choices
Nw = sum(1 .- oc_vec .- of_vec) / N  # Share of hired labor
Nf = sum(of_vec) / N  # Share of food farm operators
Nc = sum(oc_vec) / N  # Share of cash farm operators
LAB_f = (sum(of_vec .* nf_vec) / N) + Nf  # Total share of labor in food farms (hired + operators)
LAB_c = (sum(oc_vec .* nc_vec) / N) + Nc  # Total share of labor in cash farms (hired + operators)

# Crop Outputs - based on occupational choices
Yc = sum(oc_vec .* yc_vec) / N  # Total output of cash crop farms
Yf = sum(of_vec .* yf_vec) / N  # Total output of food crop farms

# Ratio of cash-to-food labor productivities (INCLUDES OPERATORS)
YNR_cf_model = (Yc / LAB_c) / (Yf / LAB_f)

# Clear land and labor markets
f1 = (sum(of_vec .* lf_vec) / N) + (sum(oc_vec .* lc_vec) / N) - LN  # Land market clearing condition
f2 = (sum(of_vec .* nf_vec) / N) + (sum(oc_vec .* nc_vec) / N) - Nw  # Labor market clearing condition
f = [f1, f2]  # Construct the vector f from f1 and f2




#--------------------------------------------------------------------------
#Production choices for all individuals (including non-operators)
#--------------------------------------------------------------------------

# Farm size vector (land input)
l_vec = of_vec .* lf_vec + oc_vec .* lc_vec

# Hired labor vector
n_vec = of_vec .* nf_vec + oc_vec .* nc_vec

# Output vector
y_vec = of_vec .* yf_vec + oc_vec .* yc_vec




#--------------------------------------------------------------------------
#Calculate distributions and statistics of interest
#--------------------------------------------------------------------------

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




#---------------CONDITIONAL (ON OPERATING) DISTRIBUTIONS OVER SPECIFIED BINS---------------------- 

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




###########################################################################
#Figures
###########################################################################

# BINS (horizontal axis points)
Ha = [1, 2, 3, 4, 5, 6, 7]

#---------DISTRIBUTIONS IN THE DATA----------------------------------------
#Value added per worker across specified bins - DATA
VApw_data = [0.426269805, 0.474383459, 0.572924063, 0.967545082, 1.475254571, 1.530611618, 1.389143127]
#Hired labor to land ratio across specified bins - DATA 
HLph_data = [1.202307204, 1.042192279, 0.848612502, 0.783360938, 1.060585854, 1.048177806, 1.36807215]

#--------AVERAGES IN MODEL-------------------------------------------------
#Average output per worker (operators+hired) - MODEL
AggVApw = sum(y_model)/sum(n_model)
#Average hired labor per hectare (hired labor to land ratio) - MODEL
AggHLph = sum(h_model)/sum(l_model)

#-------CONDITIONAL DISTRIBUTIONS RELATIVE TO AVERAGES IN MODEL------------
#Conditional (on operating) distribution of OUTPUT PER WORKER(OPERATORS+HIRED) for specified bins relative to average - MODEL
VApw_model = (y_model./n_model)/AggVApw
#Conditional (on operating) distribution of HIRED LABOR PER HECTARE(HIRED LABOR TO LAND RATIO) for specified bins relative to average - MODEL
HLph_model = (h_model./l_model)/AggHLph

begin
    farm_pdf_85panel_matrix = reshape(farm_pdf_85panel, :, 1)
    data = hcat(reshape(l_model, :, 1), farm_pdf_85panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p1 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 0.4),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Farms")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "1985 Survey Data"),
             bar_width=width, color=:auto)
    end

    savefig(p1, "BE_Farm Size Distribution Across Specified Bins.png")
end

begin
    land_pdf_85panel_matrix = reshape(land_pdf_85panel, :, 1)
    data = hcat(reshape(l_model, :, 1), land_pdf_85panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p2 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 0.4),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Land")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "1985 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p2, "BE_Land Size Distribution Across Specified Bins.png")
end

begin
    VApw_model_col = reshape(VApw_model, :, 1)
    VApw_data_col = reshape(VApw_data, :, 1)
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
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "1985 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p3, "BE_Value Added per Worker Distribution Across Specified Bins.png")

end

begin
    HLph_model_col = reshape(HLph_model, :, 1)
    HLph_data_col = reshape(HLph_data, :, 1)
    data = hcat(HLph_model_col, HLph_data_col)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p4 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 2),
              xlabel="Farm Size Class in Ha", ylabel="Value added per worker")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "1985 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p4, "BE_Hired Labor per Hectare Distribution Across Specified Bins.png")

end



# Find active farm operators
indACTIVE = findall(x -> x == 1, Indic)


#Distribution of TFP
#--------------------------------------------------------------------------
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
#--------------------------------------------------------------------------
println("\nCORR of log-TFP AND log-TFPR")
CORR_lTFP_lTFPR_model = cor(lTFP_vec_ac, lTFPR_vec_ac)
println(CORR_lTFP_lTFPR_model)




#**************************************************************************
#OUTPUT / OUTPUT PER WORKER
#**************************************************************************

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

























# ReplicationPackage_KexinCHEN_JinxuMI.jl
# Replicating the article "Land Reform and Productivity: 
# A Quantitative Analysis with Micro Data" (Adamopoulos & Restuccia, 2020)

# LR_main.jl
# Government-mandated land reform.

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


#PARAMETERS OF THE LAND REFORM/REDISTRIBUTION
#ceiling
l_max       = 5
#degree of enforcement: a fraction THETA of potentially constrained farmers gets to keep their land
THETA       = 0.8
#fraction of landless that receive land
PSI0        = 0.329
#fraction of small-holders that receive land
PSI1        = 0.108


#FIND FARMERS WITH FARM SIZE BEFORE THE REDISTRIBUTION IN EXCESS OF 5 Ha
#construct indicator function that takes the value of 1 if farmer is potentially constrained and 0 otherwise
c_vec = zeros(N) 
II = findall(x -> x >= l_max, l_vec) 
c_vec[II] .= 1 

#with prob. THETA potentially constrained farmers are not hit by the ceiling and keep their original farm size
lhat_vec = (1 .- c_vec) .* l_vec .+ c_vec .* (THETA .* l_vec .+ (1 .- THETA) .* l_max)

#excess land per capita to be redistributed
extraLN     = LN - (sum(lhat_vec))/N

# A fraction PSI0 of landless (those with g < g_lb) receive land
min_index = findfirst(==(1), Indic)
print(min_index)

Indic = float(Indic)

for i = 1:(min_index - 1)
    Indic[i] = PSI0
end

print(PSI0)

# Share of landless after the reform
Nw_post = sum((1 .- Indic)) / N

# Find the smallholders that receive land
INsmall = findfirst(cdf_g .> hired_lab_sh * (1 + PSI1))
g_cutoff_small = g_vec[INsmall]
g_cutoff_Indic_small = findfirst(x -> x > g_cutoff_small, g_vec) - 1
sh_small = cdf_g[INsmall]

#how much land each landless and smallholder receives: total extra available land divided by the total number of recipients
lstar = extraLN / (sum(Indic[1:g_cutoff_Indic_small]) / N)

for ii = 1:(min_index - 1)
    lhat_vec[ii] = lstar
end

for iii = min_index:g_cutoff_Indic_small
    lhat_vec[iii] += lstar
end

#update l_vec and create new conditional distribution 
l_vec = lhat_vec
ConDist_LR = Indic / sum(Indic)  

#Post-reform average farm size
AFS_post = sum(ConDist_LR .* l_vec)
l_value = AFS_post


#Equilibrium
#--------------------------------------------------------------------------
#Once the land distribution from the government mandated land reform has been determined find the w for which the labor market clears.

oc_vec = nothing
of_vec = nothing
I = nothing

GC.gc()


#ADDITIONAL FACTORS
#A=1.1
#KAPPAf = KAPPAf*1.65
#KAPPAc = KAPPAc*1.65

#--------------------------------------------------------------------------
#Parameters
#--------------------------------------------------------------------------

# Parameters
GAMMA = 0.7
ALPHA = 0.3
A = 1 # Normalization
KAPPAf = 1 # Normalization
Pf = 1 # Normalization
Pc = 1 # Normalization

AFS = 3.7 # Average Farm Size in 1984-85 (Total Land Area/Total Number of Holdings) - from panel sample

LN = AFS * (1 - hired_lab_sh) # Choose land per capita to match AFS
KAPPAc = 1.25 # Chosen to match closely the cash-food productivity ratio


# Function LR_main_eval
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


guess  = [0.7]

# Additional parameters, assuming 'A' is needed as demonstrated in the previous explanation
params = [A] 

result = nlsolve((res, x) -> res .= LR_main_eval(x, params), guess, show_trace=true, xtol=1e-16)


x = result.zero # The solution
print(x)

# Assigning the result to 'w'
w = x


# Check for convergence
converged = result.f_converged
println("Converged: ", converged)


#Solve Farmer's Problem 
nf_vec = ((1 .- ALPHA) .* GAMMA .* (Pf ./ w) .* ((A .* KAPPAf .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA * (1 .- ALPHA)))
nc_vec = ((1 .- ALPHA) .* GAMMA .* (Pc ./ w) .* ((A .* KAPPAc .* g_vec) .^ (1 .- GAMMA)) .* (l_vec .^ (ALPHA .* GAMMA))) .^ (1 ./ (1 .- GAMMA * (1 .- ALPHA)))
yf_vec = (A .* KAPPAf .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nf_vec .^ (1 .- ALPHA)) .^ GAMMA
yc_vec = (A .* KAPPAc .* s_vec) .^ (1 .- GAMMA) .* (l_vec .^ ALPHA .* nc_vec .^ (1 .- ALPHA)) .^ GAMMA
PIf_vec = (1 .- GAMMA .* (1 .- ALPHA)) .* Pf .* yf_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pf .* Cf .* ones(length(N))
PIc_vec = (1 .- GAMMA .* (1 .- ALPHA)) .* Pc .* yc_vec .* (phi_vec .^ (1 .- GAMMA)) .- Pc .* Cc .* ones(length(N))


# Determine technology choice vector
oc_vec = zeros(Int, N) 

typeof(PIf_vec)
typeof(w)

max_value = max.(PIf_vec, w)

I = findall(PIc_vec .> max_value)

oc_vec[I] .= 1 


#Reform implied occupational choice vector
of_vec = zeros(Int, N)

for ij = 1:minimum(I) .- 1
        of_vec[ij] = 1
end

Indic = float(Indic)
of_vec = of_vec .* Indic


# Remaining hired workers (=landless)
Nww = sum((1 .- oc_vec .- of_vec)) / N
f = (sum(of_vec .* nf_vec) ./ N) .+ (sum(oc_vec .* nc_vec) ./ N) .- Nww

# Number of operators by crop
Nf = sum(of_vec) / N
Nc = sum(oc_vec) / N

# Output and hired labor vectors
y_vec = of_vec .* yf_vec + oc_vec .* yc_vec
n_vec = of_vec .* nf_vec + oc_vec .* nc_vec


# Bins for distribution
bin1 = 1
bin2 = 2
bin3 = 5
bin4 = 7
bin5 = 10
bin6 = 15


# Indices for bins
lhat1in = findall((l_vec .> 0) .& (l_vec .< bin1))
lhat2in = findall((l_vec .< bin2) .& (l_vec .>= bin1))
lhat3in = findall((l_vec .< bin3) .& (l_vec .>= bin2))
lhat4in = findall((l_vec .< bin4) .& (l_vec .>= bin3))
lhat5in = findall((l_vec .< bin5) .& (l_vec .>= bin4))
lhat6in = findall((l_vec .< bin6) .& (l_vec .>= bin5))
lhat7in = findall(l_vec .> bin6)


#Compute distributions: farms, land, output, hired labor
if isempty(lhat1in)==1
    farm1   = 0
    land1   = 0
    output1 = 0
    hirelab1=0
else
    farm1   = sum(ConDist_LR[lhat1in])
    land1   = sum(l_vec[lhat1in].*ConDist_LR[lhat1in])/l_value
    output1 = sum(y_vec[lhat1in].*ConDist_LR[lhat1in])
    hirelab1 = sum(n_vec[lhat1in].*ConDist_LR[lhat1in])
end

if isempty(lhat2in)
    farm2 = 0
    land2 = 0
    output2 = 0
    hirelab2 = 0
else
    farm2 = sum(ConDist_LR[lhat2in])
    land2 = sum(l_vec[lhat2in] .* ConDist_LR[lhat2in]) / l_value
    output2 = sum(y_vec[lhat2in] .* ConDist_LR[lhat2in])
    hirelab2 = sum(n_vec[lhat2in] .* ConDist_LR[lhat2in])
end

if isempty(lhat3in)
    farm3 = 0
    land3 = 0
    output3 = 0
    hirelab3 = 0
else
    farm3 = sum(ConDist_LR[lhat3in])
    land3 = sum(l_vec[lhat3in] .* ConDist_LR[lhat3in]) / l_value
    output3 = sum(y_vec[lhat3in] .* ConDist_LR[lhat3in])
    hirelab3 = sum(n_vec[lhat3in] .* ConDist_LR[lhat3in])
end

if isempty(lhat4in)
    farm4 = 0
    land4 = 0
    output4 = 0
    hirelab4 = 0
else
    farm4 = sum(ConDist_LR[lhat4in])
    land4 = sum(l_vec[lhat4in] .* ConDist_LR[lhat4in]) / l_value
    output4 = sum(y_vec[lhat4in] .* ConDist_LR[lhat4in])
    hirelab4 = sum(n_vec[lhat4in] .* ConDist_LR[lhat4in])
end

if isempty(lhat5in)
    farm5 = 0
    land5 = 0
    output5 = 0
    hirelab5 = 0
else
    farm5 = sum(ConDist_LR[lhat5in])
    land5 = sum(l_vec[lhat5in] .* ConDist_LR[lhat5in]) / l_value
    output5 = sum(y_vec[lhat5in] .* ConDist_LR[lhat5in])
    hirelab5 = sum(n_vec[lhat5in] .* ConDist_LR[lhat5in])
end

if isempty(lhat6in)
    farm6 = 0
    land6 = 0
    output6 = 0
    hirelab6 = 0
else
    farm6 = sum(ConDist_LR[lhat6in])
    land6 = sum(l_vec[lhat6in] .* ConDist_LR[lhat6in]) / l_value
    output6 = sum(y_vec[lhat6in] .* ConDist_LR[lhat6in])
    hirelab6 = sum(n_vec[lhat6in] .* ConDist_LR[lhat6in])
end

if isempty(lhat7in)
    farm7 = 0
    land7 = 0
    output7 = 0
    hirelab7 = 0
else
    farm7 = sum(ConDist_LR[lhat7in])
    land7 = sum(l_vec[lhat7in] .* ConDist_LR[lhat7in]) / l_value
    output7 = sum(y_vec[lhat7in] .* ConDist_LR[lhat7in])
    hirelab7 = sum(n_vec[lhat7in] .* ConDist_LR[lhat7in])
end


# Distributions of output, hired labor, farms, farmers in the model
l_model = [land1, land2, land3, land4, land5, land6, land7]
y_model = [output1, output2, output3, output4, output5, output6, output7]
h_model = [hirelab1, hirelab2, hirelab3, hirelab4, hirelab5, hirelab6, hirelab7]
f_model = [farm1, farm2, farm3, farm4, farm5, farm6, farm7]
n_model = h_model + f_model


# Benchmark prices for crops
Pf_star = 1
Pc_star = 1

# Aggregate Output per capita
VApw = (Pc_star * sum(oc_vec .* yc_vec) + Pf_star * sum(of_vec .* yf_vec)) / N
VApw_post = VApw

# Landless
landless_post = Nww

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


# For the farm model plot
begin
    farm_pdf_03panel_matrix = reshape(farm_pdf_03panel, :, 1)
    data = hcat(reshape(l_model, :, 1), farm_pdf_03panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p1 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 0.5),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Farms")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end

    savefig(p1, "LR_main_Farm Size Distribution Across Specified Bins.png")
end

begin
    land_pdf_03panel_matrix = reshape(land_pdf_03panel, :, 1)
    data = hcat(reshape(l_model, :, 1), land_pdf_03panel_matrix)
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xtick_values = (1:length(xtick_labels), xtick_labels)  # Renamed to avoid conflicts

    group_size = size(data, 2)
    width = 0.35  # Adjust width to fit bars side by side
    offsets = [-width/2, width/2]  # Position adjustments for each bar
    p2 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 0.5),
              xlabel="Farm Size Class in Ha", ylabel="Fraction of Land")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p2, "LR_main_Land Size Distribution Across Specified Bins.png")
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
    p3 = plot(size=(600, 400), legend=:topright, xticks=xtick_values, ylims=(0, 3),
              xlabel="Farm Size Class in Ha", ylabel="Value added per worker")

    # Loop to plot bars with adjusted positions for each group
    for i in 1:group_size
        bar!(xtick_values[1] .+ offsets[i], data[:, i], label=(i == 1 ? "Model" : "2003 Survey Data"),
             bar_width=width, color=:auto)
    end        
    savefig(p3, "LR_main_Value Added per Worker Distribution Across Specified Bins.png")

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
    savefig(p4, "LR_main_Hired Labor per Hectare Distribution Across Specified Bins.png")

end







# All labor (hired+operators) by crop
LAB_f = ((sum(of_vec .* nf_vec)) / N) + Nf
LAB_c = ((sum(oc_vec .* nc_vec)) / N) + Nc

# Crop Outputs
Yc = (sum(oc_vec .* yc_vec)) / N
Yf = (sum(of_vec .* yf_vec)) / N

# Labor productivity by crop
VApw_cash = Yc / LAB_c
VApw_food = Yf / LAB_f

# Total effects of land reform (Table 5 results)
println("")
println("% CHANGE IN AVERAGE FARM SIZE FOLLOWING THE REFORM")
AFS_EFF = 100 * (AFS_post - AFS_BE) / AFS_BE

println("")
println("% CHANGE IN VALUE ADDED PER WORKER FOLLOWING THE REFORM")
VApw_EFF = 100 * (VApw_post - VApw_BE) / VApw_BE

println("")
println("CHANGE IN % OF LANDLESS FOLLOWING THE REFORM")
landless_CHANGE = 100 * (landless_post - landless_BE)









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








end
