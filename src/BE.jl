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
    xticks = (1:length(xtick_labels), xtick_labels)

    p1 = bar(data, label=["Model" "1985 Survey Data"], legend=:topright, xticks=xticks, 
        ylims=(0, 0.4), xlabel="Farm Size Class in Ha", ylabel="Fraction of Farms", 
        size=(600, 400), bar_width=0.7)
    savefig(p1,"Farm Size Distribution Across Specified Bins.png")
end


begin
    land_pdf_85panel_matrix = reshape(land_pdf_85panel, :, 1)
    data = hcat(reshape(l_model, :, 1), land_pdf_85panel_matrix)
    # Define the bin labels
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xticks = (1:length(xtick_labels), xtick_labels)

    # Create the grouped bar plot
    p2 = bar(data, label=["Model" "1985 Survey Data"], legend=:topright, xticks=xticks,
        ylims=(0, 0.35), xlabel="Farm Size Class in Ha", ylabel="Fraction of Land",
        size=(600, 400), bar_width=0.7)
    savefig(p2,"Land Input Distribution Across Specified Bins.png")
end


begin
    VApw_model_col = reshape(VApw_model, :, 1)  # Making sure it's a column vector
    VApw_data_col = reshape(VApw_data, :, 1)  # Same for VApw_data
    
    data = hcat(VApw_model_col, VApw_data_col)  # Concatenate the data side by side
    
    # Define bin labels and corresponding x-tick positions
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xticks = (1:length(xtick_labels), xtick_labels)
    
    # Create the grouped bar plot
    p3 = bar(data, label=["Model" "1985 Survey Data"], legend=:topright, xticks=xticks,
        ylims=(0, 1.7), xlabel="Farm Size Class in Ha", ylabel="Value Added Per Worker (Relative to Average)",
        size=(600, 400), bar_width=0.7)    
    savefig(p3, "Value Added per Worker Distribution Across Specified Bins.png")
end

begin
    HLph_model_col = reshape(HLph_model, :, 1)  # Reshape if necessary
    HLph_data_col = reshape(HLph_data, :, 1)    # Reshape if necessary

    data = hcat(HLph_model_col, HLph_data_col)  # Concatenate data side by side

    # Define bin labels and corresponding x-tick positions
    xtick_labels = ["<1", "1-2", "2-5", "5-7", "7-10", "10-15", "15+"]
    xticks = (1:length(xtick_labels), xtick_labels)

    # Create the grouped bar plot
    p4 = bar(data, label=["Model" "1985 Survey Data"], legend=:topright, xticks=xticks,
    ylims=(0, 1.5), xlabel="Farm Size Class in Ha", ylabel="Hired Labor Per Hectare (Relative to Average)",
    size=(600, 400), bar_width=0.7)
    savefig(p4, "Hired Labor per Hectare Distribution Across Specified Bins.png")
end


# Find active farm operators
indACTIVE = findall(x -> x == 1, Indic)

#Distribution of TFP
#--------------------------------------------------------------------------
# TFP of Active Units
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