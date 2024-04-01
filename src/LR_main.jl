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
of_vec .= zeros(Int, N)

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
    p1 = plot(bar([f_model, farm_pdf_03panel], label=["Model" "2003 Survey Data"], xticks=(1:7, xticklabels), legend=:topright), 
        colormap=:summer, xlabel="Farm Size Class in Ha", ylabel="Fraction of Farms", size=(600,400))
    xlims!(0, 8)
    ylims!(0, 0.4)
    plot!(p1, title="Farm Size Distribution")
    savefig(p1,"Farm Size Distribution.png")
end

begin
    # For the land model plot
    p2 = plot(bar([l_model, land_pdf_03panel], label=["Model" "2003 Survey Data"], xticks=(1:7, xticklabels), legend=:topright), 
        colormap=:summer, xlabel="Farm Size Class in Ha", ylabel="Fraction of Land", size=(600,400))
    plot!(p2, title="Land Distribution")
    savefig(p2,"Land Distribution.png")
end

begin
    # For the value added per worker model plot
    p3 = plot(bar([VApw_model, VApw_2003data], label=["Model" "2003 Survey Data"], xticks=(1:7, xticklabels), legend=:topright),
        colormap=:summer, xlabel="Farm Size Class in Ha", ylabel="Value Added Per Worker (Relative to Average)", size=(600,400))
    plot!(p3, title="Value Added Per Worker")
    savefig(p3,"Value Added Per Worker.png")
end

begin
    # For the hired labor per hectare model plot
    p4 = plot(bar([HLph_model, HLph_data_03], label=["Model" "2003 Survey Data"], xticks=(1:7, xticklabels), legend=:topright),
        colormap=:summer, xlabel="Farm Size Class in Ha", ylabel="Hired Labor Per Hectare (Relative to Average)", size=(600,400))
    plot!(p4, title="Hired Labor Per Hectare")
    savefig(p4,"Hired Labor Per Hectare.png")
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
