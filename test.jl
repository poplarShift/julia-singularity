#!/usr/bin/env julia
print("using ...")

###I didnt check manually one by one but I'm pretty sure I use all of these
using DifferentialEquations, RecursiveArrayTools, DiffEqParamEstim
using Optimization, OptimizationBBO ,OptimizationOptimJL, ForwardDiff
using Plots, CSV, DataFrames, Random , DataFramesMeta
using Zygote, OptimizationPolyalgorithms, SciMLSensitivity, ModelingToolkit, LinearAlgebra
using StructuralIdentifiability, Interpolations, ComponentArrays

print("... all imports done.\n")
###MANUAL SETUP FOR EACH LAKE-RUN, COULDNT FIGURE HOW TO LOOP THROUGH LAKES OR RUNS (EACH RUN STARTS WITH NEW INITIAL PARAM VALUES)
####Working lake by lake, 

lake = ARGS[1] # e.g. "Horntjernet". Should also permit special characters
run = ARGS[2] # can be running number, but ideally unique identifier (e.g. UUID4)

#Import data and reshape

print("Reading data\n")
dataset = CSV.read("Output/test/data/SmallDataset.csv", DataFrame, header=[1,], drop=[1,], missingstring="NA")
select!(dataset, Not(1))

LakeInfo = CSV.read("Output/test/data/coords.csv", DataFrame, header=[1,], normalizenames=true, missingstring="NA")
LakeInfo = dropmissing(unique(LakeInfo[:, [:Lake_name, :Latitude, :Longitude]]))

JulyTemp = CSV.read("Output/test/data/July_temperature_CHELSA_for_Marieke.csv", DataFrame, header=[1,])
JulyTemp = dropmissing(unique(JulyTemp[:, [:Longitude, :Latitude, :jult, :bp]]))
JulyTemp_id = dropmissing(leftjoin(LakeInfo, JulyTemp, on=[:Latitude, :Longitude]))
JulyTemp_id.bp = JulyTemp_id.bp * 1000

#taxa=["Trees_Shrubs", "Aquatics", "Bryophytes", "Cryptogams", "Andromeda.polifolia","Dryas.octopetala","Empetrum.nigrum","Vaccinium.uliginosum","Vaccinium.vitis.idaea","Vaccinium.myrtillus","Avenella.flexuosa", "Oreojuncus.trifidus","Linnaea.borealis", "Saxifraga.oppositifolia"]
#playset = dataset[∈(["Andromeda.polifolia","Dryas.octopetala","Empetrum.nigrum","Vaccinium.uliginosum","Vaccinium.vitis.idaea","Vaccinium.myrtillus","Avenella.flexuosa","Oreojuncus.trifidus","Linnaea.borealis", "Saxifraga.oppositifolia"]).(dataset.Taxa_sedaDNA), :]

planttaxa=["Empetrum.nigrum",
    "Vaccinium.vitis.idaea","Vaccinium.myrtillus"]


playset = dataset[∈(planttaxa).(dataset.Taxa_sedaDNA), :]


###Plants
gdf_lake = unique!(sort!(playset[playset[:,"Lake_name"] .== lake, :], [order(:Taxa_sedaDNA), order(:MedCalAge, rev=true)]))
gdf = groupby(gdf_lake, :Taxa_sedaDNA)
gdf_wide_lake=unstack(gdf_lake, :MedCalAge, :Repeats_Prop_Weighted)
gdf_wide_lake=dropmissing(gdf_wide_lake)
gdf_wide_lake=dropmissing(gdf_wide_lake[sum(eachcol(gdf_wide_lake[!,4:end])).>0,:]) ####NEW
gdf_wide_lake = @rorderby gdf_wide_lake findfirst(==(:Taxa_sedaDNA), planttaxa)
data_matrix=Matrix(select(gdf_wide_lake, Not(1:3)))

###Temperature forcing data
CHELSA_T = sort!(JulyTemp_id[JulyTemp_id[:,"Lake_name"] .== lake, :], order(:bp))
#normalize time for the CHELSA data
CHELSA_T.bp = -CHELSA_T.bp/1000 .+ first(sort(playset.MedCalAge, rev=true))/1000
TempForce = sort!(CHELSA_T, :bp)


n = length(gdf_wide_lake.Taxa_sedaDNA)
taxa = Vector(gdf_wide_lake.Taxa_sedaDNA)


# Adapting time...using real very long timespan obvs causes problems, will aim to look at problem at reduced scale, /1000 seems a good startpoint
#Using the oldest sample date as T0...

# Making sure that each taxa time step are equal

@assert all(all(sort(gdf[1].MedCalAge) .== sort(df.MedCalAge)) for df in gdf[2:end]) #throws error if missing timestep for give species 

tsteps = .- sort(gdf[1].MedCalAge, rev=true)./1000 .+ first(sort(playset.MedCalAge, rev=true))/1000



# Starting parameters...

print("Starting parameters\n")

ρ0= rand(Float64, n)*5
#α0 =rand(Float64, n, n)
i0 = rand(Float64, n)*5
c0 = rand(Float64, n)*5
a0 = rand(Float64, n)*5

tspan = (Float64(first(tsteps)),Float64(last(tsteps)))

###start values averaging the first 3 observations 

u0 = zeros(n)

for i in 1:n
    u0[i]= data_matrix[i,1] + data_matrix[i,2] + data_matrix[i,3] == 0.0 ? 0.001 : (data_matrix[i,1] + data_matrix[i,2] + data_matrix[i,3])/3
end

###save the start parameter values
CSV.write(string("Output/test/",lake,"_",run,"_startparam.csv"), DataFrame(taxa=taxa, ρ0=ρ0, i0=i0, c0=c0, a0=a0))

print("Start models\n")
###Logistic model


@variables t, N(..)[1:n]

@parameters begin
    ρ[1:n], [bounds=(0, 5), tunable=true] #(..)
end

K=1

D = Differential(t) 

vec_eqs = D.(N(t)) .~ (1 .+ I*-N(t)./K) .* N(t) .* ρ

@named Log_mod = ODESystem(vec_eqs)

Log_prob = ODEProblem(structural_simplify(Log_mod), u0, tspan, ρ0)

Log_mod

###Logistic model Temp dependant

#Temperature forcing variable
interp_linear = LinearInterpolation(TempForce.bp, TempForce.jult)

@variables t, N(..)[1:n]
@parameters begin
    ρ[1:n], [bounds=(0, 5), tunable=true] #(..)
end


T(t) = interp_linear(t)
@register_symbolic T(t)

K=1
#@derivatives D'~t
D = Differential(t)



vec_eqs = D.(N(t)) .~ (1 .- I*N(t)./K) .* N(t) .* ρ .*T(t)

@named LogT_mod = ODESystem(vec_eqs)

LogT_prob = ODEProblem(structural_simplify(LogT_mod), u0, tspan, ρ0)

LogT_mod

####Competition model Pellisier et al. 2017 (no herbivory)


m=ones(n,n)
m[diagind(m)] .= 0.0

#Temperature forcing variable
interp_linear = LinearInterpolation(TempForce.bp, TempForce.jult)

@variables t, N(..)[1:n]

@parameters begin
    ρ[1:n], [bounds=(0, 5), tunable=true]
    i[1:n], [bounds=(0, 5), tunable=true]
    c[1:n], [bounds=(0, 5), tunable=true]
end

D = Differential(t)

T(t) = interp_linear(t)
@register_symbolic T(t)

vec_eqs = D.(N(t)) .~ N(t) .* (ρ .- c .*N(t)  .-  (i .* m *N(t)))


#(1 .-N(t) .-m .*α *N(t)).* ρ .*N(t)*T(t)

@named C_mod = ODESystem(vec_eqs)


C_prob = ODEProblem(structural_simplify(C_mod), u0, tspan, (ρ0, i0, c0))
C_mod

####SPECIFY WHICH MODEL WE ARE RUNNING, THE REST OF THE CODE SHOULD BE ADAPTABLE######
######################################################################################
Data = data_matrix

LBlog =repeat([0.0001], n)
UBlog = repeat([5.0], n)

LBc =repeat([0.0001], n*3)
UBc = repeat([5.0], n*3)

LBch =repeat([0.0001], n*4)
UBch = repeat([5.0], n*4)

###BBO Logistic

cost_function = build_loss_objective(Log_prob, Tsit5(), L2Loss(tsteps, Data, differ_weight=0.5,data_weight=0.5),Optimization.AutoForwardDiff(),maxiters=10_000,verbose=true)
optprob = Optimization.OptimizationProblem(cost_function, u0, lb=LBlog, ub=UBlog) 
result_bbo_log = solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=10_000) # callback=callback,
pred_bbo_log = Array(solve(ODEProblem(structural_simplify(Log_mod), u0, tspan, result_bbo_log.u), Tsit5(),saveat = tsteps))

###BBOBFGS Logistic

function loss(p)
    sol = solve(Log_prob, Tsit5(), p = p, saveat = tsteps)
    loss = sum(abs2, sol .- Data)
    return loss, sol
end

callback = function (p, l, sol)
    display(l)
    plt = plot(sol, ylim = (0, 1))
    display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, result_bbo_log.u , lb=LBlog, ub=UBlog)
result_bbobfgs_log = Optimization.solve(optprob, BFGS(initial_stepnorm=0.01))
pred_bbobfgs_log = Array(solve(ODEProblem(structural_simplify(Log_mod), u0, tspan, result_bbobfgs_log.u), Tsit5(),saveat = tsteps))

###BBO LogisticTemp

cost_function = build_loss_objective(LogT_prob, Tsit5(), L2Loss(tsteps, Data, differ_weight=0.5,data_weight=0.5),Optimization.AutoForwardDiff(),maxiters=10_000,verbose=true)
optprob = Optimization.OptimizationProblem(cost_function, u0, lb=LBlog, ub=UBlog) 
result_bbo_logt = solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=10_000) # callback=callback,
pred_bbo_logt = Array(solve(ODEProblem(structural_simplify(LogT_mod), u0, tspan, result_bbo_logt.u), Tsit5(),saveat = tsteps))

###BBOBFGS LogisticTemp

function loss(p)
    sol = solve(LogT_prob, Tsit5(), p = p, saveat = tsteps)
    loss = sum(abs2, sol .- Data)
    return loss, sol
end

callback = function (p, l, sol)
    display(l)
    plt = plot(sol, ylim = (0, 1))
    display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, result_bbo_logt.u , lb=LBlog, ub=UBlog)
result_bbobfgs_logt = Optimization.solve(optprob, BFGS(initial_stepnorm=0.01), reltol = 1e-8)
pred_bbobfgs_logt = Array(solve(ODEProblem(structural_simplify(LogT_mod), u0, tspan, result_bbobfgs_logt.u), Tsit5(),saveat = tsteps))

###BBO C

cost_function = build_loss_objective(C_prob, Tsit5(), L2Loss(tsteps, Data, differ_weight=0.5,data_weight=0.5),Optimization.AutoForwardDiff(),maxiters=10_000,verbose=true)
optprob = Optimization.OptimizationProblem(cost_function, u0, lb=LBc, ub=UBc) 
result_bbo_c = solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=10_000) # callback=callback,
pred_bbo_c = Array(solve(ODEProblem(structural_simplify(C_mod), u0, tspan, result_bbo_c.u), Tsit5(),saveat = tsteps))

###BBOBFGS C
function loss(p)
    sol = solve(C_prob, Tsit5(), p = p, saveat = tsteps)
    if size(sol) == size(Data)
        return sum(abs2, sol .- Data)
    else
        return Inf
    end
end


adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, result_bbo_c.u , lb=vec(LBc), ub=vec(UBc))
result_bbobfgs_c = Optimization.solve(optprob, BFGS(),reltol = 1e-8)
pred_bbobfgs_c = Array(solve(ODEProblem(structural_simplify(C_mod), u0, tspan, result_bbobfgs_c.u), Tsit5(),saveat = tsteps))

###Retain best optimisation method for each model based on SSError...likely very small difference


pred_log_name = sum(vec(sum((Data-pred_bbo_log).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_log).^2, dims=2))) ? "pred_bbo_log" : "pred_bbobfgs_log"
pred_logt_name = sum(vec(sum((Data-pred_bbo_logt).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_logt).^2, dims=2))) ? "pred_bbo_logt" : "pred_bbobfgs_logt"
pred_c_name = sum(vec(sum((Data-pred_bbo_c).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_c).^2, dims=2))) ? "pred_bbo_c" : "pred_bbobfgs_c"
    



pred_log = sum(vec(sum((Data-pred_bbo_log).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_log).^2, dims=2))) ? pred_bbo_log : pred_bbobfgs_log
pred_logt = sum(vec(sum((Data-pred_bbo_logt).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_logt).^2, dims=2))) ? pred_bbo_logt : pred_bbobfgs_logt
pred_c = sum(vec(sum((Data-pred_bbo_c).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_c).^2, dims=2))) ? pred_bbo_c : pred_bbobfgs_c


    
result_log = sum(vec(sum((Data-pred_bbo_log).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_log).^2, dims=2))) ? result_bbo_log : result_bbobfgs_log
result_logt = sum(vec(sum((Data-pred_bbo_logt).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_logt).^2, dims=2))) ? result_bbo_logt : result_bbobfgs_logt
result_c = sum(vec(sum((Data-pred_bbo_c).^2, dims=2))) < sum(vec(sum((Data-pred_bbobfgs_c).^2, dims=2))) ? result_bbo_c : result_bbobfgs_c
    
   
modelname=[pred_log_name ,pred_logt_name, pred_c_name]
names=pushfirst!(modelname, "model")


###Retain best optim method for each model , compare SSERROR

df=DataFrame(
    a=taxa, 
    b=vec(sum((Data-pred_log).^2, dims=2)), 
    c=vec(sum((Data-pred_logt).^2, dims=2)),
    d=vec(sum((Data-pred_c).^2, dims=2)))
    

df=permutedims(rename!(df, names),1)
df.mean = sum(eachcol(df[!,2:(n+1)]))./size(df[!,2:(n+1)])[2] #./size(tsteps)[1] #divide by tsteps to correct for sites with more or less samples 
df.lake.=lake
df=sort(df, :mean)

CSV.write(string("Output/test/",lake,"_",run,"_", "modelerror.csv"), df)
    

#plot bestmod
Predplotname = 
df[1,1] == "pred_bbobfgs_log" || df[1,1]== "pred_bbo_log" ? "Logistic model" : 
df[1,1] == "pred_bbobfgs_logt" || df[1,1]== "pred_bbo_logt" ? "Logistic model - temp dep" : 
df[1,1] == "pred_bbobfgs_c" || df[1,1]== "pred_bbo_c" ? "Simplified competition model" : "help"

Pred =
df[1,1] == "pred_bbobfgs_log" || df[1,1]== "pred_bbo_log" ? pred_log : 
df[1,1] == "pred_bbobfgs_logt" || df[1,1]== "pred_bbo_logt" ? pred_logt : 
df[1,1] == "pred_bbobfgs_c" || df[1,1]== "pred_bbo_c" ? pred_c : "help"

CSV.write(string("Output/test/",lake, "_", run,"_", Predplotname,".csv"), stack(DataFrame(hcat(taxa[1:n], Pred), ["taxa"; string.(gdf[1].MedCalAge)]), 2:size(data_matrix)[2]+1))

###Error plot for validation

plt = plot(title=lake)
scatter!(plt, tsteps, vec(sum((Data-pred_log).^2, dims=1)), label = pred_log_name, markercolor=1)
scatter!(plt, tsteps.+0.1, vec(sum((Data-pred_logt).^2, dims=1)), label = pred_logt_name, markercolor=2)
scatter!(plt, tsteps.+0.2, vec(sum((Data-pred_c).^2, dims=1)), label = pred_c_name,markercolor=3)  

savefig(plt, string("Output/test/",lake, "_", run, "_", "_SSErrorFig.png"))

###Parameters - Growth

series = [result_log.u, result_logt.u, result_c.u[1:n]]
rows = [[1:length(s);] for s in series]
df = DataFrames.flatten(DataFrame(g=[pred_log_name, pred_logt_name, pred_c_name], s=series, r=rows), [:s, :r])
df=unstack(df, :g, :s)
df.lake.=lake
df[!, "taxa"] .=  taxa[1:n]

CSV.write(string("Output/test/",lake, "_", run, "_", "r.csv"), df)

#Parameters - interactions

    
CSV.write(string("Output/test/",lake,"_",run,"_", "c_CompTrait.csv"), Tables.table(hcat(Array(taxa[1:n]), result_c.u[n+1: 2n])))
CSV.write(string("Output/test/",lake,"_",run,"_", "c_CompRate.csv"), Tables.table(hcat(Array(taxa[1:n]), result_c.u[2n+1: 3n])))
