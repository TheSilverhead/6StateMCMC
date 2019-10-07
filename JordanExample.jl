using DifferentialEquations,ParameterizedFunctions, DiffEqParamEstim #DiffEq
using RecursiveArrayTools, StatsBase,Distributions #Vector of Arrays and stats
using StatsPlots, CSV, DataFrames,Printf,Dierckx,ProgressMeter

###########################
#ToDo
#Figure out how to deal with parameter boundaries
#Include priors into MH MCMC?
###########################

#Import the MCMC solver
include("main.jl")
#Import error function
include("LossFunctionJordan.jl")

function Model!(dy,y,par,t)
  #IFN, ODE 1 parameters
  k11=0 #PR8, RIGI is assumed antagonized
  k12=par[2]
  n=3
  k13=par[3]
  k14=par[4]
  tau1=par[22]
  #IFN_env, ODE 2 parameters
  k21=par[5]
  tau2=par[23]
  #STAT, ODE 3 parameters
  r31=par[6]
  k31=par[7]
  tau3=par[24]
  #STATP, ODE 4 parameters
  k41=par[8]
  k42=par[9]
  tau4=par[25]
  #IRF7, ODE 5 parameters
  k51=par[10]
  k52=par[11]
  tau5=par[26]
  #IRF7P, ODE 6 parameters
  k61=par[12]
  tau6=par[27]
  #Target cells, ODE 7 paramters
  k71=par[13]
  #Eclipse infected cells, ODE 8 parameters
  k81=par[14]
  k82=par[15]
  #Productive infected cells, ODE 9 parameters
  k91=par[16]
  #Viral count, ODE 10 parameters
  k10_1=par[17]
  k10_2=par[18]
  k10_3=par[19]
  #TJ Constants
  #TJ describes the binding of IFN and SOCS feedback
  TJtot=0.0001
  k11_1=par[20]
  k11_2=par[21]
  TJ=TJtot*(y[2]/(k11_1+y[2])*(1/(1+k11_2))); #Eq. 11

  #ODE System
  v=y[10]/(10^4) #scale virus effects
  i1=y[8]/(10^4) #scale eclipse cell effects
  dy[1]=(k11*v)+(k12*(v^n))/(k13+(v^n))+k14*y[6]-y[1]*tau1 #IFN in cytoplasm
  dy[2]=(k21*y[1])-(y[2]*tau2) #IFN in environment
  dy[3]=r31+k31*y[4]-y[3]*tau3 #STAT in cytoplasm
  dy[4]=(k41*TJ*y[3])/(k42+y[3])-y[4]*tau4 #STATP
  dy[5]=k51*y[4]+k52*y[6]-y[5]*tau5 #IRF7
  dy[6]=k61*y[5]-y[6]*tau6 #IRF7P
  dy[7]=-k71*y[7]*y[10] #Uninfected target cells
  dy[8]=k71*y[7]*y[10]-(k81*i1)/(1+k82*y[2]) #Eclipse infected cells
  dy[9]=(k81*i1)/(1+k82*y[2])-k91*y[9] #Productive infected cells
  dy[10]=(k10_1*y[9])/(1+k10_2*y[2])-k10_3*y[10] #Virus count

end

#parNames = ["a","b","c"]
#parNum = length(parNames)
stateNames = ["IFN","IFNenv","STAT1","STATP2n","IRF7","IRF7Pn","Target","Eclipse","Productive","Virus"]
parNum = 27

#Define information for ODE model
p=rand(parNum) #Parameter values
u0 = [7.94, 0, 262.3, 12.2, 14.15, 0, 0, 250000, 0, 7.5E-2] #Initial Conditions
tspan = (0.25,24.0) #Time (start, end)

#Contruct the ODE Problem
prob = ODEProblem(Model!,u0,tspan,p)
alg = Vern7()  #ODE solver

#Solve for the "True" solution
sol = solve(prob,alg)

## Generate data

data = CSV.read("C:/Users/Portable/Documents/GitHub/6StateMCMC/PR8.csv")
t = convert(Array,data[:Time])
mock = CSV.read("C:/Users/Portable/Documents/GitHub/6StateMCMC/Control.csv")
ci = CSV.read("C:/Users/Portable/Documents/GitHub/6StateMCMC/PR8CI95.csv")
titer = CSV.read("C:/Users/Portable/Documents/GitHub/6StateMCMC/ViralTiters.csv")
t_titer=convert(Array,titer[:Time])

removeCol = [:Time :IFN_env] #Strip time and IFN_env from data
for col in removeCol
  deletecols!(data,col)
  deletecols!(mock,col)
end

removeCol = [:Time :dNS1PR8] #Strip time from viral titers
for col in removeCol
  deletecols!(titer,col)
end

data = convert(Matrix,data)
mock = convert(Matrix,mock)
titer = convert(Matrix,titer)
ci = convert(Matrix,titer)

#Data points being measured
measured = [1,3,4,5]

## Supply prior distributions
priors = fill(Uniform(0,1), parNum)
parBounds = fill([0.0, Inf], parNum) #Fill bounds with 0/inf
parBounds[14]=Float64[0.0, 1.0]; #Modify bounds with known values
parBounds[16]=Float64[0.0, 1.0];
parBounds[19]=Float64[0.0, 1.0];

lossFunc = LossLog(t,t_titer,data,measured,mock,ci,titer)

sampleNum = Int(1e5)
result = ptMCMC(prob,alg,priors,parBounds,lossFunc,sampleNum)

bestPars= dropdims(permutedims(result[1], [1, 3, 2])[argmax(result[2],dims=1),:],dims=1)
bestChains = [remake(prob;p=bestPars[i,:]) for i=1:size(bestPars,1)]
chainSols = [solve(problem,alg) for problem in bestChains]

###############
cd("C:/Users/Portable/Documents/GitHub/6StateMCMC")
ODEplots = Vector(undef,length(u0))

mockSplines = [Spline1D(t,mock[:,i];k=1) for i=1:size(mock,2)]

plotRecipe = (t,x,i,j) -> plot(t,x,label="Chain $j",title=stateNames[i],framestyle=:box,yformatter=y->"$(@sprintf("%.0e",y))",legend=false)

Lfc(c,m) = @. log2(c+1) - log2(m+1)
#loop through states
for i in eachindex(u0)

#Make the first plot for the first chain
  if i ∈ measured
      midx = findfirst(x->x==i,measured)
      logScaledSim = Lfc(chainSols[1][i,:], mockSplines[midx](chainSols[1].t))
      ODEplots[i] = plotRecipe(chainSols[1].t,logScaledSim,i,1)
      scatter!(t,data[:,findfirst(x->x==i,measured)],label="Data: "*stateNames[i],layout=(2,3))
  else
      ODEplots[i] = plotRecipe(chainSols[1].t,chainSols[1][i,:],i,1)
  end

#Loop through and add the rest of the chains
  for j=2:length(chainSols)
    if i ∈ measured
      midx = findfirst(x->x==i,measured)
      logScaledSim = Lfc(chainSols[j][i,:], mockSplines[midx](chainSols[j].t))
      plot!(chainSols[j].t,logScaledSim,label="Chain $j")
    else
      plot!(chainSols[j].t,chainSols[j][i,:],label="Chain $j")
    end
  end

end


datafitPlot = plot(ODEplots...,layout=(2,5),size=(1000,600))

savefig(datafitPlot,"JordanDataFit.pdf")


using AverageShiftedHistograms

#Posterior
HistData = [plot(ash(result[1][1:end,i,1]),title="Par $i") for i=1:parNum]
histPlot=plot(HistData...,legend=false,size=(2000,1200))

savefig(histPlot,"histPlot.pdf")

logPostPlot = plot(result[2],legend =false,title="logPosterior")
savefig(logPostPlot,"logPostPlot.pdf")

#acceptance Rate
AcceptPlot = plot(result[3],legend =false,title="Acceptance Rate")
savefig(AcceptPlot,"AcceptPlot.pdf")
