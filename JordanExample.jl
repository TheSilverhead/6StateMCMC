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


#Misc Constants
const Na=6.023E23
const C=5e5
const γ=10E9/Na
const d=0.0003

function Model!(dy,y,p,t)

  r1=p[1] #IFN induction
  k1=p[2]
  tau1=p[25] #IFN mRNA degradation
  vmax=γ*C*p[3] #IFN diffusion rate

  #STAT1 Constants
  r2=p[4]
  k3=p[5]
  tau2=p[26] #mRNA degradation

  #STATP2n Constants
  STAT=0.1 #STAT level, constant
  k4=p[6]
  K4=p[7]
  tau3=p[27] #STATp2n dephosphorylation

  #IRF7 Constants
  k51=p[15]
  k52=p[8]
  tau4=p[28] #mRNA degradation
  k6=p[16] #translation

  #NS1
  #NS1 is the primary antoginistic protein of influenza
  r3=p[9]
  n3=p[10]
  bm=p[11]
  NS=(r3*(t^n3)/((bm^n3)+(t^n3))) #NS1 protein kinetics

  #IC1/IC2
  #IC1/IC2 represent IFNB induction and nuclear mRNA antagonism,respectively
  sp=p[17] #[IC1]
  n1=p[18] #[IC1]
  sv=p[19] #[IC2]
  n2=p[20] #[IC2]
  d1=p[21] #[IC1]
  d2=p[22] #[IC2]
  IC1=(1+sp*((NS/d1)^n1))/(1+((NS/d1)^n1))
  IC2=((1+sv*((NS/d2)^n2))/(1+(NS/d2)^n2))

  #TJ Constants
  #TJ describes the binding of IFN a/b and SOCS feedback
  K2=p[13]
  TJtot=p[12]
  KTJ1=p[23]
  KTJ2=p[24]
  TJ=TJtot*((y[2])/(KTJ1+y[2]))*(1/(1+(KTJ2/d)))

  #ODE System
  dy[1]=r1+k1*y[6]*IC2-y[1]*(log(2)/tau1) #IFN in cytoplasm
  dy[2]=(vmax*y[1])/(K2+y[1])-y[2]*(log(2)/tau1) #IFN in environment
  dy[3]=(r2*IC1+k3*y[4])*IC2-y[3]*(log(2)/tau2) #STAT1
  dy[4]=k4*TJ*STAT/(2*(K4+STAT))-y[4]*log(2)/tau3 #STATP2n
  dy[5]=(k51*y[4]+k52*y[6])*IC2-y[5]*log(2)/tau4 #IRF7
  dy[6]=k6*y[5]*IC1-y[6]*log(2)/tau4 #IRF7Pn
end

#parNames = ["a","b","c"]
#parNum = length(parNames)
stateNames = ["IFN","IFNenv","STAT1","STATP2n","IRF7","IRF7Pn"]
parNum = 28

#Define information for ODE model
p=rand(parNum) #Parameter values
u0 = [7.9, 0, 262.0, 0.1, 14.0, 0.0] #Initial Conditions
tspan = (0.25,18.0) #Time (start, end)

#Contruct the ODE Problem
prob = ODEProblem(Model!,u0,tspan,p)
alg = Vern7()  #ODE solver

#Solve for the "True" solution
sol = solve(prob,alg)

## Generate data

data = CSV.read("PR8.csv")
t = convert(Array,data[:Time])
mock = CSV.read("Control.csv")

removeCol = [:Time :TLR7 :DDX58]
for col in removeCol
  deletecols!(data,col)
  deletecols!(mock,col)
end

data = convert(Matrix,data)
mock = convert(Matrix,mock)

#Data points being measured
measured = [1,3,5]

## Supply prior distributions
priors = fill(Uniform(0,1), parNum)
parBounds = fill([0.0, Inf], parNum)

lossFunc = LossLog(t,data,measured,mock)

sampleNum = Int(1e6)
result = ptMCMC(prob,alg,priors,parBounds,lossFunc,sampleNum)



bestPars= dropdims(permutedims(result[1], [1, 3, 2])[argmax(result[2],dims=1),:],dims=1)
bestChains = [remake(prob;p=bestPars[i,:]) for i=1:size(bestPars,1)]
chainSols = [solve(problem,alg) for problem in bestChains]

###############
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


datafitPlot = plot(ODEplots...,layout=(2,3),size=(1000,600))

savefig(datafitPlot,"JordanDataFit.pdf")


using AverageShiftedHistograms

#Posterior
HistData = [plot(ash(result[1][100_000:end,i,1]),title="Par $i") for i=1:parNum]
histPlot=plot(HistData...,legend=false,size=(2000,1200))

savefig(histPlot,"histPlot.pdf")

logPostPlot = plot(result[2],legend =false,title="logPosterior")
savefig(logPostPlot,"logPostPlot.pdf")

#acceptance Rate
AcceptPlot = plot(result[3],legend =false,title="Acceptance Rate")
savefig(AcceptPlot,"AcceptPlot.pdf")
