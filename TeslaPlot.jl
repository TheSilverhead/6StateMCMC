using DifferentialEquations, CSV, DataFrames, DelimitedFiles, StatsBase

function Model!(dy,y,par,t)
  #IFN, ODE 1 parameters
  k11=0 #PR8, RIGI is assumed antagonized
  n=3
  RIGI=1
  k12=par[1]
  k13=par[2]
  k14=par[3]
  #IFN_env, ODE 2 parameters
  k21=par[4]
  tau2=par[5]
  #STATP, ODE 3 parameters
  k31=par[6]
  k32=par[7]
  k33=par[8]
  #IRF7, ODE 4 parameters
  k41=par[9]
  k42=par[10]
  #IRF7P, ODE 5 parameters
  k51=par[11]
  #Cells, ODE 6 parameters
  k61=par[12]
  #Virus, ODE 7 parameters
  k71=par[13]
  k72=par[14]
  k73=par[15]

  #ODE System
  dy[1]=y[6]*(k11*RIGI*y[7]+(k12*y[7]^n)/(k13+y[7]^n)+k14*y[5])-k21*y[1]
  dy[2]=k21*y[1]-tau2*y[2]
  dy[3]=y[6]*(k31*y[2])/(k32+k33*y[2])-0.3*y[3]
  dy[4]=y[6]*(k41*y[3]+k42*y[5])-0.3*y[4]
  dy[5]=y[6]*k51*y[4]-0.3*y[5]
  dy[6]=-k61*y[6]*y[7]
  dy[7]=(k71*y[6]*y[7])/(1+k72*(y[2]*7E-5))-k73*y[7]
end

chainparams=CSV.read("chainparams.csv") #Read in all sample's parameters
stateNames = ["IFN","IFNe","STATP","IRF7","IRF7P","Live Cells","Virus"]
u0 = [0, 0, 0, 0.72205, 0, 1, 6.9e-8] #Shifted Initial Conditions
tspan = (0,48.0) #Time (start, end)
shift= [7.939, 0, 12.2075, 13.4271, 0, 0, 0] #Noise shift on data
alg = Vern7()   #ODE solver
logSpace=[1,3,4] #Which outputs should be shown in log2 space?
#Listen. This looks like trash, ok? But it takes a smooth time array, your desired time points, and makes a sorted, non-duplicated array
#tgraph = rle(sort(vcat(range(tspan[1],tspan[2],length=100),[0.25,0.5,1,1.5,2,4,6,8,12,18,24])))[1]
tgraph = [0.25,0.5,1,1.5,2,3,4,6,8,9,12,18,24]

function ProbMaker(model,u0,tspan,p)
  prob = ODEProblem(model,u0,tspan,p)
  return prob
end

start=(1) #Tesla Plot start
stop=1000000 #Tesla Plot stop

solMatrix = fill(0.0,(stop-start+1),7,length(tgraph))

for i=start:stop
  print("Processing Sample #",i)
  p=chainparams[i,:]
  sol = solve(ProbMaker(Model!,u0,tspan,p),alg)
  solMatrix[(i-start+1),:,:]=sol(tgraph)[:,:]
  #EXported PARAMeters, LONG/SHORT/BEST version
  CSV.write("exparamlong.csv",DataFrame(solMatrix[(i-start+1),:,:]),append=true)
end

#Call python. Pass solution, add chains to Tesla Plot
