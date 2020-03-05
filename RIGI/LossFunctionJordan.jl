struct LossLog{time,time_titer,data,M,control,titer,shift} <: DiffEqBase.DECostFunction
  t::time #Time points for data
  t_titer::time_titer #Time points for virus count
  data::data # LFC of data
  measured::M #Indices in model with data
  mock::control #raw data for control
  vtiter::titer #Viral titer estimates
  shift::shift #Shift in data from noise
end

mutable struct NN
  sizeLayers
  numLayers
  biases
  weights
end


function (f::LossLog)(sol::DESolution)

t = f.t
data = f.data
measured = f.measured
mock = f.mock
vtiter = f.vtiter
shift=f.shift

  if sol.retcode != :Success
      return Inf
  end

  subSol = sol(t)

  sumsq = 0.0
  maxlog2(x) = log2(max(0,x)+1)
  @inbounds for i = 1:length(subSol) #Loop through each data time point
    for (j,m) in enumerate(measured) #loop though states
      sumsq +=(maxlog2(subSol[m,i]+shift[m]) - maxlog2(mock[i,j]) - data[i,j])^2
    end
  end

subSol = sol(t_titer)

@inbounds for i = 1:length(subSol) #Loop through each virus time point
  sumsq += 1000000*(titer[i,1]-subSol[7,i])^2 #Virus ssqe, weighted
end

return sumsq
end
