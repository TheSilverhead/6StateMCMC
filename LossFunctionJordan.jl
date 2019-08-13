struct LossLog{T,D,M,Mo} <: DiffEqBase.DECostFunction
  t::T #Time points for data
  data::D # LFC of data
  measured::M #Indices in model with data
  mock::Mo #raw data for control
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

  if sol.retcode != :Success
      return Inf
  end

  subSol = sol(t)

  sumsq = 0.0
  maxlog2(x) = log2(max(0,x)+1)
  @inbounds for i = 1:length(subSol) #Loop through each time point
    for (j,m) in enumerate(measured) #loop though states
      sumsq +=(maxlog2(subSol[m,i]) - maxlog2(mock[i,j]) - data[i,j])^2
    end
  end

return sumsq
end
