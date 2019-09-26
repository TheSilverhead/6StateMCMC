struct LossLog{time,time_titer,data,M,control,ci,titer} <: DiffEqBase.DECostFunction
  t::time #Time points for data
  t_titer::time_titer #Time points for virus count
  data::data # LFC of data
  measured::M #Indices in model with data
  mock::control #raw data for control
  interval::ci #95% Confidence intervals
  vtiter::titer #Viral titer estimates
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
interval = f.interval
vtiter = f.vtiter

  if sol.retcode != :Success
      return Inf
  end

  subSol = sol(t)

  sumsq = 0.0
  maxlog2(x) = log2(max(0,x)+1)
  @inbounds for i = 1:length(subSol) #Loop through each data time point
    for (j,m) in enumerate(measured) #loop though states
      foo=(maxlog2(subSol[m,i]) - maxlog2(mock[i,j]))
      refmax = data[i,j]+interval[i,j]*0
      refmin = data[i,j]-interval[i,j]*0
      sumsq +=(maxlog2(subSol[m,i]) - maxlog2(mock[i,j]) - data[i,j])^2
      if (foo>refmax)
        sumsq +=(foo-refmax)^2 #Above CI error
      elseif (foo<refmin)
        sumsq +=(foo-refmin)^2 #Below CI error
      end

    end
  end

subSol = sol(t_titer)

@inbounds for i = 1:length(subSol) #Loop through each virus time point
  sumsq += 10*(log10(titer[i,1])-log10(subSol[10,i]))^2 #Virus ssqe, weighted
end

return sumsq
end
