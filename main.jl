###########################
#Main function that iterates through sampling new parameters
###########################
function ptMCMC(prob::DiffEqBase.DEProblem,alg,priors,parBounds,lossFunc,iterations = 1000)
#This is a Metropolis-Hastings MCMC algorithm designed to optimize ODE parameters to fit a given dataset.

    ###########################
    #Basic setup for the run
    ###########################
    parNum = length(priors) #Number of parameters
    chainNum = 1 #number of chains
    βmin, βmax = [0.1, 1.0] #Bounds for chain temperatures
    β = chainNum==1 ? 1 : range(βmin, βmax, length=chainNum) #assign each chain a temperature

    #Initialize objects to store chains and probabilities
    parCurr = [rand.(priors) for i = 1:chainNum] #Starting parameters
    chain = zeros(iterations,parNum,chainNum) #Initialize chain
    logProb = zeros(iterations,chainNum)  #log-posterior chain

    #Storage for when the proposals are accepted
    acceptRate = zeros(iterations,chainNum) #Acceptance ratio
    acceptNum = zeros(Int64,chainNum) #Number of acceptances

    #Generate the log likelihood function (need to include std dev of data)
    SumSqErr = build_loss_objective(prob,alg,lossFunc)
    LogLikeLossFunc(pars) = -0.5*SumSqErr(pars)
    logProbCurr = LogLikeLossFunc.(parCurr)

    #Create proposal distributions for each parameter to general new guesses
    σ = 10 # "Tuning parameter", i.e. how big are your jumps

    #Swap Information
    swapRatio = zeros(iterations) #Acceptance ratio of swaps
    swapNum = 0 #Number of swap acceptances

    #Input the first sample for the chain and log probability
    chain[1,:,:] = hcat(parCurr...) #hcat: array of vectors → 2D array
    logProb[1,:] = logProbCurr

    ###########################
    #Main Loop: Go through number of samples desired
    ###########################
    @showprogress 1 "Computing MCMC Samples..." for i = 2:iterations
        ###########################
        #Apply MH_MCMC step to each chain
            #(similar to lapply in R or map in python)
        ###########################
        #for pmap I'll have to be careful about rand...
        chainResult = MH_MCMC.(parCurr,logProbCurr,acceptNum,β,σ,[parBounds],LogLikeLossFunc)

        #Store the results from the MCMC step
        for (chnIndex, chn) in enumerate(chainResult)
            #Update the current parameters
            parCurr[chnIndex] = chn[1]
            logProbCurr[chnIndex] = chn[2]
            acceptNum[chnIndex] += chn[3]

            #Update the chains storage containers
            chain[i,:,chnIndex] = parCurr[chnIndex]
            logProb[i,chnIndex] = logProbCurr[chnIndex]
            acceptRate[i,chnIndex] = acceptNum[chnIndex] / (i-1)
        end

        ###########################
        #Let's swap some chains (going from highest to lowest temp chain)
        ###########################
    if chainNum > 1
        for chnIndex = chainNum:-1:2
            #Calculate the error and temperature change between chains
            ΔE = logProbCurr[chnIndex] - logProbCurr[chnIndex-1]
            Δβ = β[chnIndex] - β[chnIndex-1]

            #Check if the swap is g0od enough
            if ΔE*Δβ > rand(LogNormal(0,1))
                #swapIndex = chnIndex:-1:chnIndex-1
                swapIndex = reverse(chnIndex-1:chnIndex) #idk if more readable
                #Swap the parameters
                parCurr[reverse(swapIndex)] = parCurr[swapIndex]
                #Swap the error/energy
                logProbCurr[reverse(swapIndex)] = logProbCurr[swapIndex]

                #Count the swap
                swapNum += 1
            end
        end
        #Save swap percent
        swapRatio[i] = swapNum / (i-1)
    end

end #iterations

#return the final results of the chains
return chain, logProb, acceptRate

end #parSample2


###########################
#MH algorithm for one chain
###########################
function MH_MCMC(parCurr,logProbCurr,acceptNum,β,σ,parBounds,LogLikeLossFunc)
    #Propose a new parameter set
    parPropose = NewParameter!(parCurr,parBounds,σ)
    #Calculate change in error
    logProbPropose = LogLikeLossFunc(parPropose)
    ΔE = logProbPropose - logProbCurr
    #Calculate the acceptance ratio correction
    ΔC = cdfCalculator(parPropose,parBounds,σ) - cdfCalculator(parCurr,parBounds,σ)

    #Calculate the acceptance ratio
    acceptRatio = min(1,exp(β*ΔE + ΔC))


    if rand() < acceptRatio #If the acceptance passes...
        accept = 1
        #Update the current parameters with the proposed parameters
        return parPropose, logProbPropose, accept
    else
        accept = 0
        #Rejected, return the same values inputted
        return parCurr, logProbCurr, accept
    end

end #MH_MCMC


###########################
#Given the current parameter set, generate a new parameter set
###########################
function NewParameter!(parameters,parBounds,σ=0.3)
    #Add random normal number to current Parameter set
    #parameters += rand(Normal(0,σ^2), length(parameters))
    parameters = [rand(TruncatedNormal(p,σ^2,lb,ub)) for (p,(lb,ub)) in zip(parameters,parBounds)]
end

###########################
#Given a parameter set, calculate the CDF
#see: https://darrenjw.wordpress.com/tag/truncate/
###########################
function cdfCalculator(pars,bounds,σ=0.3)
    p = logcdf.([TruncatedNormal(0,σ^2,lb,ub) for (lb,ub) in bounds], pars)
    return sum(p)
end
