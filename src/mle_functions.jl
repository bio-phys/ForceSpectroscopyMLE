#= Reduced negative-log-likelihood for the free-energy landscape parameters βΔE and Δx_b: =#

function reduced_likelihood(β,βΔE,Δx_b,ν,data)
    F_c = F_crit(β,βΔE,Δx_b,ν)
    N = size(data,1)
    term_1 = 0.0
    term_2 = 0.0
    @inbounds for n = 1 : N
        F = data[n,1]
        Fdot = data[n,2]
        tmp_1 = 1.0 - F/F_c
        if tmp_1 > 0.0
            tmp_2 = βΔE*(1.0 - tmp_1^(1/ν))
            term_1 += (exp(tmp_2) - 1.0)/(β*Δx_b*Fdot)
            term_2 += (1/ν - 1)*log(tmp_1) + tmp_2
        else
            term_1 += 1.0 # this is to penalize parameter choices with extremely low F_c
            term_2 -= 100000.0
        end
    end
    return log(term_1) - term_2/N
end



#= Maximum likelihood estimator for the instantaneous rate k_0: =#

function k_0(β,βΔE,Δx_b,ν,data)
    F_c = F_crit(β,βΔE,Δx_b,ν)
    N = size(data,1)
    sum = 0.0
    @inbounds for n = 1 : N
        F = data[n,1]
        Fdot = data[n,2]
        sum += (exp(βΔE*(1.0 - (1.0 - F/F_c)^(1/ν))) - 1.0)/(β*Δx_b*Fdot)
    end
    return N/sum
end



#= Compute MLE estimates for the parameters βΔE, Δx_b and k_0: =#

function MLE_estimator(data, ν, T=295, βΔE_range=(0.1,100.0), Δx_b_range=(0.001,10.0), msteps=100000, mode=:compact, psize=50, tint=60.0)
    β = 1/(1.38064852e-23*1e12*1e9*T) # (pN nm)^(-1)
    result = bboptimize(q->reduced_likelihood(β,q[1],q[2],ν,data); 
        SearchRange = [βΔE_range,Δx_b_range], 
        MaxSteps = msteps, 
        TraceMode = mode, 
        PopulationSize = psize, 
        TraceInterval = tint) # optimization in βΔE-Δx_b-space performed using BlackBoxOptim package
    if round(ν,digits=1) == 1.0
        return [ best_candidate(result)[2], k_0(β,best_candidate(result)[1],best_candidate(result)[2],ν,data) ]
    else
        return vcat([ best_candidate(result), k_0(β,best_candidate(result)[1],best_candidate(result)[2],ν,data) ]...)
    end
end



#= Evaluate the uncertainty of an estimate via data resampling: =#

function bootstrapping(data, ν, N=100, T=295, βΔE_range=(0.1,100.0), Δx_b_range=(0.001,10.0), msteps=100000, mode=:silent, psize=50, tint=60.0)
    β = 1/(1.38064852e-23*1e12*1e9*T) # (pN nm)^(-1)
    M = size(data,1)
    βΔE = zeros(N)
    Δx_b = zeros(N)
    k_0 = zeros(N)
    p = Progress(N,1)
    @threads for n = 1 : N
        rnd = rand(1:M,M)
        new_data = data[rnd,:]
        est = MLE_estimator(new_data,ν,T,βΔE_range,Δx_b_range,msteps,mode,psize,tint)
        βΔE[n] = est[1]; Δx_b[n] = est[2]; k_0[n] = est[3];
        next!(p)
    end
    return βΔE, Δx_b, k_0
end

function MLE_errors(data, ν, N=100, T=295, βΔE_range=(0.1,100.0), Δx_b_range=(0.001,10.0), msteps=100000, mode=:silent, psize=50, tint=60.0)
    βΔE, Δx_b, k_0 = bootstrapping(data,ν,N,T,βΔE_range,Δx_b_range,msteps,mode,psize,tint)
    δβΔE = sqrt(var(βΔE))
    δΔx_b = sqrt(var(Δx_b))
    δk_0 = sqrt(var(k_0))
    return [δβΔE, δΔx_b, δk_0]
end
