#= Critical force: =#

F_crit(β,βΔE,Δx_b,ν) = βΔE/(β*ν*Δx_b)



#= Mean rupture force [Eq. (5) in R. W. Friddle, PRL 100, 138302 (2008)]: =#

function meanF(β,βΔE,Δx_b,k_0,ν,Fdot)
    F_c = F_crit(β,βΔE,Δx_b,ν)
    return F_c*(1.0 - (1.0 + exp(k_0/(β*Δx_b*Fdot))*expinti(-k_0/(β*Δx_b*Fdot))/βΔE)^ν)
end



#= Variance of rupture force [Eq. (6) in R. W. Friddle, PRL 100, 138302 (2008) with improved prefactor]: =#

varF(β,βΔE,Δx_b,k_0,ν,Fdot) = π^2/6*(β*Δx_b*(1.0 + k_0/(β*Δx_b*Fdot)))^(-2)*(1.0 + exp(k_0/(β*Δx_b*Fdot))*expinti(-k_0/(β*Δx_b*Fdot))/βΔE)^(2*ν-2)



# Diffusion coefficient, calculated from Kramers rate associated with ν:
function D(βΔE,Δx_b,k_0,ν)
    if ν == 1/2
        return k_0*exp(βΔE)*Δx_b^2/(2*βΔE*sqrt(βΔE/π))
    elseif ν == 2/3
        return k_0*exp(βΔE)*π*Δx_b^2/(3*βΔE)
    else
        @warn "Invalid value of ν"
    end
end
