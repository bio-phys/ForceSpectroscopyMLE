#= Read in all datasets in directory "dir".  
Each dataset is sorted in ascending order with respect to the rupture forces (first column). =#

function read_data(dir)
    files = readdir(dir)
    files = filter(!startswith(".") ∘ basename, files)
    N = length(files)
    data = Array{Array{Float64,2},1}(undef, N)
    for i = 1 : N
        tmp_data = readdlm(string(dir,files[i]))
        q = sortperm(tmp_data[:,1])
        data[i] = tmp_data[q,:]
    end
    return data
end

function read_data_unsorted(dir)
    files = readdir(dir)
    files = filter(!startswith(".") ∘ basename, files)
    N = length(files)
    data = Array{Array{Float64,2},1}(undef, N)
    for i = 1 : N
        data[i] = readdlm(string(dir,files[i]))
    end
    return data
end



#= Remove the last i datapoints from each dataset: =#

function reduce_data(data,i)
    N = length(data)
    new_data = Array{Array{Float64,2},1}(undef, N)
    for j = 1 : N
        new_data[j] = data[j][1:end-i,:]
    end
    return new_data
end



#= Remove randomly i datapoints from each dataset: =#

function random_reduce_data(data,i)
    N = length(data)
    new_data = Array{Array{Float64,2},1}(undef, N)
    for j = 1 : N
        M = size(data[j],1)
        r = sample(1:M, i, replace = false)
        new_data[j] = data[j][1:M .∉ [r],:]
    end
    return new_data
end



#= Exclude i = 1, 2, ..., i_max datapoints from each dataset and analyze the reduced data: =#

function data_exclusion(data, ν, i_max, T=295, βΔE_range=(0.1,100.0), Δx_b_range=(0.001,10.0), msteps=100000, mode=:silent, psize=50, tint=60.0)
    N = length(data)
    all_data = vcat(vcat(data...)...)
    parameters = Array{Float64}(undef, 4, 0)
    p = Progress(i_max+1, 1)
    for i = 0 : i_max
        new_data = vcat([ vcat(reduce_data(data[j],i)...) for j = 1 : N ]...)
        est = MLE_estimator(new_data,ν,T,βΔE_range,Δx_b_range,msteps,mode,psize,tint)
        percentage = 100*(size(all_data,1) - size(new_data,1))/size(all_data,1)
        parameters = hcat([parameters, vcat([percentage,est]...)]...)
        next!(p)
    end
    return parameters
end

function random_data_exclusion(data, ν, i_max, T=295, βΔE_range=(0.1,100.0), Δx_b_range=(0.001,10.0), msteps=100000, mode=:silent, psize=50, tint=60.0)
    N = length(data)
    all_data = vcat(vcat(data...)...)
    parameters = Array{Float64}(undef, 4, 0)
    p = Progress(i_max+1, 1)
    for i = 0 : i_max
        new_data = vcat([ vcat(random_reduce_data(data[j],i)...) for j = 1 : N ]...)
        est = MLE_estimator(new_data,ν,T,βΔE_range,Δx_b_range,msteps,mode,psize,tint)
        percentage = 100*(size(all_data,1) - size(new_data,1))/size(all_data,1)
        parameters = hcat([parameters, vcat([percentage,est]...)]...)
        next!(p)
    end
    return parameters
end
