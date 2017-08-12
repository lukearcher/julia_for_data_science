using Graphs

function dist(X::Array{Float64,2}, x::Array{Float64,2}, i::Int64)
    d = Array(Float64, i)

    for j = 1:i
        d[j] = norm(X[j,:] - x)
    end

    return d
end

function dmatrix(X::Array{Float64,2})
    N = size(X,1)
    D = zeros(N,N)

    for i = 1:(N-1)
        i1 = i + 1
        D[i,i1:end] = dist(X[i1:end,:], X[i,:], (N - i))
        D[i1:end,i] = D[i,i1:end]
    end
    
    return D
end

function MSTofClass(D)
    # distance matrix D doesn't have to be square, since the last row isn't used anywhere
    n = size(D,2)
    N = convert(Int64, n*(n-1) / 2)
    w = Array(Float64, N)    
    c = 0

    for i = 1:(n-1)
        for j = (i+1):n
            c += 1
            w[c] = D[i,j]
        end
    end

    g = simple_complete_graph(n, is_directed=false)
    return kruskal_minimum_spantree(g, w)
end

function MSTC(P, T, PT) # Minimum Spanning Tree Classifier        
    Q = sort(unique(T))
    q = length(Q)
    C = Array(Any, q)
    n = length(T)
    t = eltype(T)
    nn = size(PT, 1)
    z = collect(1:n)
    W = Array(Float64, q)
    Wx = Array(Float64, q)
    y = Array(t, nn) # outputs (predicted labels)
	p = Array(Float64, nn) # probabilities (confidence of classification)
    ind = Array(Any, q)
    z = Array(Int64, q)
    D = dmatrix(P)

    for i = 1:q        
        temp = (T .== (ones(n)*Q[i])).*(1:n)
        ind[i] = temp[temp .> 0]
        a, b = MSTofClass(D[ind[i],ind[i]])
        W[i] = mean(b)
        z[i] = length(ind[i])
    end

    for i = 1:nn
        for j = 1:q
            d = sqrt(sum((P[ind[j], :] - repmat(PT[i,:], z[j], 1)).^2, 2)) # distance of test point i from elements of class j
            DD = hcat(D[ind[j], ind[j]], d)		
            a, b = MSTofClass(DD)
            Wx[j] = mean(b)
        end		
		
        d = W - Wx # improvement in average edge length (weights) of MST
		index = indmax(d)
        y[i] = Q[index]
		d_ = max(d, zeros(q))
		
		if d_[index] == 0
			p[i] = 1 / q
		else
			p[i] = d_[index] / sum(d_)
		end
    end
    
    return y, p
end
