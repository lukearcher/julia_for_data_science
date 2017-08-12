function normalize{T<:Real}(x::Array{T,1}, method = "stat")
	# normalize a variable to the interval [0, 1] or (0, 1)
	
	method = lowercase(method)
	
	if method == "linear"
		m = minimum(x)
		y = (x - m) / (maximum(x) - m)
	elseif method == "stat"
		mx = mean(x)
		sx = std(x)
		y = (x - mx) / sx
	else
		error("invalid normalization method! Please select either ""stat"" or ""linear""!")
	end

	return y
end

function normalize{T<:Real}(X::Array{T,2}, method = "stat")
	N, n = size(X)
	Y = Array(Float64, N, n)
	
	for i = 1:n
		x = X[:,i]
		y = normalize(x, method)
		Y[:,i] = y
	end
	
	return Y
end
