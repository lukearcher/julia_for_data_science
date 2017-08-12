# Sampling function (stratified sampling)

function nc{T<:Any}(O::Array{T,1})
	Q = unique(O)
	q = length(Q)
	C = Array(Any, q, 1)
	N = Array(Int64, q, 1)
	Z = collect(1:length(O))
	
	for i = 1:q
		# C[i] = findin(O, [Q[i]])
		C[i] = Z[O .== Q[i]]
		N[i] = length(C[i])
	end
	
	return q, Q, C, N
end

function sample_{TI<:Any, TO<:Any}(I1::Array{TI,2}, O1::Array{TO,1}, s::Int64 = 0)
	N, nf = size(I1)
	c, cv, C, cvp = nc(O1) # cvp = class values population
	d = 0 # counter
	sc = zeros(Int64, c) # Samples of each Class
	if s == 0; s = round(Int64, N/10); end
	T1 = typeof(I1[1,1])
	T2 = typeof(O1[1])
	I = Array(T1, s, nf)
	O = Array(T2, s)
	S = Array(Int64, s) # Sample indexes
	
	for i = 1:c
		z = cvp[i]
		sc[i] = round(Int, s*z/N)
		x = copy(C[i])
		
		for j = 1:sc[i]
			d += 1
			r = rand(1:z)
			S[d] = x[r]
			deleteat!(x, r)
			z -= 1
		end
	end
	
	I[1:s, :] = I1[S,:]
	O[1:s] = O1[S]
	return I, O
end

function sample{TI<:Any, TO<:Any}(I1::Array{TI,1}, O1::Array{TO,1}, s::Int64 = 0)
	I2 = I1[:,:]
	return sample(I2, O1, s)
end

function sample{TI<:Any, TO<:Any}(I::Array{TI,2}, O::Array{TO,1}, s::Int64 = 0)
	N, nf = size(I)
	c, whatever, C, cvp = nc(O) # cvp = class values population
	S = Array(Int64, s) # Sample indexes
	z = 0
	
	for i = 1:c
		sc = round(Int64, s*cvp[i]/N) # samples of class i to be taken
		S[(z+1):(z+sc)] = C[i][randperm(cvp[i])][1:sc]
		z += sc
	end
	
	return I[S,:], O[S]
end
