include("sample.jl")

function DID{T<:Real}(I::Array{T,2}, O::Array{Any,1}, s::Int64 = 0)
	I = float(I)	
	N, na = size(I)
	
	if N < na
		N, na = na, N
	end
	
	if s == 0
		s = N
	end
	
	s = min(10000, N)
	
	if na > 1
		for i = 1:na
			X = I[:,i]
			mp = (maximum(X) + minimum(X)) / 2
			if !within(mp, 0.45, 0.55)  # feature is not normalized
				I[:,i] = normalize(X)
			end
		end
	end
	
	q, ~, C, n = nc(O)
	c = zeros(Float64, q, na) # Centers of hyperspheres
	R = zeros(Float64, q) # Radii of hyperspheres
	ICD = zeros(q,q) # Inter Class Discernibility
	sc = zeros(Int64, q) # Samples of each Class
	S = Array(Any, q) # actual Samples

	for i = 1:q
		c[i,:] = mean(I[C[i],:],1) # centers of classes
		cc = rm(c[i,:],n[i])
		D = sum( (cc - I[C[i],:]).^2 , 2) # distances to center of each class
		R[i] = sqrt( maximum(D) ) # Radius of each hypersphere
		sc[i] = round(Int64, s*n[i]/N) # size of sample of class i
		x = copy(C[i]) # indexes of class i
		z = n[i] # number of patters in class i
		temp = zeros(sc[i], na)
		
		r = sort(unique(rand(1:sc[i], z)))
		lr = length(r)	
		temp[1:lr,:] = I[x[r],:]
		deleteat!(x, r)
		z -= lr
		
		for j = (lr+1):sc[i]
			r = ceil(Int64,z*rand())
			temp[j,:] = I[x[r],:]
			deleteat!(x, r)
			z -= 1
		end
		
		S[i] = copy(temp)
	end

	R += eps()
	
	for i = 1:(q-1)		
		for j = (i+1):q
			TEMP = 0.0
			
			for k = 1:sc[i]
				temp = sqrt(sum( (S[i][k,:] - c[j,:]).^2, 2)) / R[j]
				TEMP += min(temp, 1)
			end
			
			for k = 1:sc[j]
				temp = sqrt(sum( (S[j][k,:] - c[i,:]).^2, 2)) / R[i]
				TEMP += min(temp, 1)
			end
			
			ICD[i,j] = TEMP[1] /  (sc[i] + sc[j])
		end
	end
	
	ICD = ICD + ICD'
	y = sum(ICD) / (q^2 - q)

	return y, ICD
end

function DID{T<:Real}(I::Array{T,1}, O::Array{Any,1}, s::Int64 = 0)
	J = I[:,:]
	return DID(J,O,s)
end
