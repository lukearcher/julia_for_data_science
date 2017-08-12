include("sample.jl")

function KFCV{T1<:Any, T2<:Any}(I::Array{T1,2}, O::Array{T2,1}, k::Int64 = 10)
	O_ = O
	P = Array(Any, k) # training inputs
	PT = Array(Any, k) # testing inputs
	T = Array(Any, k) # training outputs
	TT = Array(Any, k) # testing outputs
	N = size(I, 1)	
	s = round(Int64, N / k)	
	IND = collect(1:N)

	for i = 1:(k-1)
		ind, trash = sample(IND, O_, s)
		ind = ind[:]
		IND = setdiff(IND, ind)	
		ind_ = setdiff(1:N, ind)
		O_ = O[IND]
		P[i] = I[ind_,:]
		T[i]	 = O[ind_]		
		PT[i] = I[ind,:]
		TT[i] = O[ind]
	end
	
	IND_ = setdiff(1:N, IND)
	P[k] = I[IND_,:]
	T[k] = O[IND_]
	PT[k] = I[IND,:]
	TT[k] = O[IND]
	return P, T, PT, TT
end

function KFCV{T1<:Any, T2<:Any}(I::Array{T1,1}, O::Array{T2,1}, k::Int64)
	J = I[:,:]
	return KFCV(J,O,k)
end
