function lu(A :: Matrix)
	n = size(A,1)
	L = Matrix{Float64}(I(n))
	U = copy(A)

	for j in 1:n-1
			for i in j+1:n
				L[i,j] = U[i,j] / U[j,j]
				U[i,:] -= L[i,j] * U[j,:]
			end
	end
	return L, U
end

function retrosub(A :: Matrix, b ::Vector)
    soma = 0

    n = size(A)[1]
    x = zeros(n)
    x[n] = b[n]/A[n,n]

    rafa = n-1:-1:1
    for i in rafa
        soma = b[i]
        
        for j in i+1:n
            soma += -A[i,j] * x[j]
        end
        
    x[i] = soma/A[i,i]
    end

    return x
end

function forwardsub(A :: Matrix, b ::Vector)
    soma = 0

    n = size(A)[1]
    x = zeros(n)
    x[1] = b[1]/A[1,1]

    rafa = 2:n
    for i in rafa
        soma = b[i]
        
        for j in n:-1:1
            soma += -A[i,j] * x[j]
        end
        
    x[i] = soma/A[i,i]
    end

    return x
end

function solveLU(L :: Matrix, U ::Matrix, b ::Vector)
	y = forwardsub(L,b)
	x = retrosub(U, y)
	return x
end