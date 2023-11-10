using LinearAlgebra
using Plots

function hilbert(n :: Integer)
    H = [1 / (i + j - 1) for i in 1:n, j in 1:n]
    b = ones(n)
    return H, b
end
