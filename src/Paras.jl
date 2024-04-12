struct Paras 
    Dr::Float64
    vel::Float64
    J::Float64
    epsilon::Float64
    cutoff::Float64
    L::Float64
    sigma::Float64

end

function initParas(;Dr=1, vel=1, J=1, epsilon=1, cutoff=0.05, L=1.0, sigma=0.001)
    return Paras(Dr, vel, J, epsilon, cutoff, L, sigma)
end
