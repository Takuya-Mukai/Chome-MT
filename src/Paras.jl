"""
Parameters of system and particles 
    Dr: rotational diffusion constant
    vel: self-propelling velocity
    J: alignment strength
    epsilon: 
    cutoff: cutoff of pariwise interaction
    L: length of simulation box
    sigma: radius of particle
    A = source rate of chemical produced
"""

mutable struct Paras 
    Dr::Float64
    vel::Float64
    J::Float64
    epsilon::Float64
    cutoff::Float64
    L::Float64
    sigma::Float64
    A::Float64
    dx::Float64
end

function initParas(;Dr=1, vel=1, J=1, epsilon=1, cutoff=0.05, L=1.0, sigma=0.001, A = 1, dx=0.2)
    return Paras(Dr, vel, J, epsilon, cutoff, L, sigma, A, dx)
end



