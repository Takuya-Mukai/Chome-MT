using Pkg ; Pkg.activate(".")

using StaticArrays
using Statistics
using LinearAlgebra
using Interpolations
using CellListMap.PeriodicSystems
using Parameters 
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
import CellListMap.wrap_relative_to
using JLD2
using ProgressMeter
using GLMakie

include("Paras.jl")
include("Interaction.jl")
include("Interaction/Field.jl")
include("Interaction/Viscek.jl")
include("Interaction/viscek_spatio.jl")
include("Interaction/Chemotaxis.jl")
include("Simulation.jl")
include("Relaxation.jl")
include("vis.jl")
include("Analysis.jl")




function vec2comp(vec)
    x = [v[1] for v in vec]
    y = [v[2] for v in vec]

    return x, y
end

angle2dir(ϕ) = SA[cos(ϕ), sin(ϕ)]



function pos2idx(pos, dx)
    idx = round.(Int64, pos / dx) .+ 1
    return idx
end


function idx2pos(idx, dx)
    pos = idx*dx .- dx
    return pos
end

function wrap_index(idx, lim)
    if idx <= 0
        idx = lim + idx
        idx = wrap_index(idx, lim)
    elseif idx > lim
        idx = idx - lim
        idx = wrap_index(idx, lim)
    else
        return idx
    end
end

function wrap_position(x, L)
    x = x - floor(x / L) * L
    return x
end
