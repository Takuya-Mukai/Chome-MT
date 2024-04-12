using Pkg ; Pkg.activate(".")

using StaticArrays
using Statistics
using LinearAlgebra
using CellListMap.PeriodicSystems
using Parameters 
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
import CellListMap.wrap_relative_to
using JLD2
using GLMakie

include("Paras.jl")
include("Interaction.jl")
include("Simulation.jl")
include("Relaxation.jl")
include("vis.jl")

