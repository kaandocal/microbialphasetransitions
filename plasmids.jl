using Random
using Distributions
using UnPack 
using StatsBase

include("sim.jl")

mutable struct Cell
    k::Int 
    t::Float64 
    t_b::Float64
    t_d::Float64
    k1::Int
end 

Cell(k::Int, t::Float64=0.) = Cell(k, t, t, Inf, -1)

get_n_offspring(cell::Cell, model) = 2

function sample_offspring(anc::Cell, t, model)
    if anc.k1 == -1 
        anc.k1 = rand(Binomial(2 * anc.k, 0.5))
        Cell(anc.k1, t)
    elseif anc.k1 >= 0
        k2 = 2 * anc.k - anc.k1
        Cell(2 * anc.k - anc.k1, t)
    end 
end 

function get_τ(cell::Cell, model)
    1 + model.β * cell.k + model.α * (cell.k == 0)
end

function simulate_agent!(cell::Cell, tmax, model)
    @assert cell.t <= tmax 

    if !isfinite(cell.t_d)
        cell.t_d = cell.t + get_τ(cell, model)
    end 
    
    cell.t = min(cell.t_d, tmax)
end 

# Utility functions 
""" Simulate a realisation of the population and record its state at regular intervals """
function simulate(N, k0, tmax, model; kwargs...)
    u = [ Cell(k0) for i in 1:N ] 
    pop = Population(u)
        
    simulate!(pop, tmax, model; kwargs...)
    
    pop
end 