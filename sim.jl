using Random
using Distributions
using DataStructures

mutable struct Population{T}
    u::Vector{T}
    t::Float64 
    log_f::Float64
    parents::Dict{T,T}
end 

Population(u::AbstractVector{T}) where {T} = Population(u, 0., 0., Dict{T,T}())

function ancestors(tree::Population, x)
    ret = typeof(x)[]

    while haskey(tree.parents, x)
        x = tree.parents[x]
        push!(ret, x)
    end 

    reverse!(ret)
end 

function save_lineage!(dict::AbstractDict{T}, tree::Population{T}, x::T) where {T}
    while haskey(tree.parents, x) && !haskey(dict, x)
        anc = tree.parents[x]
        dict[x] = anc 
        x = anc
    end
end  

function trim!(pop::Population{T}) where {T}
    parents = Dict{T,T}()

    for x in pop.u
        save_lineage!(parents, pop, x)
    end 

    pop.parents = parents 
    pop
end 

function maybe_trim!(pop::Population{T}, N) where {T}
    if length(pop.parents) > pop.t * N * 100
        trim!(pop)
    end
end

###

# This seems to make simulations much faster, but can cause issues
# Consider commenting this out if this is the case
Base.hash(x::Int64) = UInt(x)

function simulate!(pop, tmax, model; Nmax=length(pop.u), trim_every=tmax, kwargs...)
    pop.t < tmax || return

    q = PriorityQueue{Int,Float64}()

    for (i, x) in enumerate(pop.u)
        if x.t <= tmax 
            simulate_agent!(x, tmax, model)
            push!(q, i => x.t)
        end
    end

    ttrim = pop.t:trim_every:tmax 

    for t in ttrim[2:end]
        simulate_queue!(pop, q, tmax, model; Nmax, kwargs...)
        trim!(pop)
    end 
end

function simulate_queue!(pop, q, tmax, model; Nmax, save_ancestors=false)
    while !isempty(q)
        i, pop.t = popfirst!(q)
        a = pop.u[i]
        @assert pop.t == a.t

        pop.t >= tmax && break

        if a.t_d > pop.t    # individual is still alive
            simulate_agent!(a, tmax, model)
            push!(q, i => pop.u[i].t)
        else                # individual has reached the end of its lifetime
            n_offspr = get_n_offspring(a, model)

            # remove individual from array
            # simply move the last individual to slot i and shrink the population
            if i != length(pop.u)
                pop.u[i] = pop.u[end]
                delete!(q, length(pop.u))
                push!(q, i => pop.u[i].t)
            end 

            resize!(pop.u, length(pop.u) - 1)

            N_extrap = length(pop.u) + n_offspr
            add_offspring!(pop, q, n_offspr, a, model; Nmax, save_ancestors)
            #maybe_trim!(pop, N)
            pop.log_f += log(N_extrap / length(pop.u))
        end
    end
end 

function add_offspring!(pop::Population, q, n, anc, model; Nmax, save_ancestors=false)
    # If we have less than Nmax individuals in total,
    # add offspring without culling
    for i in length(pop.u)+1:Nmax 
        n == 0 && return 
        offspring = sample_offspring(anc, pop.t, model)
        push!(pop.u, offspring)
        if save_ancestors 
            pop.parents[offspring] = anc
        end 

        @assert i == length(pop.u)
        push!(q, i => pop.t)
        n -= 1
    end 

    # We have more than Nmax individuals in total; randomly discard some of them
    # If we discard an existing individual, we replace it by a newborn
    while n > 0
        @assert length(pop.u) == Nmax
        i = sample(1:Nmax+n)
        n -= 1
        i > Nmax && continue 
        offspring = sample_offspring(anc, pop.t, model)
        pop.u[i] = offspring
        if save_ancestors 
            pop.parents[offspring] = anc
        end 
        
        delete!(q, i)
        push!(q, i => pop.t)
    end
end 
