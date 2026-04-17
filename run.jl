# run.jl - run a simulation of the plasmid model in a population of fixed size
# Note: this can be *slow*.

using JLD2 
using ArgParse

include("plasmids.jl")

s = ArgParseSettings()

@add_arg_table! s begin
    "--tmax"
        help = "Simulation time"
        arg_type = Float64
	    required = true
    "--save-every"
        arg_type = Float64 
        required = true
    "--N"
        arg_type = Int 
        required = true
    "a"
        arg_type = Float64
    "b"
        arg_type = Float64
end

args = ArgParse.parse_args(s)

N = args["N"]

tmax = args["tmax"]
α = args["a"]
β = args["b"]

function load_or_create(; k0 = 20, id=1, kwargs...)
    fn = DrWatson.savename(Dict(kwargs))
    OUT = "data/plasmid$id/$fn.jld2"

    data = nothing 
    if !args["restart"]
        try 
            data = load(OUT)
            @info "Loaded previous data from file $OUT..."
        catch
            @info "Could not load from file $OUT, restarting..." 
        end 
    end 

    if data == nothing
        @unpack α, β, N = kwargs
        @info "Starting new simulation"
        data = Dict("pop" => Population([ Cell(k0) for i in 1:N ]),
                    "model" => (; τ_0 = log(2), α, β))
    end

    data, OUT
end 

@info "N = $N, α = $α, β = $β"
    
pop = Population([ Cell(100) for i in 1:N ])
model = (; α, β)

fn = "data/N=$(N)_α=$(α)_β=$(β).jld2"

@info "Running simulations..."
tt = unique([ 0:args["save-every"]:tmax; tmax ])

for t in tt
    @info "Running to time $t"
    simulate!(pop, t, model; save_ancestors=true, trim_every=10.);

    @info "Saving..."
    @save fn pop model
end 

@info "Done!"
