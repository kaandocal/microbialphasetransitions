using JLD2 
using CairoMakie

include("plasmids.jl")

αα = 0.01:0.0025:1
ββ = logrange(1e-4, 1e-2, length=65)
N = 10000

ΛΛ = zeros(length(αα), length(ββ))
kk_m = zeros(length(αα), length(ββ))

for (j, α) in enumerate(αα)
    for (i, β) in enumerate(ββ)
        fn = "data/N=$(N)_α=$(α)_β=$(β).jld2"
        pop = load(fn, "pop")

        ΛΛ[j,i] = (log(length(pop.u) / N) + pop.log_f) / pop.t
        kk = [ cell.k for cell in ancestors(pop, pop.u[1]) ]
        kk_m[j,i] = mean(kk)
    end
end 

let 
    fig = Figure()

    axΛ = Axis(fig[1,1], 
               xlabel="metabolic burden β", ylabel="addiction penalty α", 
               title="growth rate Λ")

    axk = Axis(fig[1,3], xlabel="metabolic burden β", ylabel="addiction penalty α",
               title="mean plasmid number k̄")

    hmΛ = heatmap!(axΛ, log10.(ββ), αα, ΛΛ, interpolate=false)
    hmk = heatmap!(axk, log10.(ββ), αα, kk_m, interpolate=false)

    Colorbar(fig[1,2], hmΛ)
    Colorbar(fig[1,4], hmk)

    fig 
end 