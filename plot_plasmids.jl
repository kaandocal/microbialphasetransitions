using JLD2 
using CairoMakie

include("plasmids.jl")

N = 10000
α = 0.05
ββ = logrange(1e-4, 1e-2, length=65)

ΛΛ = zeros(length(ββ))
kk_m = zeros(length(ββ))

for (i, β) in enumerate(ββ)
    fn = "data/N=$(N)_α=$(α)_β=$(β).jld2"
    pop = load(fn, "pop")

    ΛΛ[i] = (log(length(pop.u) / N) + pop.log_f) / pop.t
    kk = [ cell.k for cell in ancestors(pop, pop.u[1]) ]
    kk_m[i] = mean(kk)
end 

let 
    fig = Figure()

    axΛ = Axis(fig[1,1], xscale=log10, ylabel="growth rate Λ")
    axs = Axis(fig[2,1], xscale=log10, ylabel="selection differential σ")
    axkm = Axis(fig[3,1], xscale=log10, ylabel="mean plasmid number k̄",
                xlabel="metabolic burden β")

    scatter!(axΛ, ββ, ΛΛ)

    ss = -ββ .* ΛΛ .* kk_m

    scatter!(axs, ββ, ss)
    scatter!(axkm, ββ, kk_m)

    fig 
end 
