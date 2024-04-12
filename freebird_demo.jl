### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ad6b7485-5af8-4cd3-867a-dc3e53cf9cc0
begin
	# begin
	#     import Pkg
	#     # careful: this is _not_ a reproducible environment
	#     # activate the global environment
	#     Pkg.activate("/Users/myang/Work/git/Demo-FreeBird/Demo")
	# using Plots, PlutoUI, LinearAlgebra
	# using Unitful
	
	# import PlutoUI: combine
	
	# plotly()
	
	# end
	import Pkg; Pkg.add(url="https://github.com/wexlergroup/FreeBird.jl",rev="feature/analysis-and-debug")
	using Plots, PlutoUI, LinearAlgebra
	using Unitful
	using FreeBird

	plotly()
end

# ╔═╡ 82211f2c-c971-48f8-a7a3-e980da7e53c0
html"""<style>
main {
    max-width: 900px;
	font-size: large;
}
"""

# ╔═╡ 4e860050-5b3a-494a-8ce8-e7b5a3ebedf2
md"""
# FreeBird.jl Demo
"""

# ╔═╡ e32b01e1-9c4d-4dd9-85e9-286fd16ed0fa
md"""
# Outline
- Lennard-Jones Potential
- What is a walker?
- Monte Carlo and Nested sampling
"""

# ╔═╡ 42a7c952-e8e5-43ff-8247-520d6fad4be8
md"""
## Lennard-Jones Potential
"""

# ╔═╡ 8941699d-714b-4072-aadd-f15fb2e4adb4
md"""
$V_\mathrm{LJ}(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$
"""

# ╔═╡ c61a733f-8a45-4481-a866-7b07d069cfd3
md"""
- ``\epsilon`` = $(@bind epsilon Slider(1:0.1:15)), defines how deep your potential well is
- ``\sigma`` = $(@bind sigma Slider(1:0.1:4)), at which the potential energy is zero ('size of the particle')
- `cutoff` = $(@bind cutoff Slider(1:0.1:4,default=4)) ``\sigma``, the potential energy is set to zero after the cutoff
- `shift`? $(@bind shift CheckBox()) If true, the whole potential will be shifted at the cutoff
"""

# ╔═╡ 3187ac88-e235-48fb-ad52-12b3fd7709ea
lj_potential = LJParameters(epsilon=epsilon,sigma=sigma,cutoff=cutoff,shift=shift)

# ╔═╡ ff1cf1c9-c2c9-4218-8b9d-a8c826c7d898
plot((0.8:0.01:4)*u"Å", r -> lj_energy(r,lj_potential), ylimits=(-epsilon-0.5,2))

# ╔═╡ e8ad89be-bbc3-4c24-800c-b3dd1400c82d
md"""
```julia
struct LJParameters
    epsilon::typeof(1.0u"eV")
    sigma::typeof(1.0u"Å")
    cutoff::Float64
    shift::typeof(0.0u"eV")
end
```
"""

# ╔═╡ 25148d81-f849-4379-bb47-a0f26af0e257
md"## Let's build some walkers"

# ╔═╡ fe5a1609-c37e-4d75-bbe5-813454091d83
md"""
But before that, let's define a new set of Lennard-Jones parameters.
"""

# ╔═╡ 2216fc0e-d8c3-4cfa-a003-f54701fcea66
lj = LJParameters()

# ╔═╡ 161a24f5-e5f1-4174-8236-ccb57d75063c
md"""
Why does this line work?
"""

# ╔═╡ 257840d3-1004-418b-b732-200bc65fb3a0
md"""
```julia
function LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
    if isa(shift,Bool) && shift
        shiftenergy = lj_energy(Float64(epsilon)*u"eV", Float64(sigma)*u"Å", cutoff*Float64(sigma)*u"Å")  
    else
        shiftenergy = Float64(shift)*u"eV"
    end
    return LJParameters(epsilon*u"eV", sigma*u"Å", cutoff, shiftenergy)
end
```
"""

# ╔═╡ 0d37c8c6-f858-4e58-bf41-aae0a6bee7ad
md"""
## Now let's define some walkers
"""

# ╔═╡ c74d16b3-8424-4ba3-bb0c-d7434cabddcc
md"""
First, let's generate some random initial configurations, with 6 hydrogen atoms in a box:
"""

# ╔═╡ 9e91b9d0-6f3c-4a12-8bab-9c8b6f28857e
hydrogens = generate_initial_configs(5::Int64, 100.0::Float64, 6::Int64; particle_type=:H) 

# ╔═╡ de82ef79-26da-4863-bbee-82e65d01be80
md"""
They are configurations, not walkers yet...
"""

# ╔═╡ 391aa43e-a4b3-4e23-a4eb-55c34af2cd51
md"""
## Now let's define some walkers
"""

# ╔═╡ 96978dc5-fa1e-4947-8cab-b969a67468a6
md"""
Warpping the configurations into the `AtomWalker` objects:
"""

# ╔═╡ 205e7cf1-c435-4433-8f9a-9f833197b5e4
walkers = [AtomWalker(config) for config in hydrogens] # comprehension

# ╔═╡ 552a77eb-f578-4497-a580-674ad9702d44
md"""
## Define a liveset
"""

# ╔═╡ ca3cd8bd-3057-4624-b834-35d4e76545b6
md"""
We have previosuly defined a set of walkers, and a set of Lennard-Jones parameters. Let's combine them to make a `liveset`.
"""

# ╔═╡ 3ff30da5-b855-4b74-9d6f-81c7f5f610c4
liveset = LJAtomWalkers(walkers,lj)

# ╔═╡ 9418bea4-cd5b-47eb-97fd-d3e03ced0a8c
md"""
Now we have everything needed for a nested sampling calculation!
"""

# ╔═╡ 30095fa5-cc4b-4260-892e-f5defa326a3c
md"""
## Nested Sampling!
"""

# ╔═╡ eed6a262-f804-419d-ba96-9db596202508
mc = MCRandomWalkClone()

# ╔═╡ 8a2920e2-9058-4e29-8560-8a249530178f
ns_params=NestedSamplingParameters(200::Int64, 0.1::Float64, 0.01::Float64, 1e-5::Float64, 1.0::Float64, 0::Int64, 200::Int64)

# ╔═╡ f8d5538f-0621-4c11-875d-729af53123da
n_iters = 20_000

# ╔═╡ 6999a230-82d1-480d-a0fc-a7bdada6fdb6
save = SaveEveryN(n=100) # no saving until the last step

# ╔═╡ 354c14c8-df26-460e-bb9b-c3945561d81d
md"""
## Run the code
"""

# ╔═╡ 68f4f490-f80a-4269-9000-e43ae269544a
begin
	if isfile("output_df.csv")
		run(`rm output_df.csv output.ls.extxyz output.traj.extxyz`) # delete previous output
	end
	new_walkers = AtomWalker.(generate_initial_configs(120::Int64, 100.0::Float64, 6::Int64))
	new_ls = LJAtomWalkers(new_walkers, lj)
	energies, live, _ = nested_sampling_loop!(new_ls, ns_params, n_iters, mc, save)
end

# ╔═╡ 1beb8de4-bf20-44d3-8ff8-19a5dcfa6f94
energies

# ╔═╡ f2394846-280f-4c1e-9a82-5c78711452d7
new_ls

# ╔═╡ 48b0f871-293d-44e2-9ffa-12481f187237
md"""## Analysis"""

# ╔═╡ 34f93745-4823-4b50-9b1b-7de0df5d7ef3
gi = gamma_factors(energies.iter::Vector{Int64}, 120::Int64) 

# ╔═╡ 0eed5c28-2152-4294-877a-008912486b98
partition_function(0.1::Float64, gi::Vector{Float64}, energies.emax::Vector{Float64}) 

# ╔═╡ 7e267448-9bfa-4897-a117-ed75ea5e433e
kb = 8.617333262e-5 # eV/K

# ╔═╡ 601be34c-11fd-4693-9098-0eb1e4d8b0c1
ts = 10:10:10000 

# ╔═╡ 4f2a5238-b644-445d-ac33-c521d402ed41
beta = 1 ./(kb.*ts)

# ╔═╡ c58d8710-8f71-4b52-a0cf-b14760e19afb
dof = 18

# ╔═╡ 1aac25a8-37cc-469d-ada4-39a027c7373c
ei = energies.emax .- minimum(energies.emax)

# ╔═╡ 729b7408-e46e-4134-8cf6-9e7735d94840
zs = [partition_function(b, gi, ei) for b in beta]

# ╔═╡ 22169a82-0b5a-4012-8a85-64051e4ed6e9
plot(ts, zs, xlabel="Temperature (K)", ylabel="Z")

# ╔═╡ 72329d32-e194-4090-9f41-0950adad6748
u = [internal_energy(b, gi, ei) for b in beta]

# ╔═╡ 7756264e-7e71-4d12-8a09-baeab1ee4077
plot(ts, u, xlabel="Temperature (K)", ylabel="Internal energy")

# ╔═╡ 1b7d86b4-9d84-44a6-852f-cceea593b2e9
cvs = [cv(b, gi, ei, dof) for b in beta]

# ╔═╡ 39386f04-3539-4b1f-8dec-612b081ffa46
plot(ts, cvs, xlabel="Temperature (K)",ylabel="Heat Capacity")

# ╔═╡ 6071da62-0458-4fa3-8059-f3d337aaced2
md"""
![LJ6]("LJ6.gif")
"""

# ╔═╡ 0bfc63b6-6dc3-4ca4-85a9-7f2028dcaeb1
md"""
## How about surfaces?
"""

# ╔═╡ 6a11ff64-0ad8-4b57-a8db-06765238bfac
surf_energies = read_output("08_surf.csv")

# ╔═╡ 77f8c9fd-9634-4921-a3f9-c0aa5bc6520f
surf_gi = gamma_factors(surf_energies.iter, 640)

# ╔═╡ de6e92b6-d971-4567-b1bf-615f267bd5c8
surf_ei = surf_energies.emax .- minimum(surf_energies.emax)

# ╔═╡ 26a6ea9f-6d59-452f-a8b8-6cb5fd91b7e1
surf_ts = 1:2500

# ╔═╡ f2f8953e-15ba-4261-822f-61b048637880
surf_beta = 1 ./(kb.*surf_ts)

# ╔═╡ df5e0581-bbe3-4522-86e0-7e736a4c3d00
surf_cvs = [cv(b, surf_gi, surf_ei, 3*8) for b in surf_beta]

# ╔═╡ 27ff4b7e-a10f-4d2c-accf-cbcd4aa2b24e
plot(surf_ts, surf_cvs, xlabel="Temperature (K)",ylabel="Heat Capacity")

# ╔═╡ 5cb4cb4e-1d02-4057-bccd-7e96e0f68d85
md"""
##
![surf](https://github.com/yangmr04/FreeBird-demo/blob/main/08.gif?raw=true)
"""

# ╔═╡ Cell order:
# ╟─ad6b7485-5af8-4cd3-867a-dc3e53cf9cc0
# ╟─82211f2c-c971-48f8-a7a3-e980da7e53c0
# ╟─4e860050-5b3a-494a-8ce8-e7b5a3ebedf2
# ╟─e32b01e1-9c4d-4dd9-85e9-286fd16ed0fa
# ╟─42a7c952-e8e5-43ff-8247-520d6fad4be8
# ╟─8941699d-714b-4072-aadd-f15fb2e4adb4
# ╟─c61a733f-8a45-4481-a866-7b07d069cfd3
# ╟─3187ac88-e235-48fb-ad52-12b3fd7709ea
# ╟─ff1cf1c9-c2c9-4218-8b9d-a8c826c7d898
# ╟─e8ad89be-bbc3-4c24-800c-b3dd1400c82d
# ╟─25148d81-f849-4379-bb47-a0f26af0e257
# ╟─fe5a1609-c37e-4d75-bbe5-813454091d83
# ╠═2216fc0e-d8c3-4cfa-a003-f54701fcea66
# ╟─161a24f5-e5f1-4174-8236-ccb57d75063c
# ╟─257840d3-1004-418b-b732-200bc65fb3a0
# ╟─0d37c8c6-f858-4e58-bf41-aae0a6bee7ad
# ╟─c74d16b3-8424-4ba3-bb0c-d7434cabddcc
# ╠═9e91b9d0-6f3c-4a12-8bab-9c8b6f28857e
# ╟─de82ef79-26da-4863-bbee-82e65d01be80
# ╟─391aa43e-a4b3-4e23-a4eb-55c34af2cd51
# ╟─96978dc5-fa1e-4947-8cab-b969a67468a6
# ╠═205e7cf1-c435-4433-8f9a-9f833197b5e4
# ╟─552a77eb-f578-4497-a580-674ad9702d44
# ╟─ca3cd8bd-3057-4624-b834-35d4e76545b6
# ╠═3ff30da5-b855-4b74-9d6f-81c7f5f610c4
# ╟─9418bea4-cd5b-47eb-97fd-d3e03ced0a8c
# ╟─30095fa5-cc4b-4260-892e-f5defa326a3c
# ╠═eed6a262-f804-419d-ba96-9db596202508
# ╠═8a2920e2-9058-4e29-8560-8a249530178f
# ╠═f8d5538f-0621-4c11-875d-729af53123da
# ╠═6999a230-82d1-480d-a0fc-a7bdada6fdb6
# ╟─354c14c8-df26-460e-bb9b-c3945561d81d
# ╠═68f4f490-f80a-4269-9000-e43ae269544a
# ╠═1beb8de4-bf20-44d3-8ff8-19a5dcfa6f94
# ╠═f2394846-280f-4c1e-9a82-5c78711452d7
# ╟─48b0f871-293d-44e2-9ffa-12481f187237
# ╠═34f93745-4823-4b50-9b1b-7de0df5d7ef3
# ╠═0eed5c28-2152-4294-877a-008912486b98
# ╠═7e267448-9bfa-4897-a117-ed75ea5e433e
# ╠═601be34c-11fd-4693-9098-0eb1e4d8b0c1
# ╠═4f2a5238-b644-445d-ac33-c521d402ed41
# ╠═c58d8710-8f71-4b52-a0cf-b14760e19afb
# ╠═1aac25a8-37cc-469d-ada4-39a027c7373c
# ╠═729b7408-e46e-4134-8cf6-9e7735d94840
# ╠═22169a82-0b5a-4012-8a85-64051e4ed6e9
# ╠═72329d32-e194-4090-9f41-0950adad6748
# ╠═7756264e-7e71-4d12-8a09-baeab1ee4077
# ╠═1b7d86b4-9d84-44a6-852f-cceea593b2e9
# ╠═39386f04-3539-4b1f-8dec-612b081ffa46
# ╠═6071da62-0458-4fa3-8059-f3d337aaced2
# ╟─0bfc63b6-6dc3-4ca4-85a9-7f2028dcaeb1
# ╠═6a11ff64-0ad8-4b57-a8db-06765238bfac
# ╠═77f8c9fd-9634-4921-a3f9-c0aa5bc6520f
# ╠═de6e92b6-d971-4567-b1bf-615f267bd5c8
# ╠═26a6ea9f-6d59-452f-a8b8-6cb5fd91b7e1
# ╠═f2f8953e-15ba-4261-822f-61b048637880
# ╠═df5e0581-bbe3-4522-86e0-7e736a4c3d00
# ╠═27ff4b7e-a10f-4d2c-accf-cbcd4aa2b24e
# ╟─5cb4cb4e-1d02-4057-bccd-7e96e0f68d85
