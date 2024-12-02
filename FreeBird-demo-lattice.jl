### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ dc468f4a-b0ca-11ef-021b-91e7c71a8637
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
	import Pkg; Pkg.add(url="https://github.com/wexlergroup/FreeBird.jl",rev="feature/multi-component-lattice")
	using Plots, PlutoUI, LinearAlgebra
	using Unitful
	using FreeBird

	plotly()
end

# ╔═╡ dde5eba8-70d0-43f5-86e8-ffa15b29d939
html"""<style>
main {
    max-width: 900px;
	font-size: large;
}
"""

# ╔═╡ fe9f5532-4b07-46a6-a19a-3c67bcc739c4
md"""
# FreeBird Demo -- Lattice
"""

# ╔═╡ 923c934b-e9d9-44ca-9415-0fb0b6cd142c
md"""
## Outline
Brief demonstrations on how to set up various lattice systems implemented in `FreeBird`:
- `SLttice` -- single-component lattice object
- `MLattice` -- multi-component lattice object
and how to do samplings with them.
"""

# ╔═╡ 8ae1cf69-76e4-4d04-9214-046bec3ea593
md"## `AbstractLattice`"

# ╔═╡ 60a05192-9663-4150-af02-f585ec059998
md"""
In `Julia`, there is no concept of *objects*, but data with different *types*.
"""

# ╔═╡ 75ab78ff-576c-49bb-9d05-2d9e4bf4bccb
AbstractLattice

# ╔═╡ f56a52fe-b8dd-4357-ad7d-fe717a5609c3
md"""
## `SLattice` -- single-component Lattices
"""

# ╔═╡ 0831efa6-2319-4799-b741-6b9bd92d56cc
md"""
```julia
mutable struct SLattice{G} <: AbstractLattice
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    occupations::Vector{Bool}
    cutoff_radii::Vector{Float64}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}
end
```
"""

# ╔═╡ dc74a5c8-1abc-42fd-ae6f-8089a59744e7
md"""
Here, `SLattice` is a parameterized type. It takes a parameter `G` which define the lattice geometry.
"""

# ╔═╡ cff5edfe-75aa-4139-9df4-2ffbf277ac50
LatticeGeometry

# ╔═╡ e60eae60-8dcd-41f0-b1b4-d2ab80228f4e
md"""
### Let's build a `SLattice`
"""

# ╔═╡ 4d423c0d-6840-48b4-a461-4fb4d8e13b34
sl = SLattice{SquareLattice}(supercell_dimensions=(3,3,1), occupations=[1,2,3,4])

# ╔═╡ fa9d9c94-433e-4fcd-b871-2911f8a28b7d
sl |> typeof

# ╔═╡ db930975-3d35-4e66-b794-7e5b2e8a6136
md"""
## Why does it work? Why don't I need to give all info needed?
"""

# ╔═╡ 5a0bee0b-202a-44b3-9bca-209fa5940ae0
md"Because we define a (outer) `constructor` function:"

# ╔═╡ fd153a30-0a1a-45a1-898c-8a33f6102065
md"""
```julia
function SLattice{SquareLattice}(;lattice_constant::Float64=1.0,
                                      basis::Vector{Tuple{Float64, Float64, Float64}}=[(0.0, 0.0, 0.0)],
                                      supercell_dimensions::Tuple{Int64, Int64, Int64}=(4, 4, 1),
                                      periodicity::Tuple{Bool, Bool, Bool}=(true, true, true),
                                      cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                      occupations::Union{Vector{Int}, Symbol}=[1, 2, 3, 4],
                                      adsorptions::Union{Vector{Int}, Symbol}=:full,
                                      )

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 lattice_constant 0.0; 0.0 0.0 1.0]
    dim = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3]
    lattice_occupations = zeros(Bool, dim*length(basis))
    lattice_adsorptions = zeros(Bool, dim*length(basis))

    if occupations == :full
        lattice_occupations = [true for i in 1:dim*length(basis)]
    else
        for i in occupations
            lattice_occupations[i] = true
        end
    end

    if adsorptions == :full
        lattice_adsorptions = [true for i in 1:dim*length(basis)]
    else
        for i in adsorptions
            lattice_adsorptions[i] = true
        end
    end

    return SLattice{SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_occupations, lattice_adsorptions, cutoff_radii)
end
```
"""

# ╔═╡ cd73ddcb-d001-4754-84e0-89e7b779caf8
md"""
## `MLattice` -- multi-component lattices
"""

# ╔═╡ 77452b81-9fdf-4091-905c-54b68a7751b2
md"""
```julia
mutable struct MLattice{C,G} <: AbstractLattice
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    components::Vector{Vector{Bool}}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}
end
```
"""

# ╔═╡ db87de34-3e97-429c-853f-dbd9c41f7a9b
md"""
Here, `MLattice` is also a parameterized type. It takes a parameters, `G` which define the lattice geometry, and `C` which is the number of components.
"""

# ╔═╡ 095ae60a-646d-4401-92c4-35131e1cfa9b
md"""
### Let's build a `MLattice`
"""

# ╔═╡ 07ec57e9-b029-41a5-aefd-69e88efdad71
ml = MLattice{2,SquareLattice}(supercell_dimensions=(3,3,1))

# ╔═╡ 0c549200-c96a-476f-a91e-437897aae4cc
ml |> typeof

# ╔═╡ 506b0cff-f170-4f71-ac1c-01480ce2dfc6
md"## Now, how to run a sampling calculation?"

# ╔═╡ be90097f-d7c2-47b3-bdb2-87c2abee158f
md"Let's define how a lattice energy should be calculated, via `ClassicalHamiltonian`s."

# ╔═╡ 2490969e-f3e3-457d-a425-3690fd6e3fb7
sham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

# ╔═╡ 0fea9062-dc16-484b-abc4-ec73aa1f4e68
sham |> typeof

# ╔═╡ d87285b9-76a4-4e20-b1de-229772f12d12
initial_lattice = SLattice{SquareLattice}(supercell_dimensions=(3,3,1), occupations=[1,2,3,4])

# ╔═╡ 3084250d-7b36-44aa-8b31-b6430a12186a
md"That's all we need to do exact enumeration:"

# ╔═╡ 516dc32b-359a-493d-9c18-c1f46c565622
sdf, sls = exact_enumeration(initial_lattice, sham)

# ╔═╡ 79fb837c-971a-46d8-80a2-ad147114b9c8
md"""
## Try with `MLattice`
"""

# ╔═╡ 3852aada-1778-45ac-95c8-e24289c9c7c9
initial_mlattice = MLattice{2,SquareLattice}(supercell_dimensions=(3,3,1))

# ╔═╡ e7d55611-9c51-40c1-abf2-df53410f0205
mham = MLatticeHamiltonian(2,[GenericLatticeHamiltonian(-0.04, [-0.01*i^2, -0.0025], u"eV") for i in 1:3])

# ╔═╡ b3910021-3a47-4a4b-8450-bf0005ba659e
# mdf, mls = exact_enumeration(initial_mlattice, mham)

# ╔═╡ 60c9f7c4-b371-402a-9e89-cebf1d53c267
# mls

# ╔═╡ 4d27ec54-24e7-44b5-98c5-f5832f9c93ce
# configs = [walker.configuration for walker in mls.walkers]

# ╔═╡ Cell order:
# ╟─dc468f4a-b0ca-11ef-021b-91e7c71a8637
# ╟─dde5eba8-70d0-43f5-86e8-ffa15b29d939
# ╟─fe9f5532-4b07-46a6-a19a-3c67bcc739c4
# ╟─923c934b-e9d9-44ca-9415-0fb0b6cd142c
# ╟─8ae1cf69-76e4-4d04-9214-046bec3ea593
# ╟─60a05192-9663-4150-af02-f585ec059998
# ╠═75ab78ff-576c-49bb-9d05-2d9e4bf4bccb
# ╟─f56a52fe-b8dd-4357-ad7d-fe717a5609c3
# ╟─0831efa6-2319-4799-b741-6b9bd92d56cc
# ╟─dc74a5c8-1abc-42fd-ae6f-8089a59744e7
# ╠═cff5edfe-75aa-4139-9df4-2ffbf277ac50
# ╟─e60eae60-8dcd-41f0-b1b4-d2ab80228f4e
# ╠═4d423c0d-6840-48b4-a461-4fb4d8e13b34
# ╠═fa9d9c94-433e-4fcd-b871-2911f8a28b7d
# ╟─db930975-3d35-4e66-b794-7e5b2e8a6136
# ╟─5a0bee0b-202a-44b3-9bca-209fa5940ae0
# ╟─fd153a30-0a1a-45a1-898c-8a33f6102065
# ╟─cd73ddcb-d001-4754-84e0-89e7b779caf8
# ╟─77452b81-9fdf-4091-905c-54b68a7751b2
# ╟─db87de34-3e97-429c-853f-dbd9c41f7a9b
# ╟─095ae60a-646d-4401-92c4-35131e1cfa9b
# ╠═07ec57e9-b029-41a5-aefd-69e88efdad71
# ╠═0c549200-c96a-476f-a91e-437897aae4cc
# ╟─506b0cff-f170-4f71-ac1c-01480ce2dfc6
# ╟─be90097f-d7c2-47b3-bdb2-87c2abee158f
# ╠═2490969e-f3e3-457d-a425-3690fd6e3fb7
# ╠═0fea9062-dc16-484b-abc4-ec73aa1f4e68
# ╠═d87285b9-76a4-4e20-b1de-229772f12d12
# ╟─3084250d-7b36-44aa-8b31-b6430a12186a
# ╠═516dc32b-359a-493d-9c18-c1f46c565622
# ╟─79fb837c-971a-46d8-80a2-ad147114b9c8
# ╠═3852aada-1778-45ac-95c8-e24289c9c7c9
# ╠═e7d55611-9c51-40c1-abf2-df53410f0205
# ╠═b3910021-3a47-4a4b-8450-bf0005ba659e
# ╠═60c9f7c4-b371-402a-9e89-cebf1d53c267
# ╠═4d27ec54-24e7-44b5-98c5-f5832f9c93ce
