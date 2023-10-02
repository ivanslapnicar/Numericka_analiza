### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 34d20879-1b40-498a-83fe-022ec52cdadc
using PlutoUI, Random, LinearAlgebra

# ╔═╡ 732d0a51-841e-4ca7-aca6-32fe38ef90f2
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 54201b21-914a-4152-b187-3ecfe9982e65
md"""
# Norme


Općenito, __norma__ na vektorskom prostoru $X$ je svaka funkcija $\| \phantom{x} \| : X\to \mathbb{R}$ sa sljedećim svojstvima:

1.  $\| x\|=0 \Leftrightarrow x=0$
2.  $\| \lambda x\|=|\lambda| \|x\|$
3.  $\| x+y\| \leq \|x\|+\|y\|  \qquad$ (nejednakost trokuta)

"""

# ╔═╡ 92d26690-7f76-4a73-ac98-38c42f9c46f4
md"""
## Vektorske norme

Za $X=\mathbb{R}^n$ imamo

$$\|x\|_p=\big(\sum_{i=1}^n |x_i|^p\big)^{1/p}$$

Posebno:

*  $\|x\|_1=\sum_{i=1}^n |x_i|$
*  $\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}$
*  $\|x\|_\infty = \max\limits_{i=1,\ldots,n} |x_i|$
"""

# ╔═╡ 1ad65234-9ba8-4eaf-af4d-b881f8d5724d
begin
	Random.seed!(1244)
	x=rand(-9:9,5)
end

# ╔═╡ c1ffcca3-d537-4c4f-86b1-87215325c115
norm(x,1), norm(x), norm(x,Inf)

# ╔═╡ 31a749ca-adb3-4e0e-9484-4506e1586a74
md"""
## Matrične norme

Iz svake vektorske norme možemo izvesti matričnu normu (__inducirane norme__):

$$\|A\| = \max\limits_{x\neq 0} \frac{\|Ax\|}{\|x\|}=\max\limits_{\|x\|=1} \|Ax\|$$

Posebno:

*  $\|A\|_1=\max\limits_{j=1:n} \sum_{i=1}^n |a_{ij}|$  - najveća 1-norma stupca
*  $\|A\|_{\infty}=\max\limits_{i=1:n} \sum_{j=1}^n |a_{ij}|$ - najveća 1-norma retka
*  $\|A\|_2$ - najveća singularna vrijednost  matrice $A$

__Frobeniusova__ ili __Euklidska__ norma

$$\|A\|_F =\sqrt{\sum_{i,j=1}^n a_{ij}^2}$$

nije inducirana norma.

Matrične norme još imaju i svojstvo

$$
\|A\cdot B\|\leq \|A\| \cdot \| B\|.$$
"""

# ╔═╡ 365a36f4-e332-40ea-9113-d322dc47d480
A=rand(-4:4,5,5)

# ╔═╡ 10bc88f1-3675-4446-b31a-033848b6160e
norm(A,1), norm(A), norm(A,2), norm(A,Inf), opnorm(A),
maximum(svdvals(A)), opnorm(A,1), opnorm(A,Inf)

# ╔═╡ 196d75ee-c304-4501-8ee3-b685d2fc2a65
md"""
# Skalarni produkt, norma i ortogonalnost funkcija


__Skalarni produkt__ na vektorskom prostoru $X$ je svako preslikavanje
$\cdot : X\times X \to \mathbb{R}$ sa sljedećim svojstvima:

1.  $x\cdot x\geq 0$
1.  $x\cdot x=0 \Leftrightarrow x=0$
2.  $x\cdot y=y\cdot x$
3.  $(\alpha x)\cdot y =\alpha (x\cdot y)$
3.  $(x+y)\cdot z=x\cdot z+y \cdot z$

Ukoliko je na vektorskom prostoru definiran skalarni produkt, normu možemo definirati kao

$$\|x\|=\sqrt{x\cdot x}.$$

Također, ako je $x \cdot y=0$ kažemo da su vektori $x$ i $y$ __međusobno ortogonalni (okomiti)__.

Na primjer, standardna vektorska norma

$$\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}$$

je definirana pomoću skalarnog produkta vektora,

$$x\cdot y=\sum_{i=1}^n  x_i y_i,$$

a vektori $x$ i $y$ su ortogonalni, odnosno $x\perp y$, ako je
$x\cdot y=0$.

Skalarni produkt funkcija definiramo pomoću određenog integrala:

$$f\cdot g = \int_a^b f(x)g(x) \, dx.$$

Ostale definicije ostaju iste:

$$\| f\|_2= \sqrt{f\cdot f} = \sqrt{\int_a^b [f(x)]^2 \, dx},$$

$$f\perp g \Longleftrightarrow f\cdot g =0.$$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╠═34d20879-1b40-498a-83fe-022ec52cdadc
# ╠═732d0a51-841e-4ca7-aca6-32fe38ef90f2
# ╟─54201b21-914a-4152-b187-3ecfe9982e65
# ╟─92d26690-7f76-4a73-ac98-38c42f9c46f4
# ╠═1ad65234-9ba8-4eaf-af4d-b881f8d5724d
# ╠═c1ffcca3-d537-4c4f-86b1-87215325c115
# ╟─31a749ca-adb3-4e0e-9484-4506e1586a74
# ╠═365a36f4-e332-40ea-9113-d322dc47d480
# ╠═10bc88f1-3675-4446-b31a-033848b6160e
# ╟─196d75ee-c304-4501-8ee3-b685d2fc2a65
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
