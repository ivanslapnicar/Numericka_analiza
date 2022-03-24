### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 26644809-5e85-46be-a15b-dbe3f99b14ca
using PlutoUI, Random, LinearAlgebra

# ╔═╡ daf29e3a-f99e-4df8-95f9-d5a6c71e1271
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 1a406352-1739-4709-85bc-6ca3ecb19253
md"""
# Iterativne metode


Za velike sustave, a posebno za sustave s malom ispunom (malo elemenata različitih od nule), te ukoliko je matrica sustava _strogo dijagonalno dominantna_ , rješenje se može brzo naći __iterativnim metodama__
(vidi [Numerička matematika, poglavlje 3.8](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)):

__Definicija.__ Funkcija $F:\mathbb{R}^n\to \mathbb{R}^n$ je __kontrakcija__ ako postoji broj $q<1$ za koji vrijedi

$$
\| F(x)-F(y)\| < q\|x-y\|\qquad \forall x,y.$$

__Banachov teorem o fiksnoj točki.__
Ako je $F$ kontrakcija, onda niz definiran s

$$x_{k+1}=F(x_k)$$

konvergira prema jedinstvenom vektoru $\tilde x$ za kojeg vrijedi

$$
\tilde x = F(\tilde x).$$

Broj $\tilde x$ se zove __fiksna točka__ funkcije $F$. Za pogrešku u $k$-tom koraku vrijede ocjene

$$
\|x_k- \tilde x\| \leq \frac{q}{1-q} \|x_k-x_{k-1}\|$$

i 

$$
\|x_k- \tilde x\| \leq \frac{q^k}{1-q} \|x_1-x_{0}\|,$$

pri čemu je druga ocjena bolja. Brzina konvergencije je __linearna__,

$$
\|x_{k+1}-\tilde x\| \leq q\| x_k-\tilde x\|.$$

"""

# ╔═╡ 610ef7a4-f0a6-42c8-a2cc-1a03cb155a22
md"""

# Jacobijeva i Gauss-Seidelova metoda

Neka je 

$$F(x)=Bx+c,$$

pri čemu je $B$ kvadratna matrica. Tada je

$$
\| F(x)-F(y)\|=\| Bx+c-(By+c)\|=\|B(x-y)\| \leq \|B\| \|x-y\|,$$

pa je $F$ kontrakcija ako je

$$
 \|B\|=q<1.$$

Neka je zadan sustav  $Ax=b$. Matricu $A$ rastavimo kao

$$
A=D\,(L+I+U)$$

pri čemu je $D$ dijagonalna matrica, $L$ strogo donje trokutasta matrica i $U$ strogo gornje trokutasta matrica.

## Jacobijeva metoda 

Neka je 

$$
B=-(L+U), \quad c=D^{-1}b.$$


Ako je matrica $A$ __strogo dijagonalno dominantna__, 

$$
\| B\|_{\infty} = \max_i \sum_{{j=1} \atop {j\neq i}}^n \frac{|a_{ij}|}{|a_{ii}|}<1$$

onda je preslikavanje $F$ kontrakcija (moguće je uzeti i druge norme) pa niz

$$
x_{k+1}=-(L+U)x_k+c$$

konvergira prema rješenju sustava $x$.

## Gauss-Seidelova metoda 

Neka je 

$$
B=-(I+L)^{-1}U, \quad c=(I+L)^{-1}\, D^{-1}b.$$

Bez dokaza navodimo sljedeću tvrdnju: ako je matrica $A$ strogo dijagonalno dominantna,
onda je preslikavanje $F$ kontrakcija pa niz

$$
x_{k+1}=-(I+L)^{-1}Ux_k+(I+L)^{-1}D^{-1}b,$$

odnosno

$$
x_{k+1}=-Lx_{k+1}-Ux_k+D^{-1}b,$$

konvergira prema rješenju sustava $x$.
"""

# ╔═╡ b2ed8cb0-1386-11eb-364f-ab9a25428f4b
md"
Pogledajmo kako izgleda rastav na faktore $A=D(L+I+U)$: 
"

# ╔═╡ 583e17b2-e189-49ae-8c80-0014d53c40c2
begin
	Random.seed!(123)
	n=5
	A=randn(n,n)
	# Napravimo matricu dijagonalno dominantnom
	A=A+n*I
	b=randn(n)
end

# ╔═╡ f1da40d0-1386-11eb-1102-698c5a0ac84a
A

# ╔═╡ 0489e230-1387-11eb-152e-2b64f289cfd7
D=Diagonal(A)

# ╔═╡ 69d74380-1387-11eb-127b-c5fcad190c08
inv(D)*A

# ╔═╡ 402ac6b0-1387-11eb-0672-337698b01aa9
L=inv(D)*tril(A,-1)

# ╔═╡ 664c0240-1379-11eb-2f93-1f3bb9532c19
U=inv(D)*triu(A,1)

# ╔═╡ 5aba7e24-0424-45f2-9716-3b32a71fc610
function Jacobi(A::Array,b::Array,x::Array)
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    tol=1000*eps()
    d=1.0
    B=-(L+U)
    c=inv(D)*b
    q=norm(B,Inf)
    # @show q
	println()
    while d>tol
        y=B*x+c
        d=norm(x-y,Inf)
        @show d
        x=y
    end
    x,d
end

# ╔═╡ 213d2b7b-b742-4274-9bb0-e029aec6f892
# Početni vektor
x₀=rand(n)

# ╔═╡ 91b52c67-de20-4bfc-9da6-3e04ed73b990
# x je rješenje, d je norma razlike dvije zadnje iteracije
x,d=Jacobi(A,b,x₀)

# ╔═╡ 0d07f057-9012-42ad-bf52-31a0f14614df
# Rezidual
r=A*x-b

# ╔═╡ 2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# Provjerimo i normu relativnog reziduala
norm(r)/(norm(A)*norm(x))

# ╔═╡ adbd72cb-4dcf-490b-bbcb-1d681358c455
function GaussSeidel(A::Array,b::Array,x::Array)
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    tol=1000*eps()
    d=1.0
    # B=-inv(I+L)*U
    B=-(I+L)\U
    c=(I+L)\(inv(D)*b)
    # @show norm(U,Inf)
    # y=Vector{Float64}(undef,n)
    while d>tol
        y=B*x+c
        d=norm(x-y)
		@show d
        x=y
    end
    x,d
end

# ╔═╡ 9ae4a166-f3ee-4850-89ff-4c0a41bae48c
xᵧ,dᵧ=GaussSeidel(A,b,x₀)

# ╔═╡ 4a832be5-26c0-46fa-bcbb-3afb4eb0cf2e
# Rezidual
A*xᵧ-b

# ╔═╡ 92d8992a-f55e-4576-a1aa-c8b0de9e5806
md"""
Izmjerimo brzinu za veće matrice:
"""

# ╔═╡ 54e5e9c4-5ad2-45f9-a3ed-0aaced61663d
begin
	n₁=1024
	A₁=rand(n₁,n₁)+n₁*I
	b₁=rand(n₁)
	# Početni vektor
	x₁=rand(n₁)
end

# ╔═╡ 8794b613-ddab-4fa9-82e6-eec4192705dd
@time GaussSeidel(A₁,b₁,x₁);

# ╔═╡ 02d632a3-3ab8-4b02-ba94-709880df6313
@time A\b;

# ╔═╡ aba2f7f6-8690-4ede-8282-ad925e4aae8d
md"""
__Zadatak.__ Analizirajte korištenje memorije i probajte preraditi programe tako da alociraju manje memorije.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.15"
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

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "f6532909bf3d40b308a0f360b6a0e626c0e263a8"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.1"

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
deps = ["Libdl", "libblastrampoline_jll"]
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

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "633f8a37c47982bff23461db0076a33787b17ecd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.15"

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

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═26644809-5e85-46be-a15b-dbe3f99b14ca
# ╟─daf29e3a-f99e-4df8-95f9-d5a6c71e1271
# ╟─1a406352-1739-4709-85bc-6ca3ecb19253
# ╟─610ef7a4-f0a6-42c8-a2cc-1a03cb155a22
# ╟─b2ed8cb0-1386-11eb-364f-ab9a25428f4b
# ╠═583e17b2-e189-49ae-8c80-0014d53c40c2
# ╠═f1da40d0-1386-11eb-1102-698c5a0ac84a
# ╠═0489e230-1387-11eb-152e-2b64f289cfd7
# ╠═69d74380-1387-11eb-127b-c5fcad190c08
# ╠═402ac6b0-1387-11eb-0672-337698b01aa9
# ╠═664c0240-1379-11eb-2f93-1f3bb9532c19
# ╠═5aba7e24-0424-45f2-9716-3b32a71fc610
# ╠═213d2b7b-b742-4274-9bb0-e029aec6f892
# ╠═91b52c67-de20-4bfc-9da6-3e04ed73b990
# ╠═0d07f057-9012-42ad-bf52-31a0f14614df
# ╠═2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# ╠═adbd72cb-4dcf-490b-bbcb-1d681358c455
# ╠═9ae4a166-f3ee-4850-89ff-4c0a41bae48c
# ╠═4a832be5-26c0-46fa-bcbb-3afb4eb0cf2e
# ╟─92d8992a-f55e-4576-a1aa-c8b0de9e5806
# ╠═54e5e9c4-5ad2-45f9-a3ed-0aaced61663d
# ╠═8794b613-ddab-4fa9-82e6-eec4192705dd
# ╠═02d632a3-3ab8-4b02-ba94-709880df6313
# ╟─aba2f7f6-8690-4ede-8282-ad925e4aae8d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
