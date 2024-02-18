### A Pluto.jl notebook ###
# v0.19.36

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

## Kontrakcija

__Definicija.__ Funkcija $F:\mathbb{R}^n\to \mathbb{R}^n$ je __kontrakcija__ ako postoji broj $q<1$ za koji vrijedi

$$
\| F(x)-F(y)\| < q\|x-y\|\qquad \forall x,y.$$

## Teorem o fiksnoj točki

__Banachov teorem o fiksnoj točki.__
Ako je $F$ kontrakcija, onda niz definiran s

$$x_{k+1}=F(x_k)$$

konvergira prema jedinstvenom vektoru $\tilde x$ za kojeg vrijedi

$$
\tilde x = F(\tilde x).$$

Vektor $\tilde x$ se zove __fiksna točka__ funkcije $F$. Za pogrešku u $k$-tom koraku vrijede ocjene

$$
\|x_k- \tilde x\| \leq \frac{q}{1-q} \|x_k-x_{k-1}\|$$

i 

$$
\|x_k- \tilde x\| \leq \frac{q^k}{1-q} \|x_1-x_{0}\|,$$

pri čemu je druga ocjena bolja. Brzina konvergencije je __linearna__,

$$
\|x_{k+1}-\tilde x\| \leq q\| x_k-\tilde x\|.$$

__Napomena.__ Dokaz Banachovog teorema za jednostavniji slučaj $n=1$, dat ćemo prilikom razmatranja iterativne metoda za računanje nul-točaka funkcija.

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
function Jacobi(A::Array,b::Array,x::Array,tol::Float64=1000.0*eps())
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    d=1.0
    B=-(L+U)
    c=inv(D)*b
    q=opnorm(B,Inf)
    @show q
	println()
	y=Vector{Float64}(undef,n)
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

# ╔═╡ ac75bd58-59ca-4c44-a6a5-1208e568cc12
cond(A)

# ╔═╡ 2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# Provjerimo i normu relativnog reziduala
norm(r)/(norm(A)*norm(x))

# ╔═╡ adbd72cb-4dcf-490b-bbcb-1d681358c455
function GaussSeidel(A::Array,b::Array,x::Array,tol::Float64=1000.0*eps())
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    d=1.0
    B=-(I+L)\U
    c=(I+L)\(inv(D)*b)
	q=opnorm(B,Inf)
    @show q
	println()
    y=Vector{Float64}(undef,n)
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
	Random.seed!(3345)
	n₁=1024
	A₁=rand(n₁,n₁)+n₁*I
	b₁=rand(n₁)
	# Početni vektor
	x₁=rand(n₁)
end

# ╔═╡ 8794b613-ddab-4fa9-82e6-eec4192705dd
@time x₂,d₂=Jacobi(A₁,b₁,x₁);

# ╔═╡ f19969e6-8e49-4c00-987a-3d16425cd8d5
# Rezidual
norm(A₁*x₂-b₁)

# ╔═╡ 02d632a3-3ab8-4b02-ba94-709880df6313
@time A₁\b₁;

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
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "519c88b955a16a6f52e4beee9c744049f942c2fe"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═26644809-5e85-46be-a15b-dbe3f99b14ca
# ╠═daf29e3a-f99e-4df8-95f9-d5a6c71e1271
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
# ╠═ac75bd58-59ca-4c44-a6a5-1208e568cc12
# ╠═2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# ╠═adbd72cb-4dcf-490b-bbcb-1d681358c455
# ╠═9ae4a166-f3ee-4850-89ff-4c0a41bae48c
# ╠═4a832be5-26c0-46fa-bcbb-3afb4eb0cf2e
# ╟─92d8992a-f55e-4576-a1aa-c8b0de9e5806
# ╠═54e5e9c4-5ad2-45f9-a3ed-0aaced61663d
# ╠═8794b613-ddab-4fa9-82e6-eec4192705dd
# ╠═f19969e6-8e49-4c00-987a-3d16425cd8d5
# ╠═02d632a3-3ab8-4b02-ba94-709880df6313
# ╟─aba2f7f6-8690-4ede-8282-ad925e4aae8d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
