### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 4eb20c31-d9d7-41c4-886e-e1200c6b62b2
using PlutoUI, SparseArrays, LinearAlgebra, DelimitedFiles

# ╔═╡ 65f12579-3c3f-452b-852a-c26d45ec595f
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ c1dc1294-a3a7-431f-8ecf-ea3c8448c44a
md"""
# PageRank
"""

# ╔═╡ 85992dca-d40f-461a-8984-ac57ddfff970
md"""
## Doba pretraživanja

google (i ostali)


* [50 milijardi stranica](http://www.worldwidewebsize.com/), [3.5 milijardi pretraga dnevno](http://www.internetlivestats.com/google-search-statistics/)
* __PageRank__
* povijest, kontekst - kolačići, spremanja podataka (o Vama), [200+ parametara](http://backlinko.com/google-ranking-factors)

Nakon doba slijedilo je __doba preporuka__, a od prošle godine imamo i __doba AI__ (ima li bolje ime?). 

 __age of search__ $\to$ __age of recommendations__ $\to$ __age of AI__

# Matrica prijelaza i graf

* Teorija grafova i linearna algebra
* [C. Moler, Google PageRank](https://www.mathworks.com/moler/exm/chapters/pagerank.pdf)


Neki programi:

* [https://github.com/purzelrakete/Pagerank.jl](https://github.com/purzelrakete/Pagerank.jl)

* [https://gist.github.com/domluna/2b9358ccc89fee7d5e26](https://gist.github.com/domluna/2b9358ccc89fee7d5e26)

Probat ćemo primjer iz Molerovog članka.
"""

# ╔═╡ d761f999-c402-49c4-bb4f-b52c23475db1
begin
	i = vec([ 2 6 3 4 4 5 6 1 1])
	j = vec([ 1 1 2 2 3 3 3 4 6])
end

# ╔═╡ daf1a60f-9827-418b-adbe-06d2046b1c96
G=sparse(i,j,1.0)

# ╔═╡ f874ba0c-9d9e-437b-9709-723ce34756c3
Matrix(G)

# ╔═╡ f9abf219-8880-4d45-a556-ba20ba114a56
begin
	G₁=similar(G)
	c=sum(G,dims=1)
	n=size(G,1)
	for j=1:n
	    if c[j]>0
	        G₁[:,j]=G[:,j]/c[j]
	    end
	end
end

# ╔═╡ 726e69a4-1606-466f-a927-2effba6bcaf2
Matrix(G₁)

# ╔═╡ 7c96e1d5-c076-4638-8c95-fe5bc8ab8936
md"""
* vjerojatnost da pratimo neki link je $p$
* vjerojatnost da posjetimo slučajnu stranicu je $1-p$
* google koristi $p=0.85$ ?
"""

# ╔═╡ ae67506f-65b2-4b67-af2e-3fc07f5044e9
begin
	p=0.85
	δ=(1-p)/n
end

# ╔═╡ ef68750f-be24-4d53-9d60-9715d4064ecf
z = ((1-p)*(c.!=0) + (c.==0))/n

# ╔═╡ 420af1cd-becc-4005-8534-96a6aaddbde5
A=p*G₁+ones(n)*z

# ╔═╡ 9d3b7528-4035-4304-bee9-a9407bced36f
sum(A,dims=1)

# ╔═╡ 90edb7b5-c882-41a6-a48b-ba15373f2283
md"""
# Slučajna šetnja i stabilno stanje

Započnimo slučajnu šetnju iz vektora $x_0=\begin{bmatrix} 1/n \cr 1/n \cr \vdots \cr 1/n \end{bmatrix}$.

Sljedeći vektori su  

$$\begin{aligned}
x_1&=A\cdot x_0 \cr
x_2&=A\cdot x_1 \cr
x_3&=A\cdot x_2\cr
& \ \vdots
\end{aligned}$$

Preslikavanje $A(x)=Ax$ nije kontrakcija u smislu Banachovog teorema o fiksnoj točci jer je $\|A\|_1=1$, ali se može pokazati da ima fiksnu točku. Također, ako je $x\geq 0$ (po komponentama), onda je $\|Ax\|_1=\|x\|_1$.

Kada se vektor __stabilizira__:

$$
A\cdot x\approx x,$$

tada je $x[i]$ __rang stranice__ $i$.
"""

# ╔═╡ fece3020-0f09-11eb-0f69-237286bd58af
function PageRank(G₁::SparseMatrixCSC{Float64,Int64},steps::Int)
	G=copy(G₁)
	p=0.85
	c=sum(G,dims=1)/p
	n=size(G,1)
	for i=1:n
	    G.nzval[G.colptr[i]:G.colptr[i+1]-1]./=c[i]
	end
	e=ones(n)
	x=e/n
	z = vec(((1-p)*(c.!=0) + (c.==0))/n)
	for j=1:steps
	    x=G*x.+(z⋅x)
	end
	return x
end

# ╔═╡ 9ddf938f-7313-4b69-9d48-98c1d342cd27
md"
__Potrebno je razumijeti CSC format i njegovo korištenje.__
"

# ╔═╡ 739c238c-03db-4ee6-9fb7-f8e5b93282f8
fieldnames(typeof(G))

# ╔═╡ 6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
G

# ╔═╡ cb04da5e-0f08-11eb-21b9-8fdaea539145
Matrix(G)

# ╔═╡ 870c8ccd-7e2c-489c-957f-fc34651bb65f
G.colptr

# ╔═╡ 738462b0-62a9-4aed-8ca7-687fb51d52e2
G.nzval

# ╔═╡ 2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
G.rowval

# ╔═╡ 023e6d22-6505-4211-8fec-55ae732405bc
# Početni vektor
x=ones(n)/n

# ╔═╡ fae1bfcb-ef52-4cd8-a066-cf138c8697f8
PageRank(G,15)

# ╔═╡ 5a02f8ad-3f97-4201-b903-9ed789721f81
md"""
## [Stanford web graph](http://snap.stanford.edu/data/web-Stanford.html)

Nešto veći testni problem.
"""

# ╔═╡ 47394960-0f02-11eb-1ddf-cb6b81f096b4
W=readdlm("./files/web-Stanford.txt",Int,comments=true)

# ╔═╡ 86ca6760-1382-11eb-1e44-0ff4051234a0
#?sparse

# ╔═╡ 573a625a-ad1c-4133-bae3-342a7501b492
S=sparse(W[:,2],W[:,1],1.0)

# ╔═╡ 0b9a7940-8d5d-11eb-372b-b7a6c4fc2596
length(S.nzval)

# ╔═╡ 83fdc63f-4aac-45c0-a226-87a4830f697e
@time x100=PageRank(S,100)

# ╔═╡ b473be54-b1a2-4f34-80d8-742386c9535b
x101=PageRank(S,101);

# ╔═╡ d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
maximum(abs,(x101-x100)./x101)

# ╔═╡ 12eb4b3f-6fe1-4bba-9556-4378eab6e191
# Ranks
sort(x100,rev=true)

# ╔═╡ a445d5e2-adef-4d13-b0b3-0af37f7039d6
# Pages
sortperm(x100,rev=true)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
DelimitedFiles = "~1.9.1"
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "935e84cd056970124a240dbef905262476a913d0"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

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
# ╠═4eb20c31-d9d7-41c4-886e-e1200c6b62b2
# ╠═65f12579-3c3f-452b-852a-c26d45ec595f
# ╟─c1dc1294-a3a7-431f-8ecf-ea3c8448c44a
# ╟─85992dca-d40f-461a-8984-ac57ddfff970
# ╠═d761f999-c402-49c4-bb4f-b52c23475db1
# ╠═daf1a60f-9827-418b-adbe-06d2046b1c96
# ╠═f874ba0c-9d9e-437b-9709-723ce34756c3
# ╠═f9abf219-8880-4d45-a556-ba20ba114a56
# ╠═726e69a4-1606-466f-a927-2effba6bcaf2
# ╟─7c96e1d5-c076-4638-8c95-fe5bc8ab8936
# ╠═ae67506f-65b2-4b67-af2e-3fc07f5044e9
# ╠═ef68750f-be24-4d53-9d60-9715d4064ecf
# ╠═420af1cd-becc-4005-8534-96a6aaddbde5
# ╠═9d3b7528-4035-4304-bee9-a9407bced36f
# ╟─90edb7b5-c882-41a6-a48b-ba15373f2283
# ╠═fece3020-0f09-11eb-0f69-237286bd58af
# ╟─9ddf938f-7313-4b69-9d48-98c1d342cd27
# ╠═739c238c-03db-4ee6-9fb7-f8e5b93282f8
# ╠═6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
# ╠═cb04da5e-0f08-11eb-21b9-8fdaea539145
# ╠═870c8ccd-7e2c-489c-957f-fc34651bb65f
# ╠═738462b0-62a9-4aed-8ca7-687fb51d52e2
# ╠═2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
# ╠═023e6d22-6505-4211-8fec-55ae732405bc
# ╠═fae1bfcb-ef52-4cd8-a066-cf138c8697f8
# ╟─5a02f8ad-3f97-4201-b903-9ed789721f81
# ╠═47394960-0f02-11eb-1ddf-cb6b81f096b4
# ╠═86ca6760-1382-11eb-1e44-0ff4051234a0
# ╠═573a625a-ad1c-4133-bae3-342a7501b492
# ╠═0b9a7940-8d5d-11eb-372b-b7a6c4fc2596
# ╠═83fdc63f-4aac-45c0-a226-87a4830f697e
# ╠═b473be54-b1a2-4f34-80d8-742386c9535b
# ╠═d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
# ╠═12eb4b3f-6fe1-4bba-9556-4378eab6e191
# ╠═a445d5e2-adef-4d13-b0b3-0af37f7039d6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
