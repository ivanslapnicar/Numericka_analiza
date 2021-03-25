### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ aa7480bd-e15a-4f07-bcea-f3f2c9f3f5fb
using LinearAlgebra

# ╔═╡ 6075a4cb-5931-49d1-987c-ffc0f40ebb12
using DelimitedFiles

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

## PageRank

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
begin
	using SparseArrays
	G=sparse(i,j,1.0)
end

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
## Ideja

Započnimo slučajnu šetnju iz vektora $x_0=\begin{bmatrix} 1/n \cr 1/n \cr \vdots \cr 1/n \end{bmatrix}$.

Sljedeći vektori su  

$$\begin{aligned}
x_1&=A\cdot x_0 \cr
x_2&=A\cdot x_1 \cr
x_3&=A\cdot x_2\cr
& \ \vdots
\end{aligned}$$

Preslikavanje $A(x)=Ax$ nije kontrakcija u smislu Banachovog teorema o fiksnoj točci jer je $\|A\|_1=1$, ali se može pokazati da ima fiksnu točku. Također, ako je $x\geq 0$ (po komponentama), onda je $\|Ax\|_1=\|x\|_1$.

Kada se vektor _stabilizira_:

$$
A\cdot x\approx x,$$

tada je $x[i]$ _rang stranice_ $i$.
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
W=readdlm("../files/web-Stanford.txt",Int,comments=true)

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

# ╔═╡ 2188e21d-f3bb-41f3-8a9c-c425ba4f0887


# ╔═╡ Cell order:
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
# ╠═aa7480bd-e15a-4f07-bcea-f3f2c9f3f5fb
# ╠═fece3020-0f09-11eb-0f69-237286bd58af
# ╠═739c238c-03db-4ee6-9fb7-f8e5b93282f8
# ╠═6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
# ╠═cb04da5e-0f08-11eb-21b9-8fdaea539145
# ╠═870c8ccd-7e2c-489c-957f-fc34651bb65f
# ╠═738462b0-62a9-4aed-8ca7-687fb51d52e2
# ╠═2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
# ╠═023e6d22-6505-4211-8fec-55ae732405bc
# ╠═fae1bfcb-ef52-4cd8-a066-cf138c8697f8
# ╟─5a02f8ad-3f97-4201-b903-9ed789721f81
# ╠═6075a4cb-5931-49d1-987c-ffc0f40ebb12
# ╠═47394960-0f02-11eb-1ddf-cb6b81f096b4
# ╠═86ca6760-1382-11eb-1e44-0ff4051234a0
# ╠═573a625a-ad1c-4133-bae3-342a7501b492
# ╠═0b9a7940-8d5d-11eb-372b-b7a6c4fc2596
# ╠═83fdc63f-4aac-45c0-a226-87a4830f697e
# ╠═b473be54-b1a2-4f34-80d8-742386c9535b
# ╠═d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
# ╠═12eb4b3f-6fe1-4bba-9556-4378eab6e191
# ╠═a445d5e2-adef-4d13-b0b3-0af37f7039d6
# ╠═2188e21d-f3bb-41f3-8a9c-c425ba4f0887
