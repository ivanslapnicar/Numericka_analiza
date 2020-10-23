### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 01708c02-1542-49de-b9ce-ef19b589d76e
begin
	using Polynomials
	using Plots
end

# ╔═╡ 68cc76f3-e124-4129-9339-98a21093bf1f
begin
	# Generirajmo slučajne točke
	using Random
	Random.seed!(125)
	n=6
	x=rand(n)
	y=rand(n)
	a=minimum(x)
	b=maximum(x)
end

# ╔═╡ 724b6f69-0f5c-4bcf-8363-42f308897070
md"""
# Interpolacijski polinomi

Neka je zadana $n+1$ točka

$$T_i=(x_i,y_i), \quad i=0,1,\ldots,n,\quad x_i\neq x_j.$$

## Standardna baza

Kroz zadane točke prolazi __interpolacijski polinom__ $p_n(x)$. Koeficijenti polinoma zadovoljavaju 
sustav linearnih jednadžbi $p_n(x_i)=y_i$, $i=0,\ldots,n$, odnosno

$$\begin{pmatrix} 
1 & x_0 & x_0^2 & x_0^3 & \cdots & x_0^n \cr
1 & x_1 & x_1^2 & x_1^3 & \cdots & x_1^n \cr
\vdots & & & & \vdots \cr
1 & x_n & x_n^2 & x_n^3 & \cdots & x_n^n \cr
\end{pmatrix}
\begin{pmatrix}a_0\cr a_1 \cr \vdots \cr a_n\end{pmatrix}
=\begin{pmatrix} y_0 \cr y_1 \cr \vdots \cr y_n\end{pmatrix}.$$

Matrica sustava $A$ se zove __Vandermondeova matrica__. Njena determinanta dana je formulom

$$\mathop{\mathrm{det}}(A)= \prod_{0\leq j<i\leq n}(x_i-x_j).$$

Kako su sve apscise različite ($x_i\neq x_j$ za $i\neq j$), vrijedi $\mathop{\mathrm{det}}(A)\neq 0$ pa je matrica $A$ regularna i zadani sustav ima jedinstveno rješenje - dakle, 

> interpolacijski polinom je __jedinstven__.

"""

# ╔═╡ 8cb77244-7845-4e64-920e-11a024f880e6
begin
	# Funkcije za manipulaciju s Vandermondeovim matricama
	import Base.getindex, Base.size
	struct Vandermonde{T} <: AbstractMatrix{T}
		c :: AbstractVector{T}
	end
	
	getindex(V::Vandermonde, i::Int, j::Int) = V.c[i]^(j-1)
	isassigned(V::Vandermonde, i::Int, j::Int) = isassigned(V.c, i)
	
	size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
	    throw(ArgumentError("Invalid dimension $r"))
	size(V::Vandermonde) = length(V.c), length(V.c)
	
	function Matrix(V::Vandermonde{T}) where T
		n=size(V, 1)
		M=Array{T}(undef,n, n)
		for i=1:n
			M[:,i] = V.c.^(i-1)
		end
		M
	end
	
end

# ╔═╡ 4f060d90-7b59-4447-ae59-856663acf3d7
A=Vandermonde(x)

# ╔═╡ d29a09b2-7cc6-477b-82fd-0591f1f0ec8f
begin
	# Vandermondeova matrica ima veliku kondiciju
	using LinearAlgebra
	cond(A)
end

# ╔═╡ 99d7f5e3-01c6-46e7-b16f-d42b46214738
c=A\y

# ╔═╡ 366d7a3e-ff6e-4bd4-a92d-935bd0d30c35
p=Polynomial(c)

# ╔═╡ f719cad8-cd2c-4013-9014-1caef03cc575
# Točke polinoma
scatter(x,y,label="Tocke")

# ╔═╡ f57091bd-a669-461f-ba5f-9ab7539a3c37
# Nacrtajmo polinom 
plot!(p,label=["Polinom"],xlims=(0,1),ylims=(-20,20))

# ╔═╡ 51668880-461a-45ce-82fe-204543677c75
begin
	# Nacrtajmo polinom s našom funkcijom
	xx=range(a,stop=b,length=100)
	pS=p.(xx)
	plot(xx,pS)
	scatter!(x,y)
end

# ╔═╡ ea329cb4-5b4b-49e3-825a-eb19d54a4e91
md"""
Za rješavanje zadanog sustava standardnim putem potrebno je $O(n^3)$ računskih operacija, no postoje metode kojima se Vandermondeovi sustavi mogu riješiti s $O(n^2)$ operacija.

Za izvrednjavanje polinoma u nekoj točki potrebno je $2n$ operacija (Hornerova shema).

Vandermondeove matrice uglavnom imaju veliku kondiciju pa ovaj način računanja koeficijenata polinoma može biti nestabilan.
Stoga se koriste i druge metode za računanje i izvredjavanje interpolacijskih polinoma.

## Lagrangeov interpolacijski polinom

Definirajmo $n+1$ polinom stupnja $n$:

$$
L_j(x)=\prod_{\displaystyle {i=0}\atop {\displaystyle i\neq j}}^n \frac{x-x_i}{x_j-x_i}.$$

Vrijedi 

$$
L_j(x_i)=\begin{cases}0, \quad i\neq j \\ 1,\quad i=j \end{cases}$$

pa je 

$$
p_n(x)=y_0\, L_0(x)+y_1 \, L_1(x)+\cdots + y_n\,  L_n(x).$$

Za računanje nazivnika polinoma prvi put je potrebno $O(n^2)$ operacija, ali se potom vrijednost 
$p_n(x)$ računa s $O(n)$ operacija (__objasnite kako!__). 

Navodimo implementaciju algoritma koja nije optimalno brza.
"""

# ╔═╡ 4195f039-de9e-4046-89c2-45e328b30478
L(t)=sum(y.*[prod(t .-x[[1:j-1;j+1:end]])/prod(x[j].-x[[1:j-1;j+1:end]]) 
        for j=1:n])

# ╔═╡ ef4cdc9e-7e9a-4b1a-b25c-c8e4266f8d6a
begin
	pL=Array{Float64}(undef,length(xx))
	for i=1:length(xx)
	    pL[i]=L(xx[i])
	end
end

# ╔═╡ 70cadf28-cb7c-4f44-836d-df7d7af666c2
begin
	plot(xx,pL)
	scatter!(x,y)
end

# ╔═╡ 0d08175e-43ab-41d4-b278-697154aa1966
norm(pS-pL,Inf)

# ╔═╡ f4766fb0-d630-4ef0-80e2-dedb057596bb
norm(abs.((pS-pL)./pL),Inf)

# ╔═╡ 65bf04af-9c28-4d16-b71d-b6d15a9371e3
md"""
## Newtonov interpolacijski polinom

Kod ovog polinoma koristi se baza

$$
1, x-x_0, (x-x_0)(x-x_1), (x-x_0)(x-x_1)(x-x_2),\ldots,(x-x_0)(x-x_1)\cdots (x-x_{n-1})$$

pa je interpolacijski polinom dan s

$$
p_n(x)=c_0 + c_1(x-x_0)+c_2(x-x_0)(x-x_1)+\cdots +c_n(x-x_0)(x-x_1)\cdots (x-x_{n-1}).$$

Koeficijenti interpolacijskog polinoma su rješenje __trokutastog__ sustava linearnih jednadžbi $Lc=y$,

$$
\begin{pmatrix} 
1 & 0 & 0 & 0 & \cdots & 0 \\
1 & x_1-x_0 & 0 & 0 & \cdots & 0 \\
1 & x_2-x_0 & (x_2-x_0)(x_2-x_1) & 0 & \cdots & 0 \\
\vdots & & & & \vdots \\
1 & x_n-x_0 & (x_n-x_0)(x_n-x_1) & (x_n-x_0)(x_n-x_1)(x_n-x_2) & \cdots & (x_n-x_0)\cdots (x_n-x_{n-1}) \\
\end{pmatrix}
\begin{pmatrix}c_0\\ c_1 \\ c_2 \\\vdots \\ a_n\end{pmatrix}
=\begin{pmatrix} y_0 \\ y_1 \\ y_2 \\ \vdots \\ y_n\end{pmatrix}.$$

Za formiranje donje trokutaste matrice $L$ potrebno je $O(n^2)$ operacija. Za računanje koeficijenata $c_0,\ldots,c_n$ potrebno je $O(n^2)$ operacija (rješavanje donje trokutastog sustava) i to rješenje je __stabilno__.

Za računanje $p_n(x)$ koristi se postupak koji je sličan Hornerovoj shemi. 
"""

# ╔═╡ d5755a6e-a155-4090-af0e-1076501397fa
# Računanje koeficijenata c
function mynewton(x,y)
    n=length(x)
    L=zeros(n,n)
    L[:,1]=ones(n)
    for i=2:n
        for j=2:i
            L[i,j]=prod([x[i]-x[k] for k=1:j-1])
        end
    end
    c=L\y
end  

# ╔═╡ 4976320f-d288-49a3-a1cf-44d2d223b488
cₙ=mynewton(x,y)

# ╔═╡ b85666e1-49c3-4f2d-b55b-9c48df27159a
# Računanje vrijednosti Newtonovog polinoma zadanog s točkama x i 
# koeficijentima c u točki t 
function evalnewton(c,x,t::Number)
    p=c[end]
    for i=length(c)-1:-1:1
        p=p*(t-x[i])+c[i]
    end
    p
end

# ╔═╡ 3f6a6282-3363-4a51-a5a8-bacc3d7e8880
begin
	pN=Array{Float64}(undef,length(xx))
	for i=1:length(xx)
	    pN[i]=evalnewton(cₙ,x,xx[i])
	end
end

# ╔═╡ 2060e3ad-c95f-4cd6-99c6-673da7dd64ac
begin
	plot(xx,pN)
	scatter!(x,y)
end

# ╔═╡ 0e0637af-4d77-442b-af67-ae3b7fe9075a
norm(abs.((pS-pN)./pN),Inf)

# ╔═╡ da067ec6-5085-467c-ae0d-143e1ba1e922
norm(abs.((pL-pN)./pN),Inf)

# ╔═╡ 477e1427-cd04-40f0-9b18-00ce68b1b837
md"""
Vidimo da su `pN` i `pL` bliže jedan drugome nego `pS` pa zaključujemo da su zaista točniji.
"""

# ╔═╡ 7433f2b7-7f03-4d70-9e7a-68df287cc14c


# ╔═╡ Cell order:
# ╟─724b6f69-0f5c-4bcf-8363-42f308897070
# ╠═01708c02-1542-49de-b9ce-ef19b589d76e
# ╠═68cc76f3-e124-4129-9339-98a21093bf1f
# ╠═8cb77244-7845-4e64-920e-11a024f880e6
# ╠═4f060d90-7b59-4447-ae59-856663acf3d7
# ╠═99d7f5e3-01c6-46e7-b16f-d42b46214738
# ╠═366d7a3e-ff6e-4bd4-a92d-935bd0d30c35
# ╠═f719cad8-cd2c-4013-9014-1caef03cc575
# ╠═f57091bd-a669-461f-ba5f-9ab7539a3c37
# ╠═51668880-461a-45ce-82fe-204543677c75
# ╠═d29a09b2-7cc6-477b-82fd-0591f1f0ec8f
# ╟─ea329cb4-5b4b-49e3-825a-eb19d54a4e91
# ╠═4195f039-de9e-4046-89c2-45e328b30478
# ╠═ef4cdc9e-7e9a-4b1a-b25c-c8e4266f8d6a
# ╠═70cadf28-cb7c-4f44-836d-df7d7af666c2
# ╠═0d08175e-43ab-41d4-b278-697154aa1966
# ╠═f4766fb0-d630-4ef0-80e2-dedb057596bb
# ╟─65bf04af-9c28-4d16-b71d-b6d15a9371e3
# ╠═d5755a6e-a155-4090-af0e-1076501397fa
# ╠═4976320f-d288-49a3-a1cf-44d2d223b488
# ╠═b85666e1-49c3-4f2d-b55b-9c48df27159a
# ╠═3f6a6282-3363-4a51-a5a8-bacc3d7e8880
# ╠═2060e3ad-c95f-4cd6-99c6-673da7dd64ac
# ╠═0e0637af-4d77-442b-af67-ae3b7fe9075a
# ╠═da067ec6-5085-467c-ae0d-143e1ba1e922
# ╟─477e1427-cd04-40f0-9b18-00ce68b1b837
# ╠═7433f2b7-7f03-4d70-9e7a-68df287cc14c
