### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 881ed68b-ac65-434a-9d62-df963fb033b1
begin
	using LinearAlgebra
	using Plots
	import Random
end

# ╔═╡ e0f03ec0-f136-411f-a270-bcc85fde0579
md"""
# Primjene QR rastava


Ortogonalne matrice imaju dva važna svojstva:

$Q^{-1}=Q^T$

$\|Qx\|_2=\|x\|_2,\quad \forall x$

Prvo svojstvo slijedi iz definicije ortogonalne matrice, jer je $Q^TQ=I$, 
a drugo svojstvo slijedi iz

$$\|Qx\|_2^2=(Qx)^T(Qx)=x^TQ^TQx=x^Tx=\|x\|_2^2.$$

## Sustav linearnih jednadžbi

QR rastav možemo koristiti za rješavanje sustava linearnih jednadžbi $Ax=b$:
množenjem jednakosti $QRx=b$ s lijeva s $Q^T$ dobijemo $Rx=Q^Tb$ pa preostaje riješiti trokutasti sustav. 

U odnosu na rješenje pomoću Gaussove eliminacije vrijedi:

* broj računskih operacija se udvostruči,
* rješenje je nešto točnije, i
* nema rasta elemenata (pivotiranje nije potrebno).
"""

# ╔═╡ d0aa1de2-f2e5-4c59-b7ce-a0022bf6dbef
begin
	Random.seed!(123)
	n=10
	A=randn(n,n)
	b=randn(n)
	Q,R=qr(A)
	c=transpose(Q)*b
	# Trokutasti sustav
	x=R\c
end

# ╔═╡ fbfa49ba-8c53-4fdf-ae26-d9b678fcce2e
# Rezidual
norm(A*x-b)

# ╔═╡ 45923270-4d63-4eaf-8d99-f9e9baab558e
md"""
## Problem najmanjih kvadrata

Programi za rješavanje problema najmanjih kvadrata uglavnom koriste QR rastav. Vrijedi

$$
\|Ax-b\|^2_2=\|QRx-b\|_2^2=\|Q(Rx-Q^Tb)\|_2=\|Rx-Q^Tb\|_2^2.$$

Neka je

$$
R=\begin{bmatrix}R_0 \\ 0\end{bmatrix},\quad Q^Tb =\begin{bmatrix}c\\ d \end{bmatrix}.$$

Tada je

$$
\|Rx-Q^T b\|_2^2 = \| R_0x-c\|_2^2+\|d\|_2^2$$

pa je rješenje trokutastog sustava

$$
R_0x=c$$
    
rješenje problema najmanjih kvadrata. Riješimo zadatak iz bilježnice o regresiji:
"""

# ╔═╡ 2e4f479e-9a3a-4495-abb9-52a3203a35c4
begin
	y₁=[1,3,2,4,3]
	A₁=transpose([1 2 3 6 7;1 1 1 1 1])
end

# ╔═╡ 743b75c7-d6a0-4fa0-bd99-6412378453e2
#?qr  # Pogledajmo strukturu rješenja

# ╔═╡ 04bb12df-5d6a-460f-a381-23058dcfe0b7
begin
	F₁=qr(A₁)
	c₁=transpose(Matrix(F₁.Q))*y₁
	x₁=F₁.R\c₁
end

# ╔═╡ c4762bb2-56fd-4e2a-9255-4922f41addfe
# Ugrađena funkcija
A₁\y₁

# ╔═╡ 881daa94-e922-4b76-bca3-6b255049ba67
begin
	# Veći primjer
	m₂=8
	n₂=5
	A₂=randn(m₂,n₂)
	b₂=randn(m₂)
	F₂=qr(A₂)
end

# ╔═╡ 16461985-6cb3-4ddb-9eab-bd13917b8b69
F₂.R

# ╔═╡ d01748e8-169d-4d39-bf79-90c0c0d2df56
# Spremanje generatora
F₂.factors

# ╔═╡ 0a0e0ea0-92c8-11eb-21d5-f3bb0e0946f0
A₂

# ╔═╡ f4b9b30e-92c7-11eb-0c68-a91a2bd974ff
begin
	a=A₂[:,1]
	v=copy(a)
	v[1]+=norm(a)
	v./=v[1]
end

# ╔═╡ ab0d603e-4e1f-4ec9-a247-c9c78c9e6488
begin
	# Rješenje
	c₂=transpose(Matrix(F₂.Q))*b₂
	x₂=F₂.R\c₂
end

# ╔═╡ f97a14a2-1e13-4a2d-81a6-c848666dbbfe
# Ugrađena funkcija
A₂\b₂

# ╔═╡ 0ad1e3dc-e5d4-416a-b3e0-21fd710ed241
md"""
## Numerička  "ortogonalizacija" polinoma

Numerička ortogonalizacija potencija vektora daje ortogonalne polinome.
"""

# ╔═╡ 24ac4e43-b99d-439e-aa6b-ba3bc494b8dd
begin
	# Standardna baza
	xₒ=range(-1,stop=1,length=101)
	# Kvazi Vandermondeova matrica
	V=[xₒ.^0 xₒ.^1 xₒ.^2 xₒ.^3 xₒ.^4 xₒ.^5]
end

# ╔═╡ 9006213a-325f-4296-a0ec-d5c02b4ec5e8
plot(xₒ,V,title="Standardna baza",legend=:bottomright,
	label=["1" "x" "x^2" "x^3" "x^4" "x^5"])

# ╔═╡ f157bb76-8445-4994-b08e-260b2a51ad51
begin
	# Ortogonalizacija s težinskom funkcijom ω(x)=1 daje normirane Legendreove polinome.
	Fₒ=qr(V)
	Qₒ=Matrix(Fₒ.Q)*sign.(Diagonal(Fₒ.R))
end

# ╔═╡ 31fb59d0-92c9-11eb-2aef-19e75eab1844
plot(xₒ,Qₒ,title="Legendreovi polinomi",label=["L₀" "L₁" "L₂" "L₃" "L₄" "L₅"])

# ╔═╡ 4b343213-f8d2-45d3-bae8-999dee675089
md"""
Dobiveni normirani vektori su vrijednosti skaliranih Legendreovih polinoma iz blježnice [NA12 Ortogonalni polinomi.ipynb](NA12%20Ortogonalni%20polinomi.ipynb).

Da bi dobili Čebiševljeve polinome, trebamo dodati težinsku funkciju $\omega(x)=\displaystyle\frac{1}{\sqrt{1-x^2}}$ i preraditi funkciju `GramSchmidtQR()` iz bilježnice [NA15 QR rastav.ipynb](NA15%20QR%20rastav.ipynb) tako da računa __težinske skalarne produkte__. Dobiveni normirani vektori su vrijednosti skaliranih Čebiševljevih polinoma.
"""

# ╔═╡ 830acec9-74f1-4b02-a5f7-b65b05f457ef
function WeightedGramSchmidtQR(A::Array,ω::Vector)
    m,n=size(A)
    R=zeros(n,n)
    Q=Array{Float64,2}(undef,m,n)
    R[1,1]=norm(A[:,1])
    Q[:,1]=A[:,1]/R[1,1]
    for k=2:n
        for i=1:k-1
            R[i,k]=Q[:,i]⋅(A[:,k].*ω)/(Q[:,i]⋅(Q[:,i].*ω))
        end
        t=A[:,k]-sum([R[i,k]*Q[:,i] for i=1:k-1])
        R[k,k]=norm(t)
        Q[:,k]=t/R[k,k]
    end
    return Q,R
end

# ╔═╡ e0bf68c5-1e30-4621-871b-47eb71bfd2bb
begin
	x₃=range(-0.99,stop=0.99,length=199)
	ω=1 ./(sqrt.(1.0.-x₃.^2))
	# Kvazi Vandermonde-ova matrica
	V₃=[x₃.^0 x₃.^1 x₃.^2 x₃.^3 x₃.^4 x₃.^5]
end

# ╔═╡ be2197e5-565e-4b13-acb9-a245f2f3780e
begin
	Q₃,R₃=WeightedGramSchmidtQR(V₃,ω)
	Q₃=Q₃*sign.(Diagonal(R₃))
end

# ╔═╡ 957a675e-a12a-47bf-b196-bd4f3431849a
plot(x₃,Q₃,title="Čebiševljevi polinomi",label=["T₀" "T₁" "T₂" "T₃" "T₄" "T₅"])

# ╔═╡ 656aa3d6-0ebc-4ca5-95d5-071af0f5c2c0
md"""
__Zadatak.__ Normirajte stupce matrice $Q₃$ tako da vektori
poprimaju sve vrijednosti u intervalu $[-1,1]$. 
"""

# ╔═╡ Cell order:
# ╟─e0f03ec0-f136-411f-a270-bcc85fde0579
# ╠═881ed68b-ac65-434a-9d62-df963fb033b1
# ╠═d0aa1de2-f2e5-4c59-b7ce-a0022bf6dbef
# ╠═fbfa49ba-8c53-4fdf-ae26-d9b678fcce2e
# ╟─45923270-4d63-4eaf-8d99-f9e9baab558e
# ╠═2e4f479e-9a3a-4495-abb9-52a3203a35c4
# ╠═743b75c7-d6a0-4fa0-bd99-6412378453e2
# ╠═04bb12df-5d6a-460f-a381-23058dcfe0b7
# ╠═c4762bb2-56fd-4e2a-9255-4922f41addfe
# ╠═881daa94-e922-4b76-bca3-6b255049ba67
# ╠═16461985-6cb3-4ddb-9eab-bd13917b8b69
# ╠═d01748e8-169d-4d39-bf79-90c0c0d2df56
# ╠═0a0e0ea0-92c8-11eb-21d5-f3bb0e0946f0
# ╠═f4b9b30e-92c7-11eb-0c68-a91a2bd974ff
# ╠═ab0d603e-4e1f-4ec9-a247-c9c78c9e6488
# ╠═f97a14a2-1e13-4a2d-81a6-c848666dbbfe
# ╟─0ad1e3dc-e5d4-416a-b3e0-21fd710ed241
# ╠═24ac4e43-b99d-439e-aa6b-ba3bc494b8dd
# ╠═9006213a-325f-4296-a0ec-d5c02b4ec5e8
# ╠═f157bb76-8445-4994-b08e-260b2a51ad51
# ╠═31fb59d0-92c9-11eb-2aef-19e75eab1844
# ╟─4b343213-f8d2-45d3-bae8-999dee675089
# ╠═830acec9-74f1-4b02-a5f7-b65b05f457ef
# ╠═e0bf68c5-1e30-4621-871b-47eb71bfd2bb
# ╠═be2197e5-565e-4b13-acb9-a245f2f3780e
# ╠═957a675e-a12a-47bf-b196-bd4f3431849a
# ╟─656aa3d6-0ebc-4ca5-95d5-071af0f5c2c0
