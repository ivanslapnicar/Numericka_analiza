### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 5aba7e24-0424-45f2-9716-3b32a71fc610
begin
	using LinearAlgebra
	function myjacobi(A::Array,b::Array,x::Array)
	    D=Diagonal(A)
	    L=inv(D)*tril(A,-1)
	    U=inv(D)*triu(A,1)
	    tol=1000*eps()
	    d=1.0
	    B=-(L+U)
	    c=inv(D)*b
	    q=norm(B,Inf)
	    # @show q
	    while d>tol
	        y=B*x+c
	        d=norm(x-y,Inf)
	        # @show d
	        x=y
	    end
	    x,d
	end
end

# ╔═╡ 1a406352-1739-4709-85bc-6ca3ecb19253
md"""
# Iterativne metode


Za velike sustave, a posebno za sustave s malom ispunom (malo elemenata različitih od nule), te ukoliko je matrica sustava _strogo dijagonalno dominantna_ , rješenje se može brzo naći _iterativnim metodama_
(vidi [Numerička matematika, poglavlje 3.8](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)):

__Definicija.__ Funkcija $F:\mathbb{R}^n\to \mathbb{R}^n$ je _kontrakcija_ ako postoji broj $q<1$ za koji vrijedi

$$
\| F(x)-F(y)\| < q\|x-y\|\qquad \forall x,y.$$

__Banachov teorem o fiksnoj točki.__
Ako je $F$ kontrakcija, onda niz definiran s

$$x_{k+1}=F(x_k)$$

konvergira prema jedinstvenom vektoru $\tilde x$ za kojeg vrijedi

$$
\tilde x = F(\tilde x).$$

Broj $\tilde x$ se zove _fiksna točka_ funkcije $F$. Za pogrešku u $k$-tom koraku vrijede ocjene

$$
\|x_k- \tilde x\| \leq \frac{q}{1-q} \|x_k-x_{k-1}\|$$

i 

$$
\|x_k- \tilde x\| \leq \frac{q^k}{1-q} \|x_1-x_{0}\|,$$

pri čemu je druga ocjena bolja. Brzina konvergencije je _linearna_ ,

$$
\|x_{k+1}-\tilde x\| \leq q\| x_k-\tilde x\|.$$

"""

# ╔═╡ 610ef7a4-f0a6-42c8-a2cc-1a03cb155a22
md"""

## Jacobijeva i Gauss-Seidelova metoda

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

### Jacobijeva metoda 

Neka je 

$$
B=-(L+U), \quad c=D^{-1}b.$$


Ako je matrica $A$ _strogo dijagonalno dominantna_, 

$$
\| B\|_{\infty} = \max_i \sum_{{j=1} \atop {j\neq i}}^n \frac{|a_{ij}|}{|a_{ii}|}<1$$

onda je preslikavanje $F$ kontrakcija (moguće je uzeti i druge norme) pa niz

$$
x_{k+1}=-(L+U)x_k+c$$

konvergira prema rješenju sustava $x$.

### Gauss-Seidel-ova metoda 

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
	import Random
	Random.seed!(123)
	n=8
	A=rand(n,n)
	# Napravimo matricu dijagonalno dominantnom
	A=A+n*I
	b=rand(n)
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

# ╔═╡ 213d2b7b-b742-4274-9bb0-e029aec6f892
# Početni vektor
x₀=rand(n)

# ╔═╡ 91b52c67-de20-4bfc-9da6-3e04ed73b990
# x je rješenje, d je norma razlike dvije zadnje iteracije
x,d=myjacobi(A,b,x₀)

# ╔═╡ 0d07f057-9012-42ad-bf52-31a0f14614df
# Rezidual
r=A*x-b

# ╔═╡ 2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# Provjerimo i normu relativnog reziduala
norm(r)/(norm(A)*norm(x))

# ╔═╡ adbd72cb-4dcf-490b-bbcb-1d681358c455
function mygaussseidel(A::Array,b::Array,x::Array)
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    tol=1000*eps()
    d=1.0
    # B=-inv(I+L)*U
    B=-(I+L)\U
    c=(I+L)\(inv(D)*b)
    # @show norm(U,Inf)
    y=Vector{Float64}(undef,n)
    while d>tol
        y=B*x+c
        d=norm(x-y)
		# @show d
        x=y
    end
    x,d
end

# ╔═╡ 9ae4a166-f3ee-4850-89ff-4c0a41bae48c
xᵧ,dᵧ=mygaussseidel(A,b,x₀)

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
@time mygaussseidel(A₁,b₁,x₁);

# ╔═╡ 02d632a3-3ab8-4b02-ba94-709880df6313
@time A\b;

# ╔═╡ aba2f7f6-8690-4ede-8282-ad925e4aae8d
md"""
__Zadatak.__ Probajte preraditi programe tako da alociraju manje memorije.
"""

# ╔═╡ 76473de1-3b58-4162-a442-eb9f189d111a


# ╔═╡ Cell order:
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
# ╠═76473de1-3b58-4162-a442-eb9f189d111a
