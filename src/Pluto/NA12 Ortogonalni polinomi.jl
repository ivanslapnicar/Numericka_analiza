### A Pluto.jl notebook ###
# v0.12.8

using Markdown
using InteractiveUtils

# ╔═╡ 934a0191-a567-497f-b8e5-61ca5ca63a30
begin
	using Polynomials
	using Plots
	using SymPy
end

# ╔═╡ fdd0be01-c2a2-42c2-9a9b-1fcca24303f9
md"""
# Ortogonalni polinomi

Neka je 

$$
L(x_0,x_1,\ldots,x_n)$$

(pot)prostor razapet linearno nezavisnim vektorima (ili funkcijama) $x_0,x_1,\ldots,x_n$.

Radi se o skupu svih linearnih kombinacija zadanih vektora. 

Koristeći __Gram-Schmidtov postupak ortogonalizacije__ možemo izračunati __ortogonalnu bazu__ tog (pot)prostora 

$$
y_0,y_1,\ldots,y_n, $$

za koju vrijedi

$$
(y_i,y_j)=0,\quad i\neq j. \tag{1}$$

Neka je 

$$\begin{aligned}
y_0&=x_0\\
y_1&=x_1-\frac{(x_1,y_0)}{(y_0,y_0)}y_0\\
y_2&=x_2-\frac{(x_2,y_0)}{(y_0,y_0)}y_0-\frac{(x_2,y_1)}{(y_1,y_1)}y_1\\
& \vdots \\
y_n&=x_n-\sum_{j=0}^{n-1} \frac{(x_n,y_j)}{(y_j,y_j)}y_j.
\end{aligned}$$

Svaki $y_j$ je linearna kombinacija od $x_0,x_1,\ldots,x_j$ pa su $y_j$ linearno nezavisni i vrijedi

$$
L(x_0,x_1,\ldots,x_n)=L(y_0,y_1,\ldots,y_n).$$

Direktnom provjerom se vidi da vrijedi (1).

__Težinski skalarni produkt__  funkcija $f$ i $g$ na intervalu $[a,b]$ s težinom $\omega(x)>0$ je

$$
(f,g)_\omega=\int_a^b f(x)g(x)\omega(x)\, dx$$

Funkcije $f$ i $g$ su _ortogonalne_ ako je $(f,g)_\omega=0$.

__Ortogonalni polinomi__ nastaju ortogonalizacijom polinoma

$$
1,x,x^2,x^3,\ldots,x^n. \tag{2}$$

Različiti odabiri težinske funkcije daju različite sustave ortogonalnih polinoma.  
"""

# ╔═╡ a9eac96e-8ca5-464f-83d1-3ed9300e9d43
md"""
## Legendreovi polinomi

Ortogonalizirajmo sustav (2) uz

$$
[a,b]=[-1,1], \quad \omega(x)=1,$$

koristeći paket `SymPy.jl` za simboličko računanje. 
"""

# ╔═╡ b999b693-8366-4d4a-8674-c080d250be3f
begin
	a=-1
	b=1
	n=8
	P=Array{Any,1}(undef,n)
	x=Sym("x")
	P[1]=x^0
	ω(x)=1
	for k=2:n
	    P[k]=x^(k-1)
	    for j=1:k-1
	        P[k]=P[k]-SymPy.integrate(x->x^(k-
					1)*P[j]*ω(x),a,b)*P[j]/SymPy.integrate(x->P[j]*P[j]*ω(x),a,b)
	    end
	end
end

# ╔═╡ bab6e737-737c-4c33-be1f-6ba60d45649f
md"""
Julia indeksiranje započima s 1 pa su svi indeksi pomaknuti, odnosno

$$
P_0(x)=P[1], \ P_1(x)=P[2], \ldots$$
"""

# ╔═╡ 8cb49fec-a93a-4456-b159-fbe178b0cb67
P[1]

# ╔═╡ 831cc7ee-31fb-42c5-be87-a75f971215fa
P[4]

# ╔═╡ 3b90d2fd-afa0-466b-8fa5-e9748608ca00
P[6]

# ╔═╡ c240885a-008a-4ba3-8bae-b7dfaf1d1591
P[7]

# ╔═╡ d2dffde7-bee4-4b06-bd6d-a548b5dc2d5b
P[8]

# ╔═╡ 6957e729-7be5-44e9-ae1c-3087bcb10a37
md"""
Polinomi $P_n$ su do na množenje konstantom jednaki __Legendreovim polinomima__

$$
L_n(x)=\frac{1}{2^n n!}\frac{d^n}{dx^n}(x^2-1)^n, \quad n=0,1,2,3,\ldots$$
"""

# ╔═╡ 0abd5a32-7398-4f1d-aa60-46e204713a75
begin
	L=Array{Any,1}(undef,n)
	L[1]=x^0
	for k=1:n-1
	    L[k+1]=expand(diff((x^2-1)^k/(2^k*factorial(k)),x,k))
	end
end

# ╔═╡ 85a8addc-73bc-4e56-a8a9-3137ab443d15
L[1], P[1]

# ╔═╡ ea0adfdc-373f-40f7-a256-08d6424b6165
L[2],P[2]

# ╔═╡ e0b11a5c-becd-48b2-9238-f8dd8caed9ca
L[4],P[4]

# ╔═╡ c494410b-62dd-4822-8413-1e159b10d5c7
L[7]

# ╔═╡ 2ac8adab-18a5-4b4e-bcb3-d4e777fd9c55
P[7]

# ╔═╡ 11f09cbf-a386-43fb-9b02-1a0303ee4e0f
L[7]*16/231

# ╔═╡ 2e40d6a7-a59c-4342-9246-1c347e715795
md"""
Pored ortogonalnosti, vrijede sljedeća svojstva:

*  $L_n(x)$ ima $n$ različitih nul-točaka na intervalu $[-1,1]$, 
* vrijedi __tročlana rekurzivna formula__: 

$$L_{n+1}(x)=\frac{2n+1}{n+1}\,x\, L_n(x)-\frac{n}{n+1} L_{n-1}(x).$$

Izračunajmo polinome numerički i nacrtajmo ih:
"""

# ╔═╡ 83b09ece-0e8c-40c9-83a0-2204bba04323
begin
	n₁=40
	Lₙ=Array{Any,1}(undef,n₁)
	Lₙ[1]=Polynomial([1])
	Lₙ[2]=Polynomial([0,1])
	for i=3:n₁
	    Lₙ[i]=(2*i-3)*Lₙ[2]*Lₙ[i-1]/(i-1)-(i-2)*Lₙ[i-2]/(i-1)
	    # @show i, length(L[i])
	end
end

# ╔═╡ dfb312e7-4c8d-456b-aa91-8813781070bd
Lₙ[7]

# ╔═╡ 650903cb-3ec5-461c-8a84-6f431721d73f
begin
	xx=range(-1,stop=1,length=301)
	# Probajte razne vrijednosti k od 1 do 40
	kk=20
	yy=Lₙ[kk].(xx)
	Plots.plot(xx,yy)
end

# ╔═╡ 1ddfa099-4b1f-4be4-9954-8616d4d7d59a
md"""
## Čebiševljevi polinomi

__Čebiševljevi polinomi__ $T_n(x)$ nastaju ortogonalizacijom sustava (2) uz

$$
[a,b]=[-1,1], \quad \omega(x)=\frac{1}{\sqrt{1-x^2}}.$$

Čebiševljevi polinomi imaju sljedeća svojstva:

* vrijedi 

$$
T_n(x)=\cos(n\arccos x),\quad n=0,1,2,3,\ldots,$$

* $T_n(x)$  ima $n$ različitih nul-točaka na intervalu $[-1,1]$, 

$$
x_k=\cos \bigg(\frac{2k-1}{n}\frac{\pi}{2} \bigg), \quad k=1,\ldots,n,$$

* vrijedi __tročlana rekurzivna formula__: 

$$\begin{aligned}
T_0(x)&=1,\\
T_1(x)&=x, \\ 
T_{n+1}(x)&=2\,x\,T_n(x)-T_{n-1}(x),\quad n=1,2,3,\ldots.
\end{aligned}$$
 
__Napomena__:

Rekurzivna formula slijedi iz _adicione formule_

$$
\cos(n+1)\varphi+\cos(n-1)\varphi=2\cos\varphi\cos n\varphi.$$

Ortogonalnost se dokazuje pomoću supstitucije 

$$
\arccos x=\varphi.$$
"""

# ╔═╡ 9f7e3056-817e-4469-9807-db2c813b9f9e
begin
	# Simbolički
	T=Array{Any,1}(missing,n)
	T[1]=x^0
	T[2]=x
	for k=2:n-1
	    T[k+1]=expand(2*x*T[k]-T[k-1])
	end
end

# ╔═╡ 46efd454-fe33-4ae2-a354-125538158172
T[3]

# ╔═╡ e0768584-e9fd-46f4-9c45-53db7123db2f
T[7]

# ╔═╡ 2d3e07c8-e41c-4312-85af-cb08fdb820bd
T[8]

# ╔═╡ a9ca0f28-0a38-4aad-9685-0a92370ea433
begin
	# Numerički
	n₂=50
	Tₙ=Array{Any,1}(undef,n₂)
	Tₙ[1]=Polynomial([1])
	Tₙ[2]=Polynomial([0,1])
	for i=3:n₂
	    Tₙ[i]=2*Tₙ[2]*Tₙ[i-1]-Tₙ[i-2]
	    # @show i, length(T[i])
	end
end

# ╔═╡ c9a60766-b6a8-494a-886f-a3452d92a6d6
begin
	xx₁=range(-1,stop=1,length=300)
	# Probajte razne vrijednosti k od 1 do 50
	k=27
	yy₁=Tₙ[k].(xx₁)
	plot(xx₁,yy₁)
end

# ╔═╡ a2eda751-ca8f-48bc-ae2a-cfab8d24c2f0
md"""
## Promjena intervala

Ortogonalni sustav funkcija $\Phi_i$ na intervalu $[-1,1]$ pomoću transformacije 

$$
\gamma :[a,b]\to [-1,1],\quad \gamma(x)=\frac{2x}{b-a}-\frac{a+b}{b-a}$$

prelazi u ortogonalni sustav funkcija na intervalu $[a,b]$

$$
\Psi_i(x)=\Phi_i(\gamma(x)).$$
"""

# ╔═╡ 3665633e-4cea-423b-81bb-b7697c75ffd0
begin
	a₁=1
	b₁=4
	xx₂=collect(range(a,stop=b,length=300))
	γ=2*xx₂/(b₁-a₁).-(b₁+a₁)/(b₁-a₁)
	# Probajte razne vrijednosti k od 1 do 50
	k₂=17
	yy₂=Tₙ[17].(γ)
	# plot(xx₂,yy₂)
end

# ╔═╡ 57b29a00-1e7d-11eb-2846-03fb58c68173
Tₙ[8].(γ)

# ╔═╡ 20f09792-7e03-4003-8261-a34d3ed1b3a5
γ

# ╔═╡ Cell order:
# ╟─fdd0be01-c2a2-42c2-9a9b-1fcca24303f9
# ╟─a9eac96e-8ca5-464f-83d1-3ed9300e9d43
# ╠═934a0191-a567-497f-b8e5-61ca5ca63a30
# ╠═b999b693-8366-4d4a-8674-c080d250be3f
# ╟─bab6e737-737c-4c33-be1f-6ba60d45649f
# ╠═8cb49fec-a93a-4456-b159-fbe178b0cb67
# ╠═831cc7ee-31fb-42c5-be87-a75f971215fa
# ╠═3b90d2fd-afa0-466b-8fa5-e9748608ca00
# ╠═c240885a-008a-4ba3-8bae-b7dfaf1d1591
# ╠═d2dffde7-bee4-4b06-bd6d-a548b5dc2d5b
# ╟─6957e729-7be5-44e9-ae1c-3087bcb10a37
# ╠═0abd5a32-7398-4f1d-aa60-46e204713a75
# ╠═85a8addc-73bc-4e56-a8a9-3137ab443d15
# ╠═ea0adfdc-373f-40f7-a256-08d6424b6165
# ╠═e0b11a5c-becd-48b2-9238-f8dd8caed9ca
# ╠═c494410b-62dd-4822-8413-1e159b10d5c7
# ╠═2ac8adab-18a5-4b4e-bcb3-d4e777fd9c55
# ╠═11f09cbf-a386-43fb-9b02-1a0303ee4e0f
# ╟─2e40d6a7-a59c-4342-9246-1c347e715795
# ╠═83b09ece-0e8c-40c9-83a0-2204bba04323
# ╠═dfb312e7-4c8d-456b-aa91-8813781070bd
# ╠═650903cb-3ec5-461c-8a84-6f431721d73f
# ╟─1ddfa099-4b1f-4be4-9954-8616d4d7d59a
# ╠═9f7e3056-817e-4469-9807-db2c813b9f9e
# ╠═46efd454-fe33-4ae2-a354-125538158172
# ╠═e0768584-e9fd-46f4-9c45-53db7123db2f
# ╠═2d3e07c8-e41c-4312-85af-cb08fdb820bd
# ╠═a9ca0f28-0a38-4aad-9685-0a92370ea433
# ╠═c9a60766-b6a8-494a-886f-a3452d92a6d6
# ╟─a2eda751-ca8f-48bc-ae2a-cfab8d24c2f0
# ╠═57b29a00-1e7d-11eb-2846-03fb58c68173
# ╠═3665633e-4cea-423b-81bb-b7697c75ffd0
# ╠═20f09792-7e03-4003-8261-a34d3ed1b3a5
