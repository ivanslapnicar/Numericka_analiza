### A Pluto.jl notebook ###
# v0.12.8

using Markdown
using InteractiveUtils

# ╔═╡ cda9b377-09e7-4a7a-8f84-f50de595a7c0
using Plots

# ╔═╡ ec0b9159-1c36-4112-8ce7-ec69c6709861
md"""
# Regresija

__Regresija__ je provlačenje funkcije $f$ koja ovisi o $n$ parametara kroz točke 

$$(x_i,y_i),\quad i=1,2,\ldots, m,$$

pri čemu je $m>n$, tako da se __minimizira norma odstupanja__:

$$
\| f(x_i)-y_i\|_{1,2,\infty}\to \min.$$

Regresija u __smislu najmanjih kvadrata__ je

$$
\| f(x_i)-y_i\|_{2}\to \min.$$

Kada je funkcija $f$ pravac,

$$
f(x)=kx+l,$$

govorimo o __linearnoj regresiji__. U tom slučaju dobije se sustav linearnih jednadžbi

$$
k x_i + l=y_i, \quad i=1,2,\ldots,m.$$

Ukoliko sve točke __ne leže na istom pravcu__, sustav nije rješiv pa se računa kvadratična prilagodba.
"""

# ╔═╡ 108def07-d919-498c-82c2-7c43c7fd27f3
md"""
## Linearna regresija

Provucimo pravac kroz točke 

$$(1,1), \ (2,3),\ (4,2), \ (6,4), \ (7,3),$$

i izračunajmo kvalitetu prilagodbe. 
"""

# ╔═╡ 64ebb9c6-9c26-42bb-9141-45f5efab3d36
begin
	x=[1,2,4,6,7]
	y=[1,3,2,4,3]
	A=[x ones(Int,length(x))]
end

# ╔═╡ 5a57b1b6-eabb-46be-9986-1c48dcaf1785
begin
	# Koeficijenti regrecijskog pravca
	using LinearAlgebra
	(k,l)=A\y
end

# ╔═╡ 748be150-4eee-4c29-93fe-bdc1058cb3e8
begin
	# Nacrtajmo točke i regresijski pravac
	f(x)=k*x+l
	scatter(x,y,label="Točke",legend=:bottomright)
	plot!(x->x,f,x[1],x[end],label="Regresijski pravac")
end

# ╔═╡ 1b894260-f180-4407-8c68-ad93ad40c55b
# Kvaliteta prilagodbe
sqrt(norm(A*[k;l]-y)/norm(y))

# ╔═╡ 986e6cb9-4a13-49e2-9891-29020262d002
md"""
## Kvadratična regresija

Kroz točke možemo provući i kvadratni polinom $y=ax^2+bx+c$. Ukoliko sve točke ne leže na istoj paraboli, sustav linearnih jednadžbi 

$$
ax_i^2+bx_i+c=y_i, \quad i=1,\ldots,m,$$

nije rješiv pa računamo kvadratičnu prilagodbu.

Provucimo kvadratni polinom kroz točke 

$$
(1,0),\ (2,1), \ (4,4),\ (5,8), \ (6,14).$$
"""

# ╔═╡ a786bb3e-e442-4deb-a7ce-ed25fbf6bd75
begin
	x₁=[1,2,4,5,6]
	y₁=[0,1,4,8,14]
	A₁=[x₁.^2 x₁ ones(Int,length(x₁))]
end

# ╔═╡ a98c4043-4f99-4aa4-977a-2f5f777e762b
# Koeficijenti regresijskog polinoma
(a,b,c)=A₁\y₁

# ╔═╡ 544441c1-f921-4767-859b-804f5fbd9508
# Nacrtajmo točke i parabolu
g(x)=a*x^2+b*x+c

# ╔═╡ 93fdb6ef-09d2-47b7-b5af-6561529f520f
begin
	scatter(x₁,y₁,label="Točke",legend=:topleft)
	plot!(x->x,g,x₁[1],x₁[end],label="Regresijski polinom")
end

# ╔═╡ 65f7a8eb-df16-4a01-91ff-0214ac092e61
# Kvaliteta prilagodbe
sqrt(norm(A₁*[a;b;c]-y₁)/norm(y₁))

# ╔═╡ b6166388-1d3b-4dce-969a-b7fce38426c1
md"""
## Rast svjetske populacije

Dosadašnji rast populacije (u milionima) da je u sljedećoj tablici (vidi [http://en.wikipedia.org/wiki/World_population](http://en.wikipedia.org/wiki/World_population)):

$$
\begin{array}{c|c|c|c|c|c|c|c|c|c}
\textrm{godina} & 1750 & 1800 & 1850 & 1900 & 1950 & 1999 & 2008 & 2010 & 2012 \\ \hline
\textrm{populacija (milijuni)} & 791 & 978 & 1262 & 1650 & 2521 & 5978 & 6707 & 6896 & 7052 
\end{array}$$




Aproksimirajmo rast populacije eksponencijalnom funkcijom 

$$
P(t)=Ce^{kt}$$

i predvidimo populaciju 2050. godine.

Sustav jednadžbi 

$$
Ce^{kt_i}=P_i, \quad i=1,2,\ldots, 9,$$

logaritmiranjem prelazi u sustav linearnih jednadžbi

$$
k \,t_i + \ln C =\ln P_i.$$

Sve točke ne leže na istoj krivulji pa sustav nije rješiv i računamo kvadratičnu prilagodbu.
"""

# ╔═╡ d86cd66b-5e87-4f06-b991-a9c734e5b9fb
begin
	nₚ=9
	t=[1750,1800,1850,1900,1950,1999,2008,2010,2012]
	P=[791,978,1262,1650,2521,5978,6707,6896,7052]
	Aₚ=[t ones(Int,length(t))]
	(kₚ,C)=Aₚ\log.(P)
end

# ╔═╡ f1b3366d-c0c2-499c-86b7-c34f86d7998a
# Vrijednosti na krivulji
P₁(t)=exp(C)*exp(kₚ*t)

# ╔═╡ 998e1154-e985-49a9-8cce-4b436a887e2f
begin
	# Nacrtajmo točke i regresijsku krivulju
	scatter(t,P,label="Populacija",legend=:topleft)
	plot!(t->t,P₁,t[1],t[end],label="Regresijska krivulja")
end

# ╔═╡ f48de770-1eb0-11eb-0d18-8f9690cd89c6
md"
__Pitanje.__ Zbog čega stvarna krivulja populacije ima lom?
"

# ╔═╡ 1ac0491b-6cc9-4a14-af99-c6c6defc3074
# Predvidimo populaciju 2050 godine
P₁(2050)

# ╔═╡ 23e34ea5-8dc6-4444-bcd6-94229abad290
md"
Izračunato predviđanje je manje od onog u tablici. Ako se ograničimo na razdoblje od 1950 godine imamo:
"

# ╔═╡ 52e78105-19f7-4830-bf59-4e336daa9272
begin
	A₂= [t[5:end] ones(5)]
	(k₂,C₂)=A₂\log.(P[5:end])
	P₂(t)=exp(C₂)*exp(k₂*t)
	scatter(t[5:end],P[5:end],label="Populacija",legend=:topleft)
	plot!(t->t,P₂,t[5],t[end],label="Regresijska krivulja")
end

# ╔═╡ e8911abc-74bd-4daa-9336-5b10270616b1
P₂(2050)

# ╔═╡ Cell order:
# ╟─ec0b9159-1c36-4112-8ce7-ec69c6709861
# ╟─108def07-d919-498c-82c2-7c43c7fd27f3
# ╠═cda9b377-09e7-4a7a-8f84-f50de595a7c0
# ╠═64ebb9c6-9c26-42bb-9141-45f5efab3d36
# ╠═5a57b1b6-eabb-46be-9986-1c48dcaf1785
# ╠═748be150-4eee-4c29-93fe-bdc1058cb3e8
# ╠═1b894260-f180-4407-8c68-ad93ad40c55b
# ╟─986e6cb9-4a13-49e2-9891-29020262d002
# ╠═a786bb3e-e442-4deb-a7ce-ed25fbf6bd75
# ╠═a98c4043-4f99-4aa4-977a-2f5f777e762b
# ╠═544441c1-f921-4767-859b-804f5fbd9508
# ╠═93fdb6ef-09d2-47b7-b5af-6561529f520f
# ╠═65f7a8eb-df16-4a01-91ff-0214ac092e61
# ╟─b6166388-1d3b-4dce-969a-b7fce38426c1
# ╠═d86cd66b-5e87-4f06-b991-a9c734e5b9fb
# ╠═f1b3366d-c0c2-499c-86b7-c34f86d7998a
# ╠═998e1154-e985-49a9-8cce-4b436a887e2f
# ╟─f48de770-1eb0-11eb-0d18-8f9690cd89c6
# ╠═1ac0491b-6cc9-4a14-af99-c6c6defc3074
# ╟─23e34ea5-8dc6-4444-bcd6-94229abad290
# ╠═52e78105-19f7-4830-bf59-4e336daa9272
# ╠═e8911abc-74bd-4daa-9336-5b10270616b1
