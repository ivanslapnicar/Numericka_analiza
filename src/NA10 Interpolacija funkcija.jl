### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 6a0417bf-00aa-47db-bf22-5842f1cc736e
begin
	using Polynomials
	using Plots
end

# ╔═╡ df89a063-e647-4bfe-91eb-167be078ac0e
md"""
# Interpolacija funkcija


Neka je zadana funkcija $f(x)$ na intervalu $[a,b]$.

Odaberimo $n+1$ točku $x_i, i=0,\ldots, n$, u intervalu $[a,b]$ tako da je $x_i\neq x_j$
te kroz točke $T_i=(x_i,f(x_i))$ provucimo interpolacijski polinom.

__Teorem.__ Za svaku točku $x\in[a,b]$ vrijedi __ocjena pogreške__ (uz pretpostavku da funkcija $f$ ima $n+1$ derivaciju) 

$$\begin{aligned}
f(x)-p_n(x)&=\frac{\omega(x)}{(n+1)!} \,f^{(n+1)}(\xi), \cr
\omega(x)&=\prod_{k=0}^n (x-x_k)=(x-x_0)(x-x_1)\cdots (x-x_n),\quad  \xi \in (a,b).
\end{aligned}$$


_Dokaz._ (Vidi [Numerička matematika, str. 23](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).) 

Za $x=x_i$,  tvrdnja je očigledna. U suprotnom, definirajmo pomoćnu funkciju 

$$g(y)= f(y)-p_n(y)-k\omega(y),$$

pri čemu je konstanta $k$ odabrana tako da je $g(x)=0$.
Na taj način funkcija $g$ ima barem $n+2$ nultočke, 
$x,x_0,x_1,\ldots,x_n$. 
Prema Rolleovom teoremu, derivacija $g'$ ima barem $n+1$ nultočku, $g''$ ima barem $n$ nultočaka, itd. Funkcija $g^{(n+1)}$ ima barem jednu nultočku $\xi\in(a,b)$.

Vrijedi $p_n^{(n+1)}(y)=0$ i $\omega^{(n+1)}(y)=(n+1)!$ (vodeći koeficijent od $\omega$ je $1$). Uvrštavanje daje 

$$0=g^{(n+1)}(\xi)=f^{(n+1)}(\xi)-k(n+1)!$$ 

pa je $k=\displaystyle\frac{f^{(n+1)}(\xi)}{(n+1)!}$ i teorem je dokazan.


## Primjer

Promotrimo funkciju 

$$
f(x)=\sin(x), \quad x\in[0,\pi].$$
"""

# ╔═╡ c88bd3a2-a46d-4931-b4cf-0b940696aa89
begin
	# Ove funkcije omogućuju manipulaciju s Vandermondeovim matricama
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

# ╔═╡ 3a599b53-077b-47e7-a1cf-17e15da6aa1f
begin
	n=6
	a=0
	b=pi
	x=range(a,stop=b,length=n)
	f(x)=sin(x)
	y=f.(x)
end

# ╔═╡ de651f77-3953-4ada-ac74-48bd58d147c1
A=Vandermonde(x)

# ╔═╡ 7f63d94f-3792-4369-b8f4-edd3a5b0cd78
c=A\y

# ╔═╡ 7c1650eb-bf62-40af-b085-1b81da0f47bf
p=Polynomial(c)

# ╔═╡ af46e285-715a-4b5c-8745-d17d61806dbf
scatter(x,y)

# ╔═╡ 6a1c5393-f4c9-4a45-8a11-0b4c884406b7
begin
	x₀=range(a,stop=b,length=100)
	p₀=p.(x₀)
	F₀=f.(x₀)
	plot!(x₀,[p₀ F₀])
end

# ╔═╡ 9e80a29c-ea79-489f-8383-60fab341f5be
begin
	# maksimalne apsolutna i relativna pogreška
	using LinearAlgebra
	norm(p₀[2:end-1]-F₀[2:end-1],Inf), 
	norm((p₀[2:end-1]-F₀[2:end-1])./F₀[2:end-1],Inf)
end

# ╔═╡ 58277216-ffc5-47e7-a8c6-b181e54ec870
md"""
## Čebiševljeve točke

__Čebiševljevi polinomi__ su polinomi stupnja $n$ dani formulom

$$
T_n(x)=\cos(n\, \arccos x), \quad n=0,1,2,\ldots$$

Vrijedi rekurzivna formula:

$$
\begin{aligned}
T_0(x)&=1,\cr 
T_1(x)&=x,\cr 
T_{n+1}(x)&=2\,x\,T_n(x)-T_{n-1}(x), \quad n=1,2,\ldots
\end{aligned}$$

Dakle, 

$$
T_2(x)=2x^2-1,\quad T_3(x)=4x^3-3x, \ldots$$

Nul-točke polinoma $T_n(x)$ su 

$$
x_k=\cos \bigg(\frac{2k-1}{2n}\pi\bigg), \quad k=1,2,\ldots,n.$$

Sve nul-točke leže unutar intervala $[-1,1]$.

Na intervalu $[-1,1]$ polinom $T_n(x)$ poprima vrijednosti u intervalu $[-1,1]$.

__Napomena.__ Rekurzivna formula slijedi iz adicione formule

$$
\cos(n+1)\varphi+\cos(n-1)\varphi=2\cos\varphi \cos n\varphi$$

uz $\varphi=\arccos x$.

### Primjer
"""

# ╔═╡ 45533a9a-caa8-401b-9c22-40e9afb46bbe
T(n,x)=cos.(n*acos.(x))

# ╔═╡ 79a264ae-3635-466a-9a2c-4ffe612f417d
x₁=range(-1,stop=1,length=100)

# ╔═╡ 89136037-71ac-4dbe-9621-90f36c7b2702
y₁=T(10,x₁)

# ╔═╡ 62239acf-92c3-4b68-ae54-5e6a1dd4fb53
plot(x₁,y₁)

# ╔═╡ 4f65fc49-4d4e-43d7-9775-d53e1e766805
xₙ=[cos((2*k-1)*pi/(2*10)) for k=10:-1:1]

# ╔═╡ 16f339d9-01ad-4409-8fec-4ff797e8fd02
yₙ=T(10,xₙ)

# ╔═╡ f2a33009-575a-4912-9c68-223e6862407e
scatter!(xₙ,yₙ)

# ╔═╡ 2ee3a5fe-8da3-4a37-b61b-63fa0282d287
md"""
### Norme funkcija

Za funkcije 

$$f,g:[a,b]\to \mathbb{R}$$

definiramo __skalarni produkt__

$$
(f,g)=\int_a^b f(x)g(x)\, dx$$

i __težinski skalarni produkt__ s __težinom__ $\omega(x)>0$

$$
(f,g)_\omega=\int_a^b f(x)g(x)\omega(x)\, dx.$$

Funkcije $f$ i $g$ su __ortogonalne__ ako je $(f,g)=0$ ili ako je $(f,g)_\omega=0$.

Sljedeće tri __norme__ su prirodna poopćenja odgovarajućih vektorskih normi:

$$
\begin{aligned}
\|f\|_2&=\sqrt{(f,f)}=\sqrt{\int_a^b f^2(x)\, dx} \cr
\|f\|_1 &= \int_a^b \big|f(x)\big|\, dx \cr
\|f\|_\infty&=\max_{x\in[a,b]} \big|f(x)\big|
\end{aligned}$$

Vrijedi sljedeći važan teorem:

__Teorem__. Od svih polinoma stupnja manjeg ili jednakog $n$ čiji je koeficijent uz najveću potenciju jednak $1$, najmanju 
$\|\cdot\|_\infty$ na intervalu $[-1,1]$ ima upravo polinom 
$\displaystyle\frac{1}{2^{n-1}}T_n(x)$ i ta norma iznosi $\displaystyle\frac{1}{2^{n-1}}$.

_Dokaz_ : (Vidi [Numerička matematika, str. 101](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).)

Koeficijent Čebiševljevog polinoma $T_n(x)$ uz potenciju $x^n$ je $2^{n-1}$.
Stoga polinom $\displaystyle\frac{1}{2^{n-1}}T_n(x)$ uz potenciju $x^n$ ima koeficijent $1$. Također, zbog svojstva Čebiševljevih polinoma vrijedi

$$\left| \displaystyle\frac{1}{2^{n-1}}T_n(x) \right| \leq \displaystyle\frac{1}{2^{n-1}}$$

pa je

$$\left\| \displaystyle\frac{1}{2^{n-1}}T_n(x) \right\|_\infty = \displaystyle\frac{1}{2^{n-1}}.$$

Pretpostavimo da za polinom 

$$p_n(x)=x^n + \alpha_{n-1}x^{n-1} +\alpha_{n-2}x^{n-2}+\cdots \alpha_1 x+\alpha_0$$

vrijedi $| p_n(x)|< \displaystyle\frac{1}{2^{n-1}}$ za svaki $x\in[-1,1]$. Neka su $\xi_0,\xi_1,\ldots,\xi_n$ točke u kojima $T_n(x)$ poprima ekstremne vrijednosti $-1$ ili $1$. Vrijedi

$$
\begin{aligned}
p_n(\xi_0) & < \frac{1}{2^{n-1}}T_n(\xi_0) = \frac{1}{2^{n-1}} \cr
p_n(\xi_1) & > \frac{1}{2^{n-1}}T_n(\xi_1) = -\frac{1}{2^{n-1}} \cr
p_n(\xi_2) & < \frac{1}{2^{n-1}}T_n(\xi_2) = \frac{1}{2^{n-1}} \cr
& \vdots
\end{aligned}$$

To bi značilo da polinom $p_n(x)-\frac{1}{2^{n-1}}T_n(x)$ stupnja $n-1$ mijenja predznak $n$ puta, što je nemoguće, pa je teorem dokazan.


Zaključujemo da će polinomna aproksimacija (1) biti najbolja ako na intervalu $[-1,1]$ 
odaberemo

$$\omega(x)=\frac{1}{2^{n}} T_{n+1}(x),$$ 

odnosno ako na intervalu $[a,b]$ za 
točke interpolacije $x_0,x_1,\ldots,x_n$ odaberemo upravo nul-točke polinoma $T_{n+1}(x)$ preslikane na interval $[a,b]$.

### Promjena intervala

Sustav ortogonalnih funkcija $\Phi_i$ na intervalu $[-1,1]$ pomoću transformacije 

$$
\gamma :[a,b]\to [-1,1],\quad \gamma(x)=\frac{2x}{b-a}-\frac{a+b}{b-a}$$

prelazi u sustav ortogonalnih funkcija na intervalu $[a,b]$

$$
\Psi_i(x)=\Phi_i(\gamma(x)).$$

Nama je potrebna inverzna transformacija:

$$
x=\frac{a+b}{2}+\frac{b-a}{2}\gamma(x).$$
"""

# ╔═╡ 9bc340ed-75a2-4e50-b209-8075508d0611
# Odaberimo za interpolaciju sinusa nultočke polinoma T(n,x)
xₜ=(a+b)/2 .+(b-a)/2*[cos((2*k-1)*pi/(2*n)) for k=n:-1:1]

# ╔═╡ 79d7ee5b-4dcc-47c4-b858-f6b6afb8b66a
yₜ=f.(xₜ)

# ╔═╡ e6ea8fd8-6fbf-4fc7-ac6a-7354b1e0be85
begin
	Aₜ=Vandermonde(xₜ)
	cₜ=Aₜ\yₜ
	pₜ=Polynomial(cₜ)
end

# ╔═╡ b2c12a8c-374a-4222-8e37-83998743988a
begin
	# x₂=range(a,stop=b,length=100)
	p₂=pₜ.(x₀)
	F₂=f.(x₀)
	plot(x₀,[p₂ F₂])
end

# ╔═╡ e147f5fa-6826-4b21-80a4-1fbd14593c8c
# maksimalne apsolutna i relativna pogreška
norm(p₂[2:end-1]-F₂[2:end-1],Inf), 
norm((p₂[2:end-1]-F₂[2:end-1])./F₂[2:end-1],Inf)

# ╔═╡ 797d58e0-6aba-4171-82b2-574c08aae214
md"""
Pogledajmo kako izgledaju stvarne pogreške u oba slučaja:
"""

# ╔═╡ af090d7e-155d-11eb-1b3b-77d83ad0b048
# Ravnomjerno raspoređene točke
plot(x₀,p₀-F₀)

# ╔═╡ bba5f710-155d-11eb-1f00-892b62a32cb4
# Čebiševljeve točke
plot(x₀,p₂-F₂)

# ╔═╡ 7aa65294-2803-4bca-aa69-ab11ed6e43a3
md"""
Vidimo da su za Čebiševljeve točke postignute manje pogreške.

__Napomena.__ Ovdje smo, radi jednostavnosti, koristili najmanje točnu varijantu računanja interpolacijskog polinoma.

### Primjer 

Napravimo još jedan zanimljiv primjer (vidi [Numerička matematika, str. 24](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)):
interpolirajmo funkciju

$$
f(x)=1-|x-1|,\quad x\in[0,2]$$

polinomom stupnja 10.
"""

# ╔═╡ c6b0381f-dd60-42f3-8f6a-df8771bbf3c6
begin
	n₃=11
	a₃=0
	b₃=2
	f₃(x)=1 .-abs.(x .-1)
	# Ravnomjerno raspoređene točke
	x₃=range(a₃,stop=b₃,length=n₃)
	y₃=f₃(x₃)
	A₃=Vandermonde(x₃)
	c₃=A₃\y₃
	p₃=Polynomial(c₃)
	xp₃=range(a₃,stop=b₃,length=100)
	pf₃=p₃.(xp₃)
	F₃=f₃(xp₃)
	scatter(x₃,y₃)
	plot!(xp₃, [pf₃ F₃])
end

# ╔═╡ bcf8a882-eca1-4739-9789-43a9de96cb94
begin
	# Čebiševljeve točke
	xt₃=(a₃+b₃)/2 .+(b₃-a₃)/2*[cos((2*k-1)*pi/(2*n₃)) for k=n₃+1:-1:1]
	yt₃=f₃(xt₃)
	At₃=Vandermonde(xt₃)
	ct₃=At₃\yt₃
	pc₃=Polynomial(ct₃)
	pCheb=pc₃.(xp₃)
	scatter(xt₃,yt₃)
	plot!(xp₃,[pCheb F₃])
end

# ╔═╡ 3760da30-1560-11eb-083e-575c20a74101
xt₃

# ╔═╡ Cell order:
# ╟─df89a063-e647-4bfe-91eb-167be078ac0e
# ╠═6a0417bf-00aa-47db-bf22-5842f1cc736e
# ╠═c88bd3a2-a46d-4931-b4cf-0b940696aa89
# ╠═3a599b53-077b-47e7-a1cf-17e15da6aa1f
# ╠═de651f77-3953-4ada-ac74-48bd58d147c1
# ╠═7f63d94f-3792-4369-b8f4-edd3a5b0cd78
# ╠═7c1650eb-bf62-40af-b085-1b81da0f47bf
# ╠═af46e285-715a-4b5c-8745-d17d61806dbf
# ╠═6a1c5393-f4c9-4a45-8a11-0b4c884406b7
# ╠═9e80a29c-ea79-489f-8383-60fab341f5be
# ╟─58277216-ffc5-47e7-a8c6-b181e54ec870
# ╠═45533a9a-caa8-401b-9c22-40e9afb46bbe
# ╠═79a264ae-3635-466a-9a2c-4ffe612f417d
# ╠═89136037-71ac-4dbe-9621-90f36c7b2702
# ╠═62239acf-92c3-4b68-ae54-5e6a1dd4fb53
# ╠═4f65fc49-4d4e-43d7-9775-d53e1e766805
# ╠═16f339d9-01ad-4409-8fec-4ff797e8fd02
# ╠═f2a33009-575a-4912-9c68-223e6862407e
# ╟─2ee3a5fe-8da3-4a37-b61b-63fa0282d287
# ╠═9bc340ed-75a2-4e50-b209-8075508d0611
# ╠═79d7ee5b-4dcc-47c4-b858-f6b6afb8b66a
# ╠═e6ea8fd8-6fbf-4fc7-ac6a-7354b1e0be85
# ╠═b2c12a8c-374a-4222-8e37-83998743988a
# ╠═e147f5fa-6826-4b21-80a4-1fbd14593c8c
# ╟─797d58e0-6aba-4171-82b2-574c08aae214
# ╠═af090d7e-155d-11eb-1b3b-77d83ad0b048
# ╠═bba5f710-155d-11eb-1f00-892b62a32cb4
# ╟─7aa65294-2803-4bca-aa69-ab11ed6e43a3
# ╠═c6b0381f-dd60-42f3-8f6a-df8771bbf3c6
# ╠═bcf8a882-eca1-4739-9789-43a9de96cb94
# ╠═3760da30-1560-11eb-083e-575c20a74101
