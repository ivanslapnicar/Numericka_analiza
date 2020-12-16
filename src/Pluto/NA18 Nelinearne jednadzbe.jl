### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 3a20c826-03cd-4c06-a0c6-7dc657067feb
using Plots

# ╔═╡ a2c79625-bb3d-4cb3-972c-14fa1763d1e6
using ForwardDiff

# ╔═╡ 1198a836-bb14-4b8e-9f30-160097bc4507
md"""
# Nelinearne jednadžbe


__Problem.__ Nađimo nul-točke funkcije $f(x)$ na zatvorenom intervalu $[a,b]$, odnosno, riješimo jednadžbu

$$
f(x)=0, \quad x\in[a,b]. \qquad\qquad (1)$$

Vrijedi sljedeće:

Ako je $f$ __neprekidna__ na $[a,b]$ i ako je $f(a)\cdot f(b)<0$, tada postoji barem jedna točka $\xi\in(a,b)$ takva da je 

$$
f(\xi)=0.$$

Ako je još i $f'(x)\neq 0$ za $x\in(a,b)$, tada je $\xi$ __jedinstvena__.

Stoga jednadžbu (1) možemo riješiti u dva koraka:

1. Nađemo interval $[a,b]$ u kojem funkcija $f$ ima jedinstvenu nultočku $\xi$,
2. Aproksimiramo točku $\xi$ s unaprijed zadanom točnošću.

Opisat ćemo četiri metode:

1. metodu bisekcije,
2. metodu jednostavih iteracija,
3. Newtonovu metodu (metodu tangente) i 
4. metodu sekante.

Sve metode, uz zadanu početnu aproksimaciju $x_0$,  generiraju niz točaka $x_n$ koji, uz određene uvjete, konvergira
prema rješenju $\xi$. 

Metoda ima __red konvergencije__ jednak $r>0$ ako postoji $A>0$ takav da je

$$
|\xi-x_{n+1}|\leq A|\xi-x_n|^r.$$

__Napomena.__ Dokazi tvrdnji i primjeri se nalaze u knjizi [Numerička matematika, poglavlje 4.1](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).

"""

# ╔═╡ 7942e1ac-57f4-493c-9d09-de4154dcd7b9
md"""
## Bisekcija

Počevši od intervala $[a,b]\equiv [a_0,b_0]$, konstruiramo niz intervala 

$$
[a_0,b_0]\supset [a_1,b_1]\supset [a_2,b_2]\supset [a_3,b_3] \supset \cdots,$$

gdje je $f(a_n)\cdot f(b_n)\leq 0$, i niz točaka

$$
x_{n+1}=\frac{a_n+b_n}{2}.$$

__Brzina konvergencije__ je __linearna__ jer je 

$$
|\xi-x_{n+1}|\leq \frac{1}{2}|\xi-x_n|,$$

a __pogreška aproksimacije__ je omeđena s

$$
|\xi-x_{n+1}|\leq \frac{1}{2}|a_n-b_n|.$$
"""

# ╔═╡ 9d707705-29a2-47fc-87b8-d39ff967efb5
function Bisection(f::Function,a::Number,b::Number,ϵ::Float64=1e-10)
    fa=f(a)
    fb=f(b)
    T=Float64
    x=T
    fx=T
    if fa*fb>zero(T)
        return "Netočan interval"
    end
    iter=0
    while b-a>ϵ && iter<1000
        x=(b+a)/2.0
        fx=f(x)
        if fa*fx<zero(T)
            b=x
            fb=fx
        else
            a=x
            fa=fx
        end
        iter+=1
        # @show x,fx
    end
    return x,fx,iter
end

# ╔═╡ 0233a80f-dfc6-420c-ae39-60828b19010a
md"""
### Primjeri

Zadane su funkcije i intervali:

$$
\begin{aligned}
f_1(x)&=e^x-x- \frac{5}{4},\quad &x\in [-2,2],\\
f_2(x)&=e^{-2x}\sin (6x)-\frac{2}{3}\,x-\frac{1}{2},\quad &x\in[-1,2],\\
f_3(x)&=x^3-6x+2,\quad &x\in[-4,4],\\
f_4(x)&=0.001\,x+0.5+\frac{\pi}{2}+\arctan(x),\quad &x\in[-1000,1000],\\
f_5(x)&=1000\,(x-4)-e^x,\quad &x\in[-10,10].
\end{aligned}$$
"""

# ╔═╡ 0be4ad52-497f-48a4-9f53-b6d1291beb58
begin
	f₁(x)=exp(x)-x-5.0/4
	(a₁,b₁)=(-1,1)
	f₂(x)=exp(-2x)*sin(6x)+2x/3-1.0/2
	(a₂,b₂)=(-1,2)
	f₃(x)=x^3-6*x+2
	(a₃,b₃)=(-4,4)
	f₄(x)=0.001x+0.5+π/2+atan(x)
	(a₄,b₄)=(-1000,1000)
	f₅(x)=1000(x-4)-exp(x)
	(a₅,b₅)=(-10,10)
end

# ╔═╡ 5a17ece9-a499-42c7-8dc8-cc9035f1b311
md"""
Pomoću grafa funkcije odredimo intervale u kojima se nalaze nul-točke, koje potom izračunamo i nacrtamo.
"""

# ╔═╡ 910d7790-488d-48b8-a951-eee68774f79f
function NulTočke(f,a,b,Intervali)
    plot(f,a,b,label="f(x)")
    for i=1:length(Intervali)
        iab=Intervali[i]
        x,y,iter=Bisection(f,iab[1],iab[2])
        scatter!([x],[y],label="Nul-točka")
    end
    scatter!()
end

# ╔═╡ 5d233a49-774b-479c-8c80-5de46817614e
# Funkcija f₁(x)
plot(f₁,a₁,b₁,label="f(x)")

# ╔═╡ aedcd7b9-8ef5-48a8-8ef8-d7dab4b733e0
Intervali₁=((-1,0),(0,1))

# ╔═╡ e7e92425-287a-482b-bb1e-8168b2b130f1
NulTočke(f₁,a₁,b₁,Intervali₁)

# ╔═╡ b00698a0-3a02-11eb-11fe-c70ec5102e23
# Računanje jedne nul-točke
x,fx,iter=Bisection(f₁,-1,0,1e-15)

# ╔═╡ 56a6bc1f-7ff8-464d-a864-ec7a19e3f9c3
# Funkcija f₂(x)
plot(f₂,a₂,b₂,label="f(x)")

# ╔═╡ 7097ec9a-b556-4be7-a356-9c9d5003dd73
Intervali₂=((-1,-0.4),(-0.4,0.2),(0.2,0.6),(0.6,1))

# ╔═╡ bdb9d5d9-861c-4fd1-89d5-e0622120f944
NulTočke(f₂,a₂,b₂,Intervali₂)

# ╔═╡ f7802590-744f-4175-ba60-5d801fa252ec
# Računanje jedne nultočke
Bisection(f₂,0.2,0.6)

# ╔═╡ d193fd54-c5a7-45f4-a61d-ae532bae4f44
# Funkcija f₃(x)
plot(f₃,a₃,b₃,label="f(x)")

# ╔═╡ 52fbcdb5-ade7-48c0-bc09-0154f3d615b8
begin
	Intervali₃=((-4,-2),(0,1),(2,3))
	NulTočke(f₃,a₃,b₃,Intervali₃)
end

# ╔═╡ e1f25e30-1ff7-4a91-b7cf-393c36c438b6
begin
	# Funkcija f₄(x)
	Intervali₄=[(-600,-400)]
	NulTočke(f₄,a₄,b₄,Intervali₄)
end

# ╔═╡ f33b37c9-2f88-422b-ac9f-40d7945ccb47
begin
	# Funkcija f₅(x)
	Intervali₅=((0,5),(5,10))
	NulTočke(f₅,a₅,b₅,Intervali₅)
end

# ╔═╡ 4c7dc70c-1832-461c-ad48-5e358b49661f
plot!(legend=:bottomright)

# ╔═╡ 69c1d80b-3d24-4ded-a226-391b337805e5
md"""
## Jednostavne iteracije

Rješavamo jednadžbu oblika 

$$
x=\varphi(x).  \qquad\qquad (2)$$

__Teorem o fiksnoj točki. (Banach)__
Neka je 

$$\varphi:[a,b]\to \mathbb{R}$$

__neprekidno derivabilna funkcija__ i neka vrijedi

$$\begin{aligned}
\varphi(x) &\in [a,b], \quad  \forall x\in [a,b], \\
|\varphi'(x)|&\leq q<1, \quad \forall x\in(a,b).
\end{aligned}\qquad\qquad (3)$$

Tada postoji jedinstvena __fiksna točka__ $\xi \in [a,b]$ za koju vrijedi
$\xi=\varphi(\xi)$. 

Nadalje, za proizvoljnu početnu točku  $x_0\in[a,b]$ niz 

$$
x_n=\varphi(x_{n-1}),\quad n=1,2,3,\ldots,$$

konvergira prema $\xi$ te vrijede __ocjene pogreške:__

$$\begin{aligned}
|\xi-x_n|&\leq \displaystyle\frac{q^n}{1-q}|x_1-x_0|, \\
|\xi-x_n|&\leq \displaystyle\frac{q}{1-q}|x_n-x_{n-1}|, \\
|\xi-x_n|&\leq q|\xi-x_{n-1}|.
\end{aligned}$$

Dakle, konvergencija je __linearna__.

Za dokaz teorema vidi [R. Scitovski, Numerička matematika, str. 73](https://www.mathos.unios.hr/nm/materijali/Num.PDF).
"""

# ╔═╡ d794399b-9c03-47a1-8a37-40f8545a0d46
function Iteration(φ::Function,x::Number,ϵ::Float64=1e-10)
    ξ=φ(x)
    iter=0
    while abs(x-ξ)>ϵ && iter<1000
        x=ξ
        ξ=φ(x)
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 5cc494f2-a789-4455-9ac9-925428d274fe
md"""
Za korištenje metode iteracije potrebno je
transformirati oblik (1) u oblik (2) i to tako da je ispunjen uvjet (3).

Za procjenu derivacije možemo koristiti paket `Calculus.jl` koji aproksimira derivaciju konačnim razlikama ili paket
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) koji koristi [automatsku diferencijaciju](https://en.wikipedia.org/wiki/Automatic_differentiation) i koji je točniji. Može se koristiti i simboličko računanje pomoću paketa `SymPy.jl`.
"""

# ╔═╡ 21495df3-1163-4d7f-afce-ad3600a5411a
varinfo(ForwardDiff.ForwardDiff)

# ╔═╡ 492a1939-ea1e-4912-b8d3-0330624f0ffd
md"""
### Primjer

Potražimo nul-točke funkcije $f_1(x)=e^x-x-\frac{5}{4}$. Iz oblika

$$
x=e^x-\frac{5}{4}\equiv \varphi(x)$$

možemo izračunati samo negativnu nul-točku, jer je u okolini pozitivne nul-točke $|\varphi'(x)|>1$.
Za $x_0=1$, niz divergira vrlo brzo, a za $x_0=0.6$, što je blizu pozitivne nul-točke, 
niz konvergira prema negativnoj nul-točki, i to bez teoretskog obrazloženja.
"""

# ╔═╡ 7dab9934-414d-44dd-94b5-96a849184687
begin
	φ(x)=exp(x)-5.0/4
	plot([f₁,φ,x->ForwardDiff.derivative(φ,x)],-2.0,2.0,label=["f(x)" "φ(x)" "φ'(x)"])
end

# ╔═╡ 1e5be7af-8236-4a8c-ab10-4f988350bd6a
Iteration(φ,-0.5)

# ╔═╡ ed701b40-1ca7-4eba-8ea4-7c2beef3e58d
Iteration(φ,1.0)

# ╔═╡ 8eb3ad86-6d79-4016-9f65-5797f51d779e
Iteration(φ,0.6)

# ╔═╡ f1079ba0-3f85-11eb-3762-c99121d890f9
md"

Pozitivnu nul-točku možemo izračunati iz oblika 

$$x=\ln \left(x+\frac{5}{4}\right)\equiv\psi(x).$$
"

# ╔═╡ 70f9b3bc-ef11-4d7c-af8c-0abc262b92fa
begin
	ψ(x)=log(x+5.0/4)
	plot([f₁,x->ForwardDiff.derivative(φ,x), x->ForwardDiff.derivative(ψ,x)],
	    -1.0,1.0,label=["f(x)" "φ'(x)" "ψ'(x)"])
end

# ╔═╡ 3ccbe741-478e-4b23-8a64-d6f3380d8b1b
scatter!([Iteration(φ,-0.5)[1],Iteration(ψ,1.0)[1]],[0,0],label="Nul-točke")

# ╔═╡ 5dee83f3-2171-478a-b811-04c3ccef1a0e
md"""
### Korijen iz 2

Izračunajmo približno $\sqrt{2}$, odnosno izračunajmo pozitivno rješenje jednadžbe 

$$
x^2-2=0.$$

Jednadžbu je moguće pretvoriti u oblik (2) kao 

$$
x=\frac{2}{x},$$

no tada je $\varphi'(x)=-\displaystyle\frac{2}{x^2}$ pa na intervalu $[1,2]$ ne vrijedi (3). Zato stavimo

$$
\frac{x}{2}=\frac{1}{x},$$

odnosno

$$
x=\frac{x}{2}+\frac{1}{x}=\frac{1}{2}\left(x+\frac{2}{x}\right)\equiv\varphi(x).$$

Točna vrijednost se postiže nakon samo 4 iteracije!
"""

# ╔═╡ b27404ef-ac18-4dce-804c-07255ab001c2
begin
	φ₁(x)=(x+2.0/x)/2.0
	plot([φ₁,x->ForwardDiff.derivative(φ₁,x)],1.0,2.0,label=["φ₁(x)" "φ₁'(x)"])
end

# ╔═╡ 25640f91-4b47-4bcf-b360-499541077c80
Iteration(φ₁,1.0,1e-15), sqrt(2)

# ╔═╡ 18d4d7c0-3a0d-11eb-087d-5db48ab8a323
# Ručni račun s racionalnim brojevima
begin
	y=1//1
	y1=(y+2//y)//2
	y2=(y1+2//y1)//2
	y3=(y2+2//y2)//2
	y4=(y3+2//y3)//2
	y5=(y4+2//y4)//2
end

# ╔═╡ 4a7aa6ae-3a0d-11eb-24b0-3f207580daa2
886731088897/627013566048-sqrt(2)

# ╔═╡ cf6ff8ea-2c1e-442a-bfc6-5eac4f7269dd
begin
	# Probajmo i sqrt(10)
	φ₂(x)=(9x+10.0/x)/10.0
	plot([φ₂,x->ForwardDiff.derivative(φ₂,x)],3.0,4.0,label=["φ₂(x)" "φ₂'(x)"])
end

# ╔═╡ 4fa06267-97f2-409f-8025-c160362014a9
Iteration(φ₂,3.0,1e-10),sqrt(10) # 1e-15

# ╔═╡ 46bebd7f-69ba-402f-91aa-6ad864add60b
begin
	# Probajmo sqrt(10) na drugi način
	φ₃(x)=(4x+10.0/x)/5.0
	plot([φ₃,x->ForwardDiff.derivative(φ₃,x)],3.0,4.0)
end

# ╔═╡ 9b237357-f6dc-4054-bee8-1cac4a5aad07
Iteration(φ₃,3.0,1e-10), sqrt(10) # 1e-15

# ╔═╡ 629f59b0-3f86-11eb-1e5d-c181f17a6209
md"
__Zadatak.__ Izvedite gornje formule za računanje $\sqrt{10}$. Kako glasi općenita formula za računanje $\sqrt{n}$, $n\in\mathbb{N}$? 
"

# ╔═╡ 4f3b3c3a-0cf4-4bcd-a8ad-9061c7653386
md"""
## Newtonova metoda

__Newtonova metoda__ ili __metoda tangente__ temelji se na sljedećoj ideji: zadanu funkciju $f(x)$ u okolini zadane početne aproksimacije $x_0$ aproksimiramo tangentom kroz točku $(x_0,f(x_0))$,

$$
f_1(x)=f(x_0)+f'(x_0)(x-x_0),$$

te za sljedeću aproksimaciju uzmemo sjecište tangente s $x$-osi. Na taj dobijemo niz aproksimacija:

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)},\quad n=0,1,2,\ldots \qquad\qquad (4)$$

Vrijedi sljedeći

__Teorem.__  Neka je zadana funkcija $f:[a,b]\to \mathbb{R}$ za koju vrijedi:

*  $f''$ je neprekidna na $(a,b)$,
*  $f(a)\cdot f(b)<0$,
*  $f'$ i $f''$ imaju stalan predznak na $(a,b)$, i 
*  $f(x_0)\cdot f''(x_0)>0$ za odabranu početnu aproksimaciju $x_0\in [a,b]$.

Tada niz (4) konvergira prema __jedinstvenom__ rješenju $\xi$ jednadžbe $f(x)=0$. 
Pri tome vrijede __ocjene pogreške__:

$$
\begin{aligned}
|\xi-x_n|&\leq \displaystyle\frac{M_2}{2m_1}(x_n-x_{n-1})^2, \\
|\xi-x_{n+1}|&\leq \displaystyle\frac{M_2}{2m_1}(\xi-x_{n})^2, \\
\end{aligned}$$

gdje je 

$$
M_2=\max_{x\in(a,b)}|f''(x)|,\quad
m_1=\min_{x\in(a,b)}|f'(x)|.$$

Dakle, konvergencija je __kvadratična__.
"""

# ╔═╡ d8d1e5e8-13c6-43cd-b2ca-aa6e96c142e3
function Newton(f::Function,x::Number,ϵ::Float64=1e-10)
    ξ=x-f(x)/(x->ForwardDiff.derivative(f,x))(x)
    iter=0
    while abs(x-ξ)>ϵ && iter<100
        x=ξ
        ξ=x-f(x)/(x->ForwardDiff.derivative(f,x))(x)
		println(ξ)
        iter+=1
    end
    ξ,iter
end

# ╔═╡ b61f1530-3f86-11eb-3fd2-0be4beedaa3c
md"

Izračunajmo nul-točke funkcije

$$f_6(x)=e^{-x}+x^2-2.$$
"

# ╔═╡ 6258bcb3-2682-47e5-8db5-0aaf03f32807
begin
	f₆(x)=exp(-x)+x^2-2
	plot(f₆,-3,4,label="f(x)")
end

# ╔═╡ f9376e6e-0595-4045-be73-c59c07e7be34
md"""
Provjerimo uvjete teorema za pozitivnu nul-točku:
"""

# ╔═╡ dca6feb2-e2e3-4c1e-a568-794c610cc611
begin
	a₀=1
	b₀=2
	x₀=1.5
	plot([f₆,x->ForwardDiff.derivative(f₆,x),
	        x->ForwardDiff.derivative(x->ForwardDiff.derivative(f₆,x),x)],a₀,b₀, 
	    label=["f(x)" "f'(x)" "f''(x)"])
end

# ╔═╡ 0cf32332-f169-465a-9315-b9704a31badf
f₆(a₀)*f₆(b₀)<0, 
f₆(x₀)*(x->ForwardDiff.derivative(
        x->ForwardDiff.derivative(f₆,x),x))(x₀)>0

# ╔═╡ 6de6ad42-cad9-40b4-8503-78f1316b3b38
Newton(f₆,x₀) # 1e-15

# ╔═╡ 9ba461d0-3b8c-4ec7-9b14-8ea720a569c3
begin
	# Negativna nul-točka
	a=-1
	b=0
	x₁=-1.0
	f₆(a)*f₆(b)<0, 
	f₆(x₁)*(x->ForwardDiff.derivative(
	        x->ForwardDiff.derivative(f₆,x),x))(x₁)>0
end

# ╔═╡ 9609f70c-e5e9-4552-8ecb-8b95cf484d29
Newton(f₆,x₁)

# ╔═╡ 518e1f14-3254-4ca6-9d87-c695008ae79e
begin
	plot(f₆,-1.0,2)
	scatter!([Newton(f₆,x₀)[1],Newton(f₆,x₁)[1]],[0,0])
end

# ╔═╡ b7109008-05de-45db-a383-78b8a6973733
md"""
__Napomena.__ Ukoliko za početne aproksimacije odaberemo vrijednosti $x_0=1$, odnosno $x_0=0$, metoda će također konvergirati prema željenim nul-točkama, premda bez teoretskog obrazloženja: 
"""

# ╔═╡ 0fbd88ae-6c19-41ac-983b-83ba8bd97965
Newton(f₆,1)

# ╔═╡ 917b4fc3-e485-4801-8d16-0e32a22649b5
Newton(f₆,0)

# ╔═╡ 1db0ba3d-e82e-4e4d-8b5d-c1100d04fd2d
md"""
## Metoda sekante

Ukoliko u formuli (4) derivaciju $f'(x_n)$ aproksimiramo konačnom razlikom (sekantom) kroz __dvije__ prethodne točke,

$$
f'(x_n)\approx \frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}, $$

dobit ćemo niz

$$
x_{n+1}=\frac{x_{n-1}f(x_n)-x_nf(x_{n-1})}{f(x_n)-f(x_{n-1})},\qquad f(x_n)\neq f(x_{n-1}), \quad n=1,2,3,\ldots.$$

Na početku trebamo odabrati __dvije__ početne aproksimacije, $x_0,x_1\in[a,b]$. Svojstva konvergencije su slična onima Newtonove metode.
"""

# ╔═╡ dce97964-5846-4642-859e-f041b8a72d76
function Sekanta(f::Function,x::Number,ζ::Number,ϵ::Float64=1e-10)
    ξ=(x*f(ζ)-ζ*f(x))/(f(ζ)-f(x))
    iter=0
    while abs(ζ-ξ)>ϵ && iter<100
        x=ζ
        ζ=ξ
        ξ=(x*f(ζ)-ζ*f(x))/(f(ζ)-f(x))
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 33402e76-139d-4062-b72d-57b45e1c3eec
Sekanta(f₆,-1,0), Sekanta(f₆,1,2)

# ╔═╡ 2d8b7b30-3f88-11eb-3add-f36e29cb5511
md"
__Zadatak.__ Nađite nul-točke funkcija $f_1,\ldots,f_5$ Newtonovom metodom i metodom sekante.
"

# ╔═╡ Cell order:
# ╟─1198a836-bb14-4b8e-9f30-160097bc4507
# ╟─7942e1ac-57f4-493c-9d09-de4154dcd7b9
# ╠═3a20c826-03cd-4c06-a0c6-7dc657067feb
# ╠═9d707705-29a2-47fc-87b8-d39ff967efb5
# ╟─0233a80f-dfc6-420c-ae39-60828b19010a
# ╠═0be4ad52-497f-48a4-9f53-b6d1291beb58
# ╟─5a17ece9-a499-42c7-8dc8-cc9035f1b311
# ╠═910d7790-488d-48b8-a951-eee68774f79f
# ╠═5d233a49-774b-479c-8c80-5de46817614e
# ╠═aedcd7b9-8ef5-48a8-8ef8-d7dab4b733e0
# ╠═e7e92425-287a-482b-bb1e-8168b2b130f1
# ╠═b00698a0-3a02-11eb-11fe-c70ec5102e23
# ╠═56a6bc1f-7ff8-464d-a864-ec7a19e3f9c3
# ╠═7097ec9a-b556-4be7-a356-9c9d5003dd73
# ╠═bdb9d5d9-861c-4fd1-89d5-e0622120f944
# ╠═f7802590-744f-4175-ba60-5d801fa252ec
# ╠═d193fd54-c5a7-45f4-a61d-ae532bae4f44
# ╠═52fbcdb5-ade7-48c0-bc09-0154f3d615b8
# ╠═e1f25e30-1ff7-4a91-b7cf-393c36c438b6
# ╠═f33b37c9-2f88-422b-ac9f-40d7945ccb47
# ╠═4c7dc70c-1832-461c-ad48-5e358b49661f
# ╟─69c1d80b-3d24-4ded-a226-391b337805e5
# ╠═d794399b-9c03-47a1-8a37-40f8545a0d46
# ╟─5cc494f2-a789-4455-9ac9-925428d274fe
# ╠═a2c79625-bb3d-4cb3-972c-14fa1763d1e6
# ╠═21495df3-1163-4d7f-afce-ad3600a5411a
# ╟─492a1939-ea1e-4912-b8d3-0330624f0ffd
# ╠═7dab9934-414d-44dd-94b5-96a849184687
# ╠═1e5be7af-8236-4a8c-ab10-4f988350bd6a
# ╠═ed701b40-1ca7-4eba-8ea4-7c2beef3e58d
# ╠═8eb3ad86-6d79-4016-9f65-5797f51d779e
# ╟─f1079ba0-3f85-11eb-3762-c99121d890f9
# ╠═70f9b3bc-ef11-4d7c-af8c-0abc262b92fa
# ╠═3ccbe741-478e-4b23-8a64-d6f3380d8b1b
# ╟─5dee83f3-2171-478a-b811-04c3ccef1a0e
# ╠═b27404ef-ac18-4dce-804c-07255ab001c2
# ╠═25640f91-4b47-4bcf-b360-499541077c80
# ╠═18d4d7c0-3a0d-11eb-087d-5db48ab8a323
# ╠═4a7aa6ae-3a0d-11eb-24b0-3f207580daa2
# ╠═cf6ff8ea-2c1e-442a-bfc6-5eac4f7269dd
# ╠═4fa06267-97f2-409f-8025-c160362014a9
# ╠═46bebd7f-69ba-402f-91aa-6ad864add60b
# ╠═9b237357-f6dc-4054-bee8-1cac4a5aad07
# ╟─629f59b0-3f86-11eb-1e5d-c181f17a6209
# ╟─4f3b3c3a-0cf4-4bcd-a8ad-9061c7653386
# ╠═d8d1e5e8-13c6-43cd-b2ca-aa6e96c142e3
# ╟─b61f1530-3f86-11eb-3fd2-0be4beedaa3c
# ╠═6258bcb3-2682-47e5-8db5-0aaf03f32807
# ╟─f9376e6e-0595-4045-be73-c59c07e7be34
# ╠═dca6feb2-e2e3-4c1e-a568-794c610cc611
# ╠═0cf32332-f169-465a-9315-b9704a31badf
# ╠═6de6ad42-cad9-40b4-8503-78f1316b3b38
# ╠═9ba461d0-3b8c-4ec7-9b14-8ea720a569c3
# ╠═9609f70c-e5e9-4552-8ecb-8b95cf484d29
# ╠═518e1f14-3254-4ca6-9d87-c695008ae79e
# ╟─b7109008-05de-45db-a383-78b8a6973733
# ╠═0fbd88ae-6c19-41ac-983b-83ba8bd97965
# ╠═917b4fc3-e485-4801-8d16-0e32a22649b5
# ╟─1db0ba3d-e82e-4e4d-8b5d-c1100d04fd2d
# ╠═dce97964-5846-4642-859e-f041b8a72d76
# ╠═33402e76-139d-4062-b72d-57b45e1c3eec
# ╟─2d8b7b30-3f88-11eb-3add-f36e29cb5511
