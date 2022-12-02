### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 3a20c826-03cd-4c06-a0c6-7dc657067feb
using PlutoUI, Plots, ForwardDiff

# ╔═╡ 537a526f-067c-4a1a-b03c-5dab6edc68a6
TableOfContents(title="📚 Sadržaj", aside=true)

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
# Bisekcija

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
## Primjeri

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
# Jednostavne iteracije

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
## Primjer

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
## Korijen iz 2

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
# Newtonova metoda

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
# Metoda sekante

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ForwardDiff = "~0.10.24"
Plots = "~1.25.2"
PlutoUI = "~0.7.23"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "1cefa896a4a7b2c69b38af8b21fc8195080a9b27"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "8fb515c5a2c8941cef957e75afb99a2c24b753f3"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "65ebc27d8c00c84276f14aaf4ff63cbe12016c70"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═3a20c826-03cd-4c06-a0c6-7dc657067feb
# ╠═537a526f-067c-4a1a-b03c-5dab6edc68a6
# ╟─1198a836-bb14-4b8e-9f30-160097bc4507
# ╟─7942e1ac-57f4-493c-9d09-de4154dcd7b9
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
