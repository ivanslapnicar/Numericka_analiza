### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 695660b9-ebbd-4eb1-8a4d-c7ee6eb8d31a
using PlutoUI, ForwardDiff, Plots, LinearAlgebra, NLsolve, Optim

# ╔═╡ 8b2f7699-fc87-44fc-9ac7-0677ed7e7d73
plotly()

# ╔═╡ 25ba3a7c-067f-413c-b04e-90e457631aba
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 94da0d5c-6dfe-487a-9711-9d101f20e97a
md"""
# Sustavi nelinearnih jednadžbi


__Problem.__ Nađimo rješenje $\xi=(\xi_1,\xi_2,\ldots,\xi_n)$ sustava od $n$ jednadžbi 

$$\begin{aligned}
f_1(x)&=0,\cr
f_2(x)&=0,\cr
&\vdots \cr
f_n(x)&=0,
\end{aligned}$$

i $n$ nepoznanica $x=(x_1,x_2,\ldots,x_n)$. Uz oznaku $f=(f_1,f_2,\ldots,f_n)^T$, ovaj sustav možemo zapisati kao 

$$
f(x)=0.$$

Opisat ćemo __Newtonovu metodu__ i tri __kvazi-Newtonove__ metode:

1. __Broydenovu__ metodu,
2. __Davidon-Fletcher-Powell__ metodu i 
3. __Broyden-Fletcher-Goldfarb-Schano__ metodu.

Sve metode, uz zadanu početnu aproksimaciju $x^{(0)}$,  generiraju niz točaka $x^{(n)}$ koji, uz određene uvjete, konvergira prema rješenju $\xi$. 

__Napomena.__ Opisi metoda i primjeri se nalaze u knjizi [Numerička matematika, poglavlje 4.4](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).
"""

# ╔═╡ 125c9230-cc68-493f-92c5-8a83a78b5863
md"""
# Newtonova metoda

__Jacobijan__ ili __Jacobijeva matrica__ funkcija $f$ u točki $x$ je matrica prvih parcijalnih derivacija

$$
J(f,x)=\begin{pmatrix} \displaystyle\frac{\partial f_1(x)}{\partial x_1} & \displaystyle\frac{\partial f_1(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_1(x)}{\partial x_n} \\
\displaystyle\frac{\partial f_2(x)}{\partial x_1} & \displaystyle\frac{\partial f_2(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_2(x)}{\partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\displaystyle\frac{\partial f_n(x)}{\partial x_1} & \displaystyle\frac{\partial f_n(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_n(x)}{\partial x_n} 
\end{pmatrix}.$$

Za zadanu početnu aproksimaciju $x^{(0)}$, računamo niz točaka

$$
x^{(k+1)}=x^{(k)}-s^{(k)}, \quad k=0,1,2,\ldots,$$

gdje je $s^{(k)}$ rješenje sustava

$$
J\big(f,x^{(k)}\big)\cdot s=f\big(x^{(k)}\big).$$

Za računanje Jacobijana koristimo paket [`ForwardDiff.jl`](http://www.juliadiff.org/ForwardDiff.jl/perf_diff.html#derivatives). Za crtanje funkcija koristimo paket `Plots.jl`.
"""

# ╔═╡ 77095a05-eb64-47e2-a2e0-5053a4fbc837
function Newton(f::Function,J::Function,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,size(x))
    ξ=x
    while norm(s)>ϵ && iter<100
        s=J(x)\f(x)
        ξ=x-s
        iter+=1
        x=ξ
    end
    ξ,iter
end

# ╔═╡ ed8971fd-1420-41af-b736-e052439a2c65
md"""
## Primjer 1

(Dennis i Schnabel, 1996) Rješenja sustava

$$\begin{aligned}
2(x+y)^2+(x-y)^2-8&=0\\
5x^2+(y-3)^2-9&=0
\end{aligned}$$

su točke $T_1=(1,1)$ i $T_2\approx(-1.18,1.59)$.
"""

# ╔═╡ 96eb13a8-ca15-447b-a252-76d481b30912
# Vektorska funkcija
f₁(x)=[2(x[1]+x[2])^2+(x[1]-x[2])^2-8,5*x[1]^2+(x[2]-3)^2-9]

# ╔═╡ 7366924f-29e6-46a6-a5f1-e1940daa15dd
f₁((1.0,2))

# ╔═╡ d7930d33-99f3-4551-8258-500bbde17aff
md"""
Nacrtajmo funkcije i konture kako bi mogli približno locirati nul-točke:
"""

# ╔═╡ 46060676-7192-4e7d-9038-9492524c1485
begin
	# Broj točaka
	m=101
	X=range(-2,stop=3,length=m)
	Y=range(-2,stop=3,length=m)
	# Prva aplikata
	surface(X,Y,(x,y)->f₁([x,y])[1],xlabel="x",ylabel="y",colorbar=false)
	# Druga aplikata
	surface!(X,Y,(x,y)->f₁([x,y])[2],seriescolor=:blues)
end

# ╔═╡ 0234d953-b83d-48a1-9c32-a0e2b8bec20a
begin
	# Odredimo rješenja pomoću kontura
	contour(X,Y,(x,y)->f₁([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->f₁([x,y])[2],contour_labels=true)
end

# ╔═╡ 7d4c9549-1233-4e03-ace4-1a1c87147036
# Jasnija slika
contour!(clims=(0,0.01),xlabel="x",ylabel="y",colorbar=:none)

# ╔═╡ e2b8b699-39e0-408a-8165-2f606ea25740
md"""
Vidimo da su nul-točke približno $T_1=(-1.2,1.5)$ i $T_2=(1,1)$. Štoviše, $T_2$ je točno jednaka $(1,1)$ (1 iteracija u trećem primjeru). Nadalje, metoda ne mora konvergirati (četvrti primjer).   
"""

# ╔═╡ 08505d4c-b799-4833-af0d-4efea3f925c1
J₁(x)=ForwardDiff.jacobian(f₁,x)

# ╔═╡ a5631cd4-c6b7-4a7b-9411-adf83f7c15b3
# Na primjer
J₁([1.0,2])

# ╔═╡ 2a26cb77-ab97-4723-b404-4eea958111e9
Newton(f₁,J₁,[-1.0,0.0])

# ╔═╡ ce3dbb73-bee6-4ffe-adb5-4bc25c067435
Newton(f₁,J₁,[0.5,1.1])

# ╔═╡ 8ee7b19c-5eaa-4188-9ea7-fe2c8dc38e4a
Newton(f₁,J₁,[1.0,1.0])

# ╔═╡ 961033d5-2a2a-4417-9523-a09164b42c2c
Newton(f₁,J₁,[0.0,0.0])

# ╔═╡ 3cfa6d61-75f0-4897-97d1-b94c8de2458e
md"""
## Primjer 2

(Dennis i Schnabel, 1996) Rješenja sustava

$$\begin{aligned}
x_1^2-x_2^2-2&=0\\
e^{x_1-1}+x_2^3-2&=0
\end{aligned}$$

su točke $T_1=(1,1)$ i $T_2\approx (-0.71,1.22)$ .
"""

# ╔═╡ 2a1b5d9b-672f-4e00-8620-daec0196d6b7
begin
	f₂(x)=[x[1]^2+x[2]^2-2,exp(x[1]-1)+x[2]^3-2]
	contour(X,Y,(x,y)->f₂([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->f₂([x,y])[2],contour_labels=true)
	contour!(clims=(0,0.01),xlabel="x",ylabel="y",colorbar=:none)
end

# ╔═╡ 6c8b54c2-d4ae-4ada-b0d9-cd83eafa9cc1
begin
	J₂(x)=ForwardDiff.jacobian(f₂,x)
	Newton(f₂,J₂,[-1.0,1]), Newton(f₂,J₂,[0.8,1.2])
end

# ╔═╡ 2f2a6ae5-755d-4f65-9c16-4db2e1aa48f4
md"""
## Primjer 3

(Dennis i Schnabel, 1996) Zadan je problem $f(x)=0$, gdje je

$$
f(x)=\begin{pmatrix}x_1 \\ x_2^2+x_2 \\ e^{x_3}-1 \end{pmatrix}.$$

Točna rješenja su $T_1=(0,0,0)$ i $T_2=(0,-1,0)$. Izračunat ćemo nul-točke s nekoliko početnih aproksimacija.
"""

# ╔═╡ c6f98f59-312d-498d-9854-4c78c969011f
begin
	f₃(x)=[x[1],x[2]^2+x[2],exp(x[3])-1]
	J₃(x)=ForwardDiff.jacobian(f₃,x)
end

# ╔═╡ ee2d1d89-d399-4d1f-8d08-c99e372d5537
Newton(f₃,J₃,[-1.0,1.0,0.0]),Newton(f₃,J₃,[1.0,1,1]),
Newton(f₃,J₃,[-1.0,1,-10]),Newton(f₃,J₃,[0.5,-1.5,0])

# ╔═╡ e195e6ad-f3ea-4a74-a91a-d16d0384a054
md"""
## Primjer 4

(Rosenbrock parabolic valley) Zadana je funkcija

$$
f(x)=100\,(x_2-x_1)^2+(1-x_1)^2.$$

Tražimo moguće ekstreme funkcije, odnosno želimo riješiti jednadžbu

$$
\mathop{\mathrm{grad}} f(x)=0.$$
"""

# ╔═╡ 6e1b5ed8-3d23-4aad-b140-1567c8dc7336
f₄(x)=100(x[2]-x[1]^2)^2+(1-x[1])^2

# ╔═╡ d9704e21-9e89-4b5d-bb4e-7dd24d6b506a
# Nacrtajmo funkciju koristeći X i Y iz Primjera 1
surface(X,Y,(x,y)->f₄([x,y]), seriescolor=:blues, xlabel="x", ylabel="y")

# ╔═╡ 008c4d22-4c2f-4087-8be7-47599600de47
begin
	# Funkcija je zahtjevna u smislu određivanje ekstrema
	g₄(x)=ForwardDiff.gradient(f₄,collect(x))
	contour(X,Y,(x,y)->g₄([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->g₄([x,y])[2],contour_labels=true)
	contour!(clims=(-0.5,0.5),xlabel="x",ylabel="y",colorbar=:none)
end

# ╔═╡ ac8f76d2-9801-4f64-8add-278c8ada3135
md"""
Iz kontura vidimo da je primjer numerički zahtjevan, dok analitički lako vidimo da je jedina nul-točka $x_1=(1,1)$.

U ovom primjeru je vektorska funkcija $f$ zadana kao gradijent skalarne funkcije pa Jacobijevu matricu računamo korištenjem funkcije `FowardDiff.hessian()` koja računa aproksimaciju matrice drugih parcijalnih derivacija polazne funkcije $f$. 
"""

# ╔═╡ 34576988-89ce-431b-9e1c-78b1a98e002e
Newton(g₄,x->ForwardDiff.hessian(f₄,x),[-1.0,2.0])

# ╔═╡ 9625fb8d-5e0e-4e07-bc49-b4cb9ab259ef
md"""
## Primjer 5

Zadana je fukcija

$$
f(x)=\sum_{i=1}^{11} \bigg(x_3 \cdot \exp\bigg(-\frac{(t_i-x_1)^2}{x_2}\bigg)-y_i\bigg)^2,
$$

gdje su brojevi $(t_i,y_i)$ zadani tablicom:

$$\begin{array}{c|c|c|c|c|c|c|c|c|c|c|c}
 i & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 \\ \hline
t_i & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 &10 \\
 y_i & 0.001 & .01 & .04 & .12 & .21 & .25 & .21 & .12 & .04 & .01 & .001 
\end{array}$$

Želimo riješiti jednadžbu 

$$
\mathop{\mathrm{grad}} f(x)=0.
$$

U prethodnim zadatcima, kondicija Jacobijena ja

$$\kappa(J)=O(10)$$ 

u zadacima (a), (b) i (c) i 

$$\kappa(J)=O(1000)$$ 

u zadatku (d). U ovom zadatku je hesseove matrice  

$$\kappa(J)>O(10^3).$$ 

Unatoč tome, metoda je netočna i ne konvergira prema točnom rješenju $x=(4.93,2.62,0.28)$.
"""

# ╔═╡ 33ca0d51-92fe-4313-ba67-7792af17b11e
begin
	t=collect(0:10)
	y=[0.001,0.01,0.04,0.12,0.21,0.25,0.21,0.12,0.04,0.01,0.001]
	f₅(x)=sum([(x[3]*exp(-(t[i]-x[1])^2/x[2])-y[i])^2 for i=1:11])
end

# ╔═╡ e5d2cee7-c2d4-45a5-901c-206be7d2ff58
begin
	# Početna točka je vrlo blizu rješenja
	x₀=[4.9,2.6,0.2]
	g₅(x)=ForwardDiff.gradient(f₅,x)
	J₅(x)=ForwardDiff.hessian(f₅,x)
	f₅(x₀), g₅(x₀), cond(J₅(x₀))
end

# ╔═╡ ab87421c-e498-44ff-af23-13ea14a0c5d3
x₅,iter₅=Newton(g₅,J₅,x₀,1e-10)

# ╔═╡ d9df727c-a96f-4af4-ab5a-60715662adce
g₅(x₅)

# ╔═╡ 82fffd1c-9ec0-459c-b033-d926641268dc
Newton(g₅,J₅,[4.9,2.62,0.28],1e-12)

# ╔═╡ e823348a-1ea3-45b7-b1ad-d51096f224dc
md"""
# Broydenova metoda

Za odabranu početnu aproksimaciju $x_0$ i matricu $B_0$, za $k=0,1,2,\ldots$, računamo redom:

$$\begin{aligned}
B_k \cdot s_k & = -f(x_k) \quad \textrm{(sustav)}\\
x_{k+1}&=x_{k}+s_k\\
y_k&=f(x_{k+1})-f(x_{k})\\
B_{k+1}&=B_k+\frac{(y_k-B_ks_k)s_k^T}{s_k\cdot s_k}
\end{aligned}$$

Na ovaj način izbjegavamo računanje Jacobijeve matrice u svakom koraku.
Možemo uzeti $B_0=J(x_0)$, ali i neku drugu matricu.
"""

# ╔═╡ 73ad34ac-e4fc-4e35-9d96-25969b37a4d6
function Broyden(f::Function,B::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<100
        s=-(B\f(x))
        ξ=x+s
        y=f(ξ)-f(x)
        B=B+(y-B*s)*(s/(s⋅s))'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 6059b96d-8ba0-475c-9929-8cc613355562
# Primjer 1
Broyden(f₁,J₁([-1.0,0.0]),[-1.0,0.0]), 
Broyden(f₁,J₁([1.0,1.5]),[1.0,1.5])

# ╔═╡ 64d397ac-ed6f-426c-9139-af21222723da
begin
	# Objasnite ponašanje metode kada za početnu matricu uzmemo jediničnu matricu!
	eye(n)=Matrix{Float64}(I,n,n)
	Broyden(f₁,eye(2),[-1.0,0.0]), Broyden(f₁,eye(2),[1.0,1.5]),
	Broyden(f₁,eye(2),[-1,1.5])
end

# ╔═╡ ba09dd4f-b114-45bb-9b6a-41c2ed9dbec5
begin
	# Primjer 2
	x0=[-1.0,1]
	x1=[0.8,1.2]
	Broyden(f₂,J₂([-1.0,1]),[-1.0,1]), Broyden(f₂,J₂([0.8,1.2]),[0.8,1.2])
end

# ╔═╡ 52506f6c-069a-4883-8160-1cf415e80c36
# Primjer 3
Broyden(f₃,J₃([-1.0,1,0]),[-1.0,1,0]), Broyden(f₃,J₃([0.5,-1.5,0]),[0.5,-1.5,0])

# ╔═╡ c0fb2bfc-7221-428d-b13d-ebb483bd0c1c
# Primjer 4
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([-1.0,2]),[-1.0,2]), # ali
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([1,2.0]),[-1.0,2]),
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([0.8,0.5]),[0.8,0.5])

# ╔═╡ 1fe960aa-f12a-4d3e-994a-1a5d5c10f29e
begin
	# Primjer 5
	x₆,iter₆=Broyden(g₅,(x->ForwardDiff.hessian(f₅,x))([4.9,2.6,0.2]), [4.9,2.6,0.2])
end

# ╔═╡ b4e246f8-72b9-4042-8160-a1eda2084e45
g₅(x₆)

# ╔═╡ 11259db5-f754-47c7-9e92-600adf6896f8
md"""
# Davidon-Fletcher-Powell (DFP) metoda

DFP je optimizacijska metoda koja traži točke ekstrema funkcije 
$F:\mathbb{R}^n \to \mathbb{R}$, u kojem slučaju je $f(x)=\mathop{\mathrm{grad}}F(x)$.

Za odabranu početnu aproksimaciju $x_0$ i matricu $H_0$, za $k=0,1,2,\ldots$, računamo redom:

$$\begin{aligned}
s_k&=-H_k f(x_k)\\
\beta_k&=\mathop{\mathrm{arg\ min}}_\beta F(x_{k}+\beta s_k) \\
s_k&=\beta_k s_k\\
x_{k+1}&=x_{k}+s_k \\
y_k&=f(x_{k+1})-f(x_{k})\\
H_{k+1}&=H_k+ \frac{s_k s_k^T}{y_k\cdot s_k}-\frac{H_k y_k y_k^T H_k}{y_k\cdot (H_k y_k)}.
\end{aligned}$$

Za matricu $H_0$ možemo uzeti jediničnu matricu, a za izvršavanje iteracije nije potrebno rješavati sustav linearnih jednadžbi, već se sva ažuriranja vrše s $O(n^2)$ operacija.

Jednodimenzionalnu minimizaciju po pravcu $x_{k}+\beta s_k$ računamo tako što metodom bisekcije nađemo nul-točke usmjerene derivacije.
"""

# ╔═╡ ff33e971-4871-4533-95ef-4fd9b14f32d1
function Bisekcija(f::Function,a::T,b::T,ϵ::Float64=1e-10) where T
    fa=f(a)
    fb=f(b)
    x=T
    fx=T
    if fa*fb>zero(T)
        # return "Netočan interval"
        if abs(fa)>abs(fb)
            return b,fb,0
        else
            return a,fa,0
        end
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
    x,fx,iter
end

# ╔═╡ f29daf79-022a-4f7d-8616-6cac0e7ea9e8
function DFP(f::Function,H::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<50
        s=-H*f(x)
        s0=s/norm(s)
        F(ζ)=f(x+ζ*s)⋅s0
        β,fx,iterb=Bisekcija(F,0.0,1.0,10*eps())
        s*=β
        ξ=x+s
        y=f(ξ)-f(x)
        z=H*y
        H=H+(s/(y⋅s))*s'-(z/(y⋅z))*z'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 7f7672da-787f-4628-83aa-e9297f7da322
md"""
__Primjer.__ Nađimo točku ekstrema funkcije 

$$
f(x,y)=(x+2y-7)^2+(2x+y-5)^2.$$

Funkcija ima minimum u točki $(1,3)$.
"""

# ╔═╡ 8f4fc906-99a7-4e06-8620-83cf427da8b9
f₆(x) = (x[1] + 2*x[2]-7)^2 + (2*x[1] + x[2]-5)^2

# ╔═╡ 7a61c7b0-5582-11eb-2359-05273b8b829f
ForwardDiff.gradient(f₆,[1,2])

# ╔═╡ b687e934-69d9-4608-b73f-0f3676c210ec
f₆([1,2])

# ╔═╡ 04b24200-5583-11eb-04bd-250c19ce9959
ForwardDiff.hessian(f₆,[1,2])

# ╔═╡ 7ebf3fb9-85c7-4ea0-a18a-a4617375082d
g₆(x)=ForwardDiff.gradient(f₆,x)

# ╔═╡ 1f6187a2-41a2-43e0-9cd3-322a069222aa
DFP(g₆,eye(2),[0.8,2.7],eps())

# ╔═╡ d1e5bfe4-1ac9-4dbe-aafb-f81822a72d3a
# Primjer 4
DFP(g₄,eye(2),[0.9,1.1])

# ╔═╡ 02ef6f4f-3259-447e-878a-b116a3d0eb4a
# Primjer 5
DFP(g₅,eye(3),[4.9,2.6,0.2])

# ╔═╡ f7483843-9cc9-4799-ad56-1e74bd0a08e7
md"""
# Broyden-Fletcher-Goldfarb-Schano (BFGS) metoda

BFGS je optimizacijska metoda koja uspješno traži točke ekstrema funkcije 
$F:\mathbb{R}^n \to \mathbb{R}$, u kojem slučaju je $f(x)=\mathop{\mathrm{grad}}F(x)$.

Metoda je slična DFP metodi, s nešto boljim svojstvima konvergencije.

Neka je zadana funkcija $F:\mathbb{R}^n \to \mathbb{R}$, čiji minimum tražimo, i neka je 
$f(x)=\mathop{\mathrm{grad}} F(x)$.

Za odabranu početnu aproksimaciju $x_0$ i matricu $H_0$, za $k=0,1,2,\ldots$, računamo redom:

$$\begin{aligned}
s_k&=-H_k f(x_k)\\
\beta_k&=\mathop{\mathrm{arg\ min}}F(x_{k}+\beta_k s_k) \\
s_k&=\beta_k s_k\\
x_{k+1}&=x_{k}+s_k \\
y_k&=f(x_{k+1})-f(x_{k})\\
H_{k+1}&=\bigg(I-\frac{s_k y_k^T}{y_k\cdot s_k}\bigg)H_k
\bigg( I-\frac{y_k s_k^T}{y_k\cdot s_k}\bigg)+\frac{s_k s_k^T}{y_k\cdot s_k}.
\end{aligned}$$

Za matricu $H_0$ možemo uzeti jediničnu matricu, a za izvršavanje iteracije nije potrebno rješavati sustav linearnih jednadžbi, već se sva ažuriranja vrše s $O(n^2)$ operacija.

Jednodimenzionalnu minimizaciju po pravcu $x_{k}+\beta s_k$ računamo tako što metodom bisekcije tražimo nul-točke usmjerene derivacije.
"""

# ╔═╡ b985704f-8f07-43b6-baa9-6b3d4541de27
function BFGS(f::Function,H::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<50
        s=-H*f(x)
        s0=s/norm(s)
        F(ζ)=f(x+ζ*s)⋅s0
        β,fx,iterb=Bisekcija(F,0.0,1.0,10*eps())
        s*=β
        ξ=x+s
        y=f(ξ)-f(x)
        z=H*y
        α=y⋅s
        s1=s/α
        H=H-s1*z'-z*s1'+s1*(y⋅z)*s1'+s1*s'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 31a1b2aa-7b5b-46b0-82ec-f75111b72a0d
BFGS(g₆,eye(2),[0.8,2.7],eps())

# ╔═╡ a0e4f6af-d4ab-4b04-9863-49de95eec947
# Primjer 4
BFGS(g₄,eye(2),[0.9,1.1])

# ╔═╡ c20593d7-885f-4f36-984c-ecded0884ecc
# Primjer 5
BFGS(g₅,eye(3),[4.9,2.6,0.2])

# ╔═╡ cef37fa4-df23-4267-bde4-c1fe6ddade23
md"""
# Julia paketi

Prethodni programi su jednostavne ilustrativne implementacije navedenih algoritama.
Julia ima paket [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) za rješavanje sustava nelineranih jednadžbi i paket [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) za nelinearnu optimizaciju.
"""

# ╔═╡ 92d27828-c97b-41d0-85e9-0d3b1738dfed
# Primjer 1
function f₁!(fvec,x)
    fvec[1] = 2(x[1]+x[2])^2+(x[1]-x[2])^2-8
    fvec[2] = 5*x[1]^2+(x[2]-3)^2-9
end

# ╔═╡ 6055ba96-dd58-420b-854e-277b60a3c1e7
nlsolve(f₁!,[-1.0,0])

# ╔═╡ ca176d00-70f9-42b4-a3f4-fd5d37c9f31a
nlsolve(f₁!,[0.5,1.1])

# ╔═╡ a8080f5e-6797-482e-84ca-7ddfcd40264d
# Primjer 2
function f₂!(fvec,x)
	fvec[1] = x[1]^2+x[2]^2-2
	fvec[2] = exp(x[1]-1)+x[2]^3-2
end

# ╔═╡ d4764a31-38b7-4002-92fb-bfd3580cf634
nlsolve(f₂!,[-1.0,1]), nlsolve(f₂!,[0.8,1.2])

# ╔═╡ 2a644f10-2683-4ea8-bd25-95343abb7e48
# Primjer 3
function f₃!(fvec,x)
	fvec[1] = x[1]
	fvec[2] = x[2]^2+x[2]
	fvec[3] = exp(x[3])-1
end

# ╔═╡ ca57bf59-f348-439a-8622-ea0edcc86f17
	nlsolve(f₃!,[-1.0,1.0,0.0]), nlsolve(f₃!,[1.0,1,1]),
	nlsolve(f₃!,[-1.0,1,-10]), nlsolve(f₃!,[0.5,-1.5,0])

# ╔═╡ bc477197-2ac5-4798-a62b-62f89791fb94
# Primjer 4
opt=optimize(f₄,[-1.0,2],Optim.BFGS())

# ╔═╡ 8ff25700-5584-11eb-1f95-abf9735d479e
opt.minimizer[1]

# ╔═╡ 90043461-3f6a-45e6-a661-57c8df3a2a93
# Primjer 5 - opet ne konvergira prema rješenju ???
opt₅=optimize(f₅,[4.0,2.5,0.27],Optim.BFGS())

# ╔═╡ a107fc1f-f0db-42cb-ad4a-ffd51726808a
opt₅.minimizer

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ForwardDiff = "~0.10.24"
NLsolve = "~4.5.1"
Optim = "~1.5.0"
Plots = "~1.25.3"
PlutoUI = "~0.7.27"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "8130e26e874d411eabac7777d37b410a37af62d6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "37b730f25b5662ac452f7bb2c50a0567cbb748d4"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
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
version = "1.0.1+0"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

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

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

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
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

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

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

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

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

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

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "35d435b512fbab1d1a29138b5229279925eba369"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.5.0"

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

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

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
git-tree-sha1 = "e4fe0b50af3130ddd25e793b471cb43d5279e3e6"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "7eda8e2a61e35b7f553172ef3d9eaa5e4e76d92e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "fed057115644d04fba7f4d768faeeeff6ad11a60"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.27"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

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
# ╠═695660b9-ebbd-4eb1-8a4d-c7ee6eb8d31a
# ╠═8b2f7699-fc87-44fc-9ac7-0677ed7e7d73
# ╠═25ba3a7c-067f-413c-b04e-90e457631aba
# ╟─94da0d5c-6dfe-487a-9711-9d101f20e97a
# ╟─125c9230-cc68-493f-92c5-8a83a78b5863
# ╠═77095a05-eb64-47e2-a2e0-5053a4fbc837
# ╟─ed8971fd-1420-41af-b736-e052439a2c65
# ╠═96eb13a8-ca15-447b-a252-76d481b30912
# ╠═7366924f-29e6-46a6-a5f1-e1940daa15dd
# ╟─d7930d33-99f3-4551-8258-500bbde17aff
# ╠═46060676-7192-4e7d-9038-9492524c1485
# ╠═0234d953-b83d-48a1-9c32-a0e2b8bec20a
# ╠═7d4c9549-1233-4e03-ace4-1a1c87147036
# ╟─e2b8b699-39e0-408a-8165-2f606ea25740
# ╠═08505d4c-b799-4833-af0d-4efea3f925c1
# ╠═a5631cd4-c6b7-4a7b-9411-adf83f7c15b3
# ╠═2a26cb77-ab97-4723-b404-4eea958111e9
# ╠═ce3dbb73-bee6-4ffe-adb5-4bc25c067435
# ╠═8ee7b19c-5eaa-4188-9ea7-fe2c8dc38e4a
# ╠═961033d5-2a2a-4417-9523-a09164b42c2c
# ╟─3cfa6d61-75f0-4897-97d1-b94c8de2458e
# ╠═2a1b5d9b-672f-4e00-8620-daec0196d6b7
# ╠═6c8b54c2-d4ae-4ada-b0d9-cd83eafa9cc1
# ╟─2f2a6ae5-755d-4f65-9c16-4db2e1aa48f4
# ╠═c6f98f59-312d-498d-9854-4c78c969011f
# ╠═ee2d1d89-d399-4d1f-8d08-c99e372d5537
# ╟─e195e6ad-f3ea-4a74-a91a-d16d0384a054
# ╠═6e1b5ed8-3d23-4aad-b140-1567c8dc7336
# ╠═d9704e21-9e89-4b5d-bb4e-7dd24d6b506a
# ╠═008c4d22-4c2f-4087-8be7-47599600de47
# ╟─ac8f76d2-9801-4f64-8add-278c8ada3135
# ╠═34576988-89ce-431b-9e1c-78b1a98e002e
# ╟─9625fb8d-5e0e-4e07-bc49-b4cb9ab259ef
# ╠═33ca0d51-92fe-4313-ba67-7792af17b11e
# ╠═e5d2cee7-c2d4-45a5-901c-206be7d2ff58
# ╠═ab87421c-e498-44ff-af23-13ea14a0c5d3
# ╠═d9df727c-a96f-4af4-ab5a-60715662adce
# ╠═82fffd1c-9ec0-459c-b033-d926641268dc
# ╟─e823348a-1ea3-45b7-b1ad-d51096f224dc
# ╠═73ad34ac-e4fc-4e35-9d96-25969b37a4d6
# ╠═6059b96d-8ba0-475c-9929-8cc613355562
# ╠═64d397ac-ed6f-426c-9139-af21222723da
# ╠═ba09dd4f-b114-45bb-9b6a-41c2ed9dbec5
# ╠═52506f6c-069a-4883-8160-1cf415e80c36
# ╠═c0fb2bfc-7221-428d-b13d-ebb483bd0c1c
# ╠═1fe960aa-f12a-4d3e-994a-1a5d5c10f29e
# ╠═b4e246f8-72b9-4042-8160-a1eda2084e45
# ╟─11259db5-f754-47c7-9e92-600adf6896f8
# ╠═ff33e971-4871-4533-95ef-4fd9b14f32d1
# ╠═f29daf79-022a-4f7d-8616-6cac0e7ea9e8
# ╟─7f7672da-787f-4628-83aa-e9297f7da322
# ╠═8f4fc906-99a7-4e06-8620-83cf427da8b9
# ╠═7a61c7b0-5582-11eb-2359-05273b8b829f
# ╠═b687e934-69d9-4608-b73f-0f3676c210ec
# ╠═04b24200-5583-11eb-04bd-250c19ce9959
# ╠═7ebf3fb9-85c7-4ea0-a18a-a4617375082d
# ╠═1f6187a2-41a2-43e0-9cd3-322a069222aa
# ╠═d1e5bfe4-1ac9-4dbe-aafb-f81822a72d3a
# ╠═02ef6f4f-3259-447e-878a-b116a3d0eb4a
# ╟─f7483843-9cc9-4799-ad56-1e74bd0a08e7
# ╠═b985704f-8f07-43b6-baa9-6b3d4541de27
# ╠═31a1b2aa-7b5b-46b0-82ec-f75111b72a0d
# ╠═a0e4f6af-d4ab-4b04-9863-49de95eec947
# ╠═c20593d7-885f-4f36-984c-ecded0884ecc
# ╟─cef37fa4-df23-4267-bde4-c1fe6ddade23
# ╠═92d27828-c97b-41d0-85e9-0d3b1738dfed
# ╠═6055ba96-dd58-420b-854e-277b60a3c1e7
# ╠═ca176d00-70f9-42b4-a3f4-fd5d37c9f31a
# ╠═a8080f5e-6797-482e-84ca-7ddfcd40264d
# ╠═d4764a31-38b7-4002-92fb-bfd3580cf634
# ╠═2a644f10-2683-4ea8-bd25-95343abb7e48
# ╠═ca57bf59-f348-439a-8622-ea0edcc86f17
# ╠═bc477197-2ac5-4798-a62b-62f89791fb94
# ╠═8ff25700-5584-11eb-1f95-abf9735d479e
# ╠═90043461-3f6a-45e6-a661-57c8df3a2a93
# ╠═a107fc1f-f0db-42cb-ad4a-ffd51726808a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
