### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 0d4fbf25-a629-416a-bb67-b087144a5862
using Plots

# ╔═╡ 69cef7ae-8c14-46fc-8583-942c1bdec167
using ODE

# ╔═╡ b5d87000-5af5-11eb-0594-57573e58e809
md"""
# Diferencijalne jednadžbe


Riješimo diferencijalnu jednadžbu

$$
\frac{d}{dx} y(x)=f(x,y(x)),$$

uz zadani __početni uvjet__ 

$$
y(x_0)=y_0.$$

__Teorem.__ Ako su funkcije $f$ i $\displaystyle \frac{\partial f}{\partial y}$ neprekidne u nekoj okolini točke $(x_0,y_0)$, onda  za neki $\varepsilon > 0$ zadani problem početnih vrijednosti ima jedinstveno rješenje $y(x)$ na intervalu $\displaystyle [x_{0}-\varepsilon ,x_{0}+\varepsilon ]$.

__Napomena.__ Nezavisnu varijablu često označavamo s $t$ (vrijeme). 
"""

# ╔═╡ 68a29859-4ba6-46ed-97ef-fdffec416003
md"""
## Eulerova metoda

Taylorovu formulu u okolini točke $x$ možemo zapisati kao (uz pretpostavku da je $y'''(x)$ omeđena funkcija) 

$$
y(x+h)=y(x)+h\cdot y'(x) +\frac{1}{2}h^2 \cdot y''(x)+ O(h^3). \tag{1}$$


Počevši od točke $x_0$ i početnog uvjeta, za niz jednako udaljenih točaka

$$
x_{k+1}=x_{k}+h,$$

__Eulerova metoda__ aproksimira vrijednost funkcije $y$ u točki s prva dva člana Taylorovog reda oko točke $x_k$:

$$
y_{k+1}=y_k+h \cdot f(x_k,y_k). \tag{2}$$

Ako je $y(x_k)=y_k$, onda je __lokalna pogreška__ metode pogreška u jednom koraku Taylorove formule, odnosno vrijedi 

$$
y(x_{k+1})-y_{k+1}=\displaystyle\frac{1}{2}h^2 \cdot y''(x_k)+ O(h^3).$$ 

__Uvjetovanost.__ Tijekom računanja je $y(x_k)\neq y_k$. Prema Lagrangeovom teoremu srednje vrijednosti je

$$
f(x_k,y_k)=f(x_k,y(x_k))+f_y(x_k,\eta)(y_k-y(x_k)) \tag{3}$$

za neki $\eta$ koji se nalazi između $y_k$ i $y(x_k)$. Koristeći Taylorovu formulu (1) za $x=x_k$ i Eulerovu formulu (2) imamo

$$
y(x_{k+1})-y_{k+1}=y(x_k)+h\cdot f(x_k,y(x_k)) +\frac{1}{2}h^2 \cdot y''(x_k)+ O(h^3)
-y_k-h \cdot f(x_k,y_k).$$

Uvrštavanje Lagrangeove formule (3) daje

$$
y(x_{k+1})-y_{k+1}=y(x_k)-y_k + h\cdot [f(x_k,y_k)-f_y(x_k,\eta)(y_k-y(x_k))] - h \cdot f(x_k,y_k)+O(h^2),$$

odnosno 

$$
y(x_{k+1})-y_{k+1}=(y(x_k)-y_k) [1+h\cdot f_y(x_k,\eta)]+O(h^2). \tag{4}$$

Zaključujemo da tijekom propagacije pogreška pada ako je $f_y<0$, a raste ako je $f_y>0$, pa je u tom slučaju problem loše uvjetovan.


__Globalna pogreška__ je pogreška u točki $x$, nakon $n$ koraka potrebnih da metoda od početne točke dođe do točke $x$. Globalna pogreška nastaje radi toga što je $y_k$ netočan zbog akumulacije prethodnih lokalnih pogrešaka te zbog pogreške Taylorove formule u $n$-tom koraku. Ako je prirast $h$, onda je očito $n=\displaystyle\frac{x-x_0}{h}$ pa očekujemo da će globalna pogreška biti (približno) proporcionalna s $h$. To možemo zaključiti i iz ocjene (4).
Stoga kažemo da je Eulerove metoda __metoda prvog reda__. 
"""

# ╔═╡ fce92627-56f1-482e-983a-083b7281fbcb
function Euler(f::Function,y₀::T,x::T1) where {T,T1}
    h=x[2]-x[1]
    y=Array{T}(undef,length(x))
    y[1]=y₀
    for k=2:length(x)
        y[k]=y[k-1]+h*f(x[k-1],y[k-1])
    end
    y
end

# ╔═╡ 934a054a-f31d-497e-bcf8-295dac6db83b
md"""
### Primjer 1

Rješenje problema početnih vrijednosti 

$$
y'=x+y,\quad y(0)=1,$$

je

$$
y=2e^x-x-1$$

(vidi [Numerička matematika, primjer 8.1](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).
"""

# ╔═╡ 4f30e023-2c8d-42b7-a124-c240e5464769
begin
	# 10 pod-intervala na intervalu [0,1]
	x₁=range(0,stop=1,length=11)
	f₁(x,y)=x+y
	y₁=Euler(f₁,1.0,x₁)
end

# ╔═╡ 6f9f8bbd-1318-4d66-bdd6-6d02f778f8b8
begin
	# Nacrtajmo točno rješenje i izračunate točke
	solution₁(x)=2*exp(x)-x-1
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Točno rješenje",legend=:topleft)
	plot!(x₁,y₁,label="Euler()")
	scatter!(x₁,y₁,label="Izračunate točke")
end

# ╔═╡ 5e6c51c1-c7a9-4542-89fd-9c7b0846437c
md"""
Uočimo da je zadani problem loše uvjetovan.
Točnu ocjenu globalne pogreške daje nam sljedeći teorem:

__Teorem.__ Neka funkcija $f(x,y)$ ima neprekidne prve parcijalne derivacije na pravokutniku $D\subseteq \mathbb{R}^2$ i neka je 

$$
K=\max_{(x,y)\in D}|f_y(x, y)|<\infty,\quad   M= \max_{(x,y)\in D}|(f_x+f\cdot f_y)(x,y)|<\infty.$$

Ako je 

$$(x_k,y_k)\in D\quad (x_k,y(x_k))\in D, \quad k=0,1,2,\ldots n,$$

i ako je $Kh<1$, onda je pogreška omeđena s

$$
|y(x)-y_n|\leq \frac{M}{2K}\left(e^{Kx}−1\right)h. \tag{5}$$

_Dokaz._ Vidi [Glenn Ledder, Error in Euler’s Method](http://www.math.unl.edu/~gledder1/Math447/EulerError).
"""

# ╔═╡ 54c26825-46b9-44f6-ae65-31f4be69bf3b
md"""
__Zadatak.__ Odredimo broj koraka $n$ da bi u Primjeru 1  izračunata vrijednost $y(1)$ imala pogrešku manju od $\epsilon=0.01$. 

Vrijedi $f_y(x,y)=1$ pa je $K=1$. Za pozitivne $x$ i $y$ je $y'>0$ pa je $y$ rastuća funkcija. Vrijednost $y(1)$ ćemo procijeniti na $4$ jer je u prethodnom računu s $10$ pod-intervala $y_{10}\approx 3.19$. Stoga je $M=5$. Uvrštavajući $K$, $M$ i $h=\displaystyle\frac{1-0}{n}$ u formulu (5), imamo

$$
\frac{5}{2}(2.7183-1)\frac{1}{n} \approx 4.3 \frac{1}{n} < 0.01,$$

pa možemo uzeti $n=500$.
"""

# ╔═╡ 22f0b250-0b37-4593-bdea-4bd13c19f2bf
begin
	# 500 pod-intervala na intervalu [0,1]
	xx₁=range(0,stop=1,length=501)
	yy₁=Euler(f₁,1.0,xx₁)
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Točno rješenje",legend=:topleft)
	plot!(xx₁,yy₁,label="Euler()")
end

# ╔═╡ bfc718ab-e392-4fb8-82c6-0c889b1977e9
# Provjera točnosti
solution₁(1)-yy₁[end]

# ╔═╡ 12ceb8e5-f976-4639-9175-047e84321c85
md"""
### Primjer 2

Rješenje problema

$$
y'=30(\sin x-y), \quad y(0)=0,$$

je

$$
y(x)=\frac{30}{901}(30\sin x-\cos x+e^{-30x})$$

(vidi [Numerička matematika, primjer 8.3](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).
"""

# ╔═╡ 1e36310b-6f42-44a1-8b04-e82391cb33aa
begin
	# 100 pod-intervala na intervalu [0,1]
	f₂(x,y)=30(sin(x)-y)
	x₂=range(0,stop=1,length=101)
	y₂=Euler(f₂,0.0,x₂)
	solution₂(x)=30(30*sin(x)-cos(x)+exp(-30x))/901
	plot(solution₂,0,1,xlabel="x",ylabel="y",
	    label="Točno rješenje",legend=:topleft)
	plot!(x₂,y₂,label="Euler()")
end

# ╔═╡ 8899d0fd-e0d9-424f-98e4-5547d029a36a
md"""
__Zadatak.__ Procijenite točnost izračunate vrijednosti $y(1)$  koristeći ocjenu (5).
"""

# ╔═╡ 2d44d600-5af7-11eb-0958-b3eab7ae938d
md"""
## Metode Runge-Kutta

Eulerova metoda je metoda prvog reda pa nema zadovoljavajuću točnost. Stoga se u praksi koriste metode višeg reda koje vrijednost funkcije $y(x)$ u točki $x_{k+1}$ aproksimiraju pomoću vrijednosti funkcije $f(x,y)$ u nekoliko odabranih točaka na intervalu 

$$[x_k,x_{k+1}]\equiv[x_k,x_k+h].$$

### Heun-ova metoda

$$\begin{aligned}
k_1&=hf(x_k,y_k),\\
k_2&=hf(x_k+h,y_k+k_1),\\
y_{k+1}&=y_k+\frac{1}{2}(k_1+k_2).
\end{aligned}$$

Heunovu metodu izvodimo na sljedeći način: integriranje zadane diferencijalne jednadžbe daje

$$
\int_{x_k}^{x_{k+1}} \frac{dy}{dx}\, dx = y(x_{k+1})-y(x_k)= \int_{x_k}^{x_{k+1}} 
f(x,y(x))\, dx.$$

Trapezna formula daje

$$
y(x_{k+1})= y(x_k) +\frac{h}{2}[f(x_k,y(x_k))+f(x_{k+1},y(x_{k+1}))].$$

Uz oznake $y(x_k)=y_k$, $y(x_{k+1})=y_{k+1}$, primjena Eulerove formule (2) na desnoj strani daje Heunove formule. 

Trapezna formula ima pogrešku reda veličine $O(h^2)$ pa je i globalna pogreška Heunove metode također reda veličine $O(h^2)$. Potrebne su dvije evaluacije funkcije $f(x,y)$ u svakom koraku.
"""

# ╔═╡ 8b933ce6-7fea-4190-940a-bd7720fc93be
md"""

### Standardna Runge-Kutta metoda

$$\begin{aligned}
k_1&=hf(x_k,y_k),\\
k_2&=hf\big(x_k+\frac{h}{2},y_k+\frac{k_1}{2}\big),\\
k_3&=hf\big(x_k+\frac{h}{2},y_k+\frac{k_2}{2}\big),\\
k_4&=hf(x_k+h,y_k+k_3),\\
y_{k+1}&=y_k+\frac{1}{6}(k_1+2k_2+2k_3+k_4).
\end{aligned}$$

Standardna Runge-Kutta metoda je metoda reda 4 - globalna pogreška je reda veličine $O(h^4)$, a potrebne su četiri evaluacije funkcije $f(x,y)$ u svakom koraku.

Lokalna pogreška standardne Runge-Kutta metode je $O(h^5)$.
"""

# ╔═╡ cb8721aa-7f4b-4672-a802-da84d9c67b4a
function RungeKutta4(f::Function,y₀::T,x::T1) where {T,T1}
    h=x[2]-x[1]
    y=Array{T}(undef,length(x))
    y[1]=y₀
    for k=2:length(x)
        ξ=x[k-1]
        η=y[k-1]
        k₁=h*f(ξ,η)
        k₂=h*f(ξ+h/2,η+k₁/2)
        k₃=h*f(ξ+h/2,η+k₂/2)
        k₄=h*f(ξ+h,η+k₃)
        y[k]=η+(k₁+2*k₂+2*k₃+k₄)/6.0
    end
    y
end

# ╔═╡ ad54ce7f-8caa-46f8-9a57-f3710ca89f40
md"""
### Primjer 3

Riješimo probleme iz Primjera 1 i 2. Za Primjer 1 se numeričko rješenje grafički preklapa s točnim rješenjem. 
Za Primjer 2 je rješenje pomoću funkcije `RungeKutta4()` za red veličine točnije od rješenja dobivenog pomoću funkcije `Euler()`. 
"""

# ╔═╡ 4adf6672-5ae1-42e1-85e9-e522d6069fab
begin
	y₃=RungeKutta4(f₁,1.0,x₁)
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Točno rješenje",legend=:topleft)
	plot!(x₁,y₃,label="RunkeKutta4()")
end

# ╔═╡ a262689c-74f3-46bb-bd55-1aa7eec5379f
begin
	x₄=range(0,stop=1,length=21)
	yEuler=Euler(f₂,0.0,x₄)
	yRK4=RungeKutta4(f₂,0.0,x₄)
	plot(solution₂,0,1,xlabel="x",ylabel="y",
	    label="Točno rješenje",legend=:topleft)
	plot!(x₄,[yEuler,yRK4],label=["Euler()" "RungeKutta4()"])
end

# ╔═╡ 481e7da8-95f2-4972-863c-fcd7a79cf417
solution₂(1), yEuler[end],yRK4[end]

# ╔═╡ d52999d7-2644-4793-a36a-413ad8011e09
md"""
### Postojeće rutine

Većina programa ima ugrađene odgovarajuće rutine za numeričko rješavanje običnih diferencijalnih jednadžbi.
Tako, na primjer, 

* Ṁatlab ima naredbe `ode*` (vidi [Matlab, Ordinary Differential Equations](https://www.mathworks.com/help/matlab/ordinary-differential-equations.html?searchHighlight=ordinary%20differential&s_tid=srchtitle)), a 
* Julia ima paket  [ODE.jl](https://github.com/JuliaODE/ODE.jl).

Standardna Runge-Kutta metoda je implementirana u funkciji `ode4()`, a Heunova metoda je implementirana u funkciji 
`ODE.ode2_heun()`. 

__Napomena__ Funkcija `ODE.ode2_heun()` nije vidljiva pozivom naredbe `varinfo()` jer nije izvezena, ali se može vidjeti u datoteci `runge_kutta.jl`.
"""

# ╔═╡ 1a0815c9-900d-4ba5-855a-7fa582e8d498
# Pogledajte!
# varinfo(ODE)

# ╔═╡ d8e15baa-34a9-447d-90d2-8dfe8065ee1b
methods(ode4)

# ╔═╡ f35f1bb2-bc2c-4bb7-9dc9-367bc5a6c15d
methods(ODE.ode2_heun)

# ╔═╡ e50cf9c2-972d-4916-a3ca-6c20a0c73519
# Riješimo problem iz Primjera 2.
# Tražene vrijednosti funkcije y su drugi element izlaza.
yode4=ode4(f₂,0.0,range(0,stop=1,length=21))[2]

# ╔═╡ f4d37a14-e6ef-4400-bb30-a3a4405bb675
yode2=ODE.ode2_heun(f₂,0.0,range(0,stop=1,length=21))[2];

# ╔═╡ a7841a3d-c7e4-4c09-9012-a3890c9e1195
# Usporedimo rješenja
[yRK4 yode4 yRK4-yode4 yode2 yode4-yode2]

# ╔═╡ 471ddc14-d7ba-4d94-a902-9ca459be5ce3
md"""
## Sustavi diferencijalnih jednadžbi

__Problem.__ Riješimo sustav od $n$ jednadžbi 

$$\begin{aligned}
y_1'(x)&=f_1(x,y_1,y_2,\ldots,y_n),\\
y_2'(x)&=f_2(x,y_1,y_2,\ldots,y_n),\\
&\vdots \\
y_n'(x)&=f_2(x,y_1,y_2,\ldots,y_n)
\end{aligned}$$

i $n$ nepoznatih funkcija $y_1,y_2,\ldots,y_n$ uz početne uvjete

$$
y_i(x_0)=\zeta_i.$$

Uz oznake 

$$
f=\begin{bmatrix}f_1\\ f_2\\ \vdots\\f_n\end{bmatrix},\quad
y=\begin{bmatrix}y_1\\ y_2\\ \vdots \\y_n\end{bmatrix},\quad
\zeta=\begin{bmatrix}\zeta_1\\ \zeta_2\\ \vdots\\ \zeta_n\end{bmatrix},$$

zadani problem možemo zapisati u vektorskom obliku kao 

$$
y'(x)=f(x,y),\quad y(x_0)=\zeta.$$

Problem se uspješno rješava Eulerovom metodom i Runge-Kutta metodama u vektorskom obliku.
"""

# ╔═╡ 705ebcb9-bfa6-4357-8610-27c43c0c3f96
md"""
### Lhotka-Volterra jednadžbe

Modeliranje sustava __lovac-plijen__ daje sustav __Lhotka-Volterra__ jednadžbi
(vidi [Matematika 2, poglavlje 5.11](http://lavica.fesb.hr/mat2/predavanja/node98.html) i [Numerička matematika, primjer 8.7](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)):

$$\begin{aligned}
\frac{dZ}{dt}&=z\,Z-a\, Z\, V = Z\,(z-a\, V), \qquad(6)\\
\frac{dV}{dt}&=-v\,V+b\, Z\, V = V\,(-v+b\, Z), \quad v,z,a,b>0,
\end{aligned}$$

uz početne uvjete 

$$
V(t_0)=V_0,\qquad Z(t_0)=Z_0.$$

__Stabilna stanja__ su stanja u kojima nema promjene, odnosno, stanja u kojima su obje derivacije jednake nuli. To su __trivijalno__ stabilno stanje, $V=Z=0$, i

$$
V=\frac{z}{a},\qquad Z=\frac{b}{v}.$$


U __faznom prostoru__, eliminacijom nezavisne varijable $t$, dobijemo jednu diferencijalnu jednadžbu prvog reda:

$$
\frac{dV}{dZ}=\frac{\displaystyle\frac{dV}{dt}}{\displaystyle\frac{dZ}{dt}}=
\frac{V\,(-v+b\, Z)}{Z\,(z-a\, V)}.$$

Ovo je jednadžba sa separiranim varijablama koja ima implicitno rješenje

$$
V^z \, Z^v = C\,  e^{aV}\, e^{bZ},\qquad 
C=\frac{V_0^{z}\, Z_0^{v}}{\displaystyle e^{a V_0}e^{b Z_0}}.  \tag{7} $$

Riješimo sustav za populacije vukova $V$ i zečeva $Z$ uz

$$
v=0.02, \quad z=0.06,\quad a=0.001,\quad b=0.00002,\quad V(0)=30, \quad Z(0)=800.$$

__Napomena.__ Funkcije `Euler()`, `RungeKutta4()` i funkcije iz paketa `ODE.jl` su već prilagođene i za rješavanje sustava.
"""

# ╔═╡ a6cd2db9-f4b6-46cb-aaa0-59d53430056c
begin
	# Zadani problem: y=[V,Z], t₀=0, y₀=y(0)=[30,800]
	# 356 dana s razmakom od 1/10 dana
	t=range(0,stop=365,length=3651)
	# Parametri sustava
	v=0.02
	z=0.06
	a=0.001
	b=0.00002
	# Početna populacija vukova
	V₀=30.0
	# Početna populacija zečeva
	Z₀=800.0
	# Početna točka
	y₀=[V₀,Z₀]
	# Vektorska funkcija
	fVZ(t,y)=[y[1]*(-v+b*y[2]),y[2]*(z-a*y[1])]
	# Rješenja
	yₗ=Euler(fVZ,y₀,t)
end

# ╔═╡ 78eae1ea-5ce9-4bda-81d3-4b056ee7df1e
begin
	# Skaliramo Z u Z/10 da graf bude čitkiji
	V=[yₗ[i][1] for i=1:length(yₗ)]
	Z=[yₗ[i][2] for i=1:length(yₗ)]
	plot(t,[V,Z/10],xlabel="t ( dani )", ylabel="Populacija", 
	    label=["Vukovi" "Zečevi / 10"])
end

# ╔═╡ 00796b0e-f2d5-4646-af33-06d81d13fc0c
begin
	# Usporedimo rješenja pomoću metoda RungeKutta4() i ode4() 
	yRK4ₗ=RungeKutta4(fVZ,y₀,t)
	yode4ₗ=ode4(fVZ,y₀,t)[2]
	[yₗ yRK4ₗ yode4ₗ]
end

# ╔═╡ 0b5ad338-3628-477e-bd0b-8cdda13d345d
# Nacrtajmo rješenje u faznom prostoru
plot(Z,V,xlabel="Z",ylabel="V", label=:none)

# ╔═╡ 5308aa28-6da6-4943-b393-c08592b98c27
md"""
### Skalirane Lhotka-Volterra jednadžbe

Crtanje egzaktnog rješenja (7) u faznom prostoru nije praktično, 
jer crtanje implicitno zadanih funkcija na velikom području traje dugo. Međutim, 
pomoću transformacija (vidi [Modeling Complex Systems, poglavlje 2.1](http://www.springer.com/us/book/9781441965615))

$$
X=\frac{b}{v}Z,\quad Y=\frac{a}{z}V,\quad \tau=\sqrt{z\cdot v}\, t,\quad \rho=\sqrt{\displaystyle\frac{z}{v}},$$

jednadžbu (6) je moguće prikazati s __bezdimenzionalnim varijablama__
u __skaliranom vremenu__ $\tau$:

$$\begin{aligned}
\frac{dX}{d\tau}&=\rho\, X\,(1-Y),\\
\frac{dY}{d\tau}&=-\frac{1}{\rho}\, Y\,(1-X). 
\end{aligned}$$

Ovaj sustav ovisi o samo __jednom__ parametru $\rho$. Sustav ima netrivijalno stabilno rješenje $X=Y=1$, a rješenje (7) u faznom prostoru je

$$
Y \, X^{1/\rho^2} = C\,  e^{Y}\, e^{X/\rho^2},\qquad 
C=\frac{Y_0\, X_0^{1/\rho^2}}{\displaystyle e^{Y_0}e^{X_0/\rho^2}}.$$


Riješimo zadani sustav u bezdimenzionalnom obliku i grafički usporedimo rješenja:

"""

# ╔═╡ 00f084ce-a5ff-4e18-913f-9d14ae346bd9
begin
	ρ=sqrt(z/v)
	τ=range(0,stop=365*sqrt(z*v),length=3651)
	Y₀=[Z₀*b/v,V₀*a/z]
	fXY(τ,y)=[ρ*y[1]*(1-y[2]),-y[2]*(1-y[1])/ρ]
	Y=Euler(fXY,Y₀,τ);
end

# ╔═╡ 0ab04df9-d9a0-406a-b45f-299a0212fd9b
begin
	X₁=[Y[i][1] for i=1:length(Y)]
	Y₁=[Y[i][2] for i=1:length(Y)]
	# Rješenja se poklapaju
	plot(t,[V,Z/10],xlabel="t ( dani )", ylabel="Populacija", 
	    label=["Vukovi" "Zečevi / 10"])
	plot!(τ/sqrt(z*v),[Y₁*z/a,X₁*v/(10b)],
	    label=["Vukovi (bezdim.)" "Zečevi / 10 (bezdim.)"])
end

# ╔═╡ e580ec4b-73f9-4522-9cb3-4d24ebe1e027
# Nacrtajmo bezdimenzionalno rješenje u faznom prostoru
plot(X₁,Y₁,xlabel="Z",ylabel="V", label=:none)

# ╔═╡ 3ada3bbf-cc3b-46b4-adbd-c3685568781a
md"""
## Diferencijalne jednadžbe višeg reda

Diferencijalna jednadžba višeg reda supstitucijama se može svesti na sustav diferencijalnih jednadžbi prvog reda.

### Primjer 4

Rješenje problema početnih vrijednosti (vidi [Matematika 2, primjer 5.28](http://lavica.fesb.hr/mat2/predavanja/node97.html#pr:ld2khv))

$$
y'''+y''=x,\qquad y(0)=0,\quad y'(0)=0,\quad y''(0)=0$$

je

$$
y(x)=-1+x+e^{-x}+\frac{x^3}{6}-\frac{x^2}{2}.$$

Supstitucije

$$
y'=u,\quad y''=v,$$

daju sustav 

$$\begin{aligned}
y'&=u \\
u'&=v \\
v'&=-v+x\\
\end{aligned}$$

uz početne uvjete

$$
y(0)=0,\quad u(0)=0,\quad v(0)=0.$$

"""

# ╔═╡ 39d0b178-70bd-4815-bee8-3738bdaea250
begin
	x₅=range(0,stop=2,length=201)
	# Vektorska funkcija
	f₅(x,y)=[y[2],y[3],-y[3]+x]
	# Početni uvjeti
	y₅=[0.0,0,0]
	# Riješimo sustav pomoću tri metode
	yEuler₅=Euler(f₅,y₅,x₅)
	yRK4₅=RungeKutta4(f₅,y₅,x₅)
	yode4₅=ode4(f₅,y₅,x₅)
	# Izračunate vrijednosti yₖ su prvi elementi izlaznih vektora
	YEuler₅=[yEuler₅[i][1] for i=1:length(yEuler₅)]
	YRK4₅=[yRK4₅[i][1] for i=1:length(yEuler₅)]
	# Izračunate vrijednosti yₖ su prvi elementi vektora
	# drugog izlaznog polja
	Yode4₅=[yode4₅[2][i][1] for i=1:length(yEuler₅)]
	# Točno rješenje
	solution₅(x)=-1+x+exp(-x)+x^3/6-x^2/2
	# Crtanje
	plot(solution₅,0,2,xlabel="x",ylabel="y",label="Točno rješenje",legend=:topleft)
	plot!(x₅,[YEuler₅,YRK4₅,Yode4₅],xlabel="x",ylabel="y",label=["Euler()" "RungeKutta4()" "ode4()"])
end

# ╔═╡ 12c6cb02-6bc0-4558-bc10-9ca1d8cd029a
# Norma pogreške u promatranim točkama
# import LinearAlgebra; LinearAlgebra.norm(solution₅.(x)-Y)
sqrt(sum((solution₅.(x₅)-YRK4₅).^2))

# ╔═╡ ae54d9f1-bed6-4f56-83a2-832e819b5068


# ╔═╡ Cell order:
# ╟─b5d87000-5af5-11eb-0594-57573e58e809
# ╟─68a29859-4ba6-46ed-97ef-fdffec416003
# ╠═0d4fbf25-a629-416a-bb67-b087144a5862
# ╠═fce92627-56f1-482e-983a-083b7281fbcb
# ╟─934a054a-f31d-497e-bcf8-295dac6db83b
# ╠═4f30e023-2c8d-42b7-a124-c240e5464769
# ╠═6f9f8bbd-1318-4d66-bdd6-6d02f778f8b8
# ╟─5e6c51c1-c7a9-4542-89fd-9c7b0846437c
# ╟─54c26825-46b9-44f6-ae65-31f4be69bf3b
# ╠═22f0b250-0b37-4593-bdea-4bd13c19f2bf
# ╠═bfc718ab-e392-4fb8-82c6-0c889b1977e9
# ╟─12ceb8e5-f976-4639-9175-047e84321c85
# ╠═1e36310b-6f42-44a1-8b04-e82391cb33aa
# ╟─8899d0fd-e0d9-424f-98e4-5547d029a36a
# ╟─2d44d600-5af7-11eb-0958-b3eab7ae938d
# ╟─8b933ce6-7fea-4190-940a-bd7720fc93be
# ╠═cb8721aa-7f4b-4672-a802-da84d9c67b4a
# ╟─ad54ce7f-8caa-46f8-9a57-f3710ca89f40
# ╠═4adf6672-5ae1-42e1-85e9-e522d6069fab
# ╠═a262689c-74f3-46bb-bd55-1aa7eec5379f
# ╠═481e7da8-95f2-4972-863c-fcd7a79cf417
# ╟─d52999d7-2644-4793-a36a-413ad8011e09
# ╠═69cef7ae-8c14-46fc-8583-942c1bdec167
# ╠═1a0815c9-900d-4ba5-855a-7fa582e8d498
# ╠═d8e15baa-34a9-447d-90d2-8dfe8065ee1b
# ╠═f35f1bb2-bc2c-4bb7-9dc9-367bc5a6c15d
# ╠═e50cf9c2-972d-4916-a3ca-6c20a0c73519
# ╠═f4d37a14-e6ef-4400-bb30-a3a4405bb675
# ╠═a7841a3d-c7e4-4c09-9012-a3890c9e1195
# ╟─471ddc14-d7ba-4d94-a902-9ca459be5ce3
# ╟─705ebcb9-bfa6-4357-8610-27c43c0c3f96
# ╠═a6cd2db9-f4b6-46cb-aaa0-59d53430056c
# ╠═78eae1ea-5ce9-4bda-81d3-4b056ee7df1e
# ╠═00796b0e-f2d5-4646-af33-06d81d13fc0c
# ╠═0b5ad338-3628-477e-bd0b-8cdda13d345d
# ╟─5308aa28-6da6-4943-b393-c08592b98c27
# ╠═00f084ce-a5ff-4e18-913f-9d14ae346bd9
# ╠═0ab04df9-d9a0-406a-b45f-299a0212fd9b
# ╠═e580ec4b-73f9-4522-9cb3-4d24ebe1e027
# ╟─3ada3bbf-cc3b-46b4-adbd-c3685568781a
# ╠═39d0b178-70bd-4815-bee8-3738bdaea250
# ╠═12c6cb02-6bc0-4558-bc10-9ca1d8cd029a
# ╠═ae54d9f1-bed6-4f56-83a2-832e819b5068