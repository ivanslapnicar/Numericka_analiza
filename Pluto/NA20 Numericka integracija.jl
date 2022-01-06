### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 7d2eceff-5dc8-4088-8504-8fa469948983
using PlutoUI, Plots, QuadGK, FastGaussQuadrature, LinearAlgebra

# ╔═╡ 456ae14f-d92c-4561-bcfd-d1b9849b8244
using FFTW

# ╔═╡ 3c1c274f-a299-4f97-a94f-5ac190e00f76
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ c42b0321-1003-4935-857e-a8422fabd485
md"""

# Numerička integracija

# Newton-Cotesove formule

Funkciju $f(x):[a,b]\to\mathbb{R}$ interpoliramo polinomom stupnja $n$ kroz $n+1$ ravnomjerno raspoređenu točku te integral aproksimiramo integralom interpolacijskog polinoma. Polinom $P_n(x)$ možemo računati u Lagrangeovom obliku  
(vidi bilježnicu [NA09 Interpolacijski polinomi.jl](https://ivanslapnicar.github.io/Numericka_analiza/NA09%20Interpolacijski%20polinomi.jl)):

$$
L_k(x)=\prod_{{i=0}\atop {i\neq k}}^n \frac{x-x_i}{x_k-x_i},$$

Tada je 

$$
f(x)\approx P_n(x)=\sum_{k=0}^n f(x_k) L_k(x),$$

pa je 

$$
\int_a^b f(x)\, dx\approx \int_a^b P_n(x) \, dx=\sum_{k=0}^n f(x_k) \int_a^b L_k(x)\, dx =(b-a)\sum_{k=0}^n \omega_k f(x_k). \qquad (1)$$

Vrijedi 

$$x_i=a+\displaystyle\frac{b-a}{n}\, i.$$

Uz supstituciju $x=a+(b-a)\,t$ vrijedi

$$
\frac{x-x_i}{x_k-x_i}=\frac{nt-i}{k-i}$$

pa su __težine__ $\omega_k$ jednake

$$
\omega_k=\frac{1}{b-a}\int_a^b L_k(x)\, dx = \int_0^1 \prod_{{i=0}\atop {i\neq k}}^n \frac{nt-i}{k-i} \, dt. \qquad(2)$$
"""

# ╔═╡ 5071e7f3-d353-4429-8885-ca9431edf5a9
md"""
## Trapezna formula

Za $n=1$ formula (2) daje 

$$\omega_0=\omega_1=\frac{1}{2}$$

pa Newton-Cotesova formula (1) daje

$$\int_a^b f(x)\, dx\approx \int_a^b P_1(x) \, dx=\frac{b-a}{2}(f(a)+f(b)),$$

što je površina trapeza s vrhovima $(a,0)$, $(a,f(a))$, $b,f(b))$ i $(b,0)$.

Točno značenje $\approx$ daje nam sljedeći teorem:

__Teorem.__ Ako je funkcija $f$ neprekidna na intervalu $[a,b]$, onda postoji $c\in(a,b)$ takav da je

$$\int_a^b f(x)\, dx= \frac{b-a}{2}(f(a)+f(b))-\frac{(b-a)^3}{12}f''(c).$$

_Dokaz:_ Prema Teoremu iz bilježnice [NA10 Interpolacija funkcija.jl](https://ivanslapnicar.github.io/Numericka_analiza/NA10%20Interpolacija%20funkcija.jl.html) za $c\in(a,b)$ vrijedi 

$$
E=\int_a^b f(x)\, dx -\int_a^b P_1(x) \, dx = \int_a^b (f(x)- P_1(x)) \, dx = \int_a^b \frac{f''(c)}{2}(x-a)(x-b)\, dx.$$

Uz supstituciju $x=a+(b-a)\,t$ vrijedi

$$
E=\frac{f''(c)}{2} (b-a)^3 \int_0^1 t(t-1)\, dt= -\frac{(b-a)^3}{12}f''(c)$$

i teorem je dokazan. 

Sada interval $[a,b]$ podjelimo na $n$ jednakih podintervala,

$$[x_{i-1},x_{i}],\quad  i=1,2,\ldots,n,$$ 

i uvedimo oznake

$$
\Delta x=\frac{b-a}{n}, \quad y_i=f(x_i).$$


Na svakom podintervalu vrijedi 

$$\int_{x_{i-1}}^{x_i} f(x)\, dx= \frac{\Delta x}{2}(y_{i-1}+y_i)
-\frac{(\Delta x)^3}{12}f''(c_i),\quad c_i\in[x_{i-1},x_i].$$

Zbrajanje daje

$$
\int_a^b f(x)\, dx=I_n+E_n,$$ 

gdje je $I_n$ __trapezna formula__,

$$
I_n=\Delta x\bigg( \frac{y_0}{2} +y_1+y_2+\cdots +y_{n-1}+\frac{y_n}{2}\bigg),$$

a za __pogrešku__ $E_n$ vrijedi

$$
E_n =-\frac{(\Delta x)^3}{12}\sum_{i=1}^n f''(c_i)=
-\frac{b-a}{12}(\Delta x)^2 f''(c), \quad c\in[a,b].$$

Koristili smo činjenicu da zbog neprekidnosti funkcije $f''$ na intervalu $(a,b)$ postoji točka $c\in[a,b]$ za koju vrijedi 

$$
\frac{1}{n}\sum_{i=1}^n f''(c_i)=f''(c).$$

Također, vrijedi 

$$
|E_n|\leq \frac{b-a}{12}(\Delta x)^2 \max_{x\in(a,b)} |f''(x)|. \tag{3}$$


__Napomena__. Izvod trapezne formule i ocjene pogreške dan je u knjigama
[Numerička matematika, poglavlje 7.1](http://www.mathos.unios.hr/pim/Materijali/Num.pdf) i [Matematika 2, poglavlje 2.7.2](http://lavica.fesb.hr/mat2/predavanja/node42.html).

"""

# ╔═╡ 04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
md"""
## Simpsonova formula

Za $n=2$ Newton-Cotesova formula (1) daje

$$\omega_0=\frac{1}{6},\quad \omega_1=\frac{2}{3},\quad \omega_2=\frac{1}{6}.$$

Interval $[a,b]$ podjelimo na paran broj $n$ jednakih podintervala, na svaki podinterval

$$[x_{2i-1},x_{2i+1}],\quad i=1,2,\ldots,\frac{n}{2},$$

primijenimo  Newton-Cotesovu formulu i zbrojimo, što daje __Simpsonovu formulu__:

$$
I_n=\frac{\Delta x}{3}\big(y_0 +4(y_1+y_3\cdots +y_{n-1})+2(y_2+y_4+\cdots+y_{n-2})+y_n\big).$$

Vrijedi

$$
\int_a^b f(x)\, dx =I_n+E_n,$$

pri čemu je __pogreška__ $E_n$ omeđena s

$$
|E_n|\leq \frac{b-a}{180}(\Delta x)^4 \max_{x\in(a,b)} |f^{(4)}(x)|.\tag{4}$$

Za detalje vidi knjige [Numerička matematika, poglavlje 7.3](http://www.mathos.unios.hr/pim/Materijali/Num.pdf) i [Matematika 2, poglavlje 2.7.3](http://lavica.fesb.hr/mat2/predavanja/node43.html).

"""

# ╔═╡ 7a11ccd0-68c1-46d3-9249-3cc9f867ef65
md"""
# Richardsonova ekstrapolacija

Ocjena pogreške pomoću formula (3) i (4) može biti složena. __Richardsonova ekstrapolacija__ nam omogućava da, uz određene uvjete, pogrešku procijenimo koristeći aproksimaciju integrala s $n/2$ točaka.
Ako se u ocjeni pogreške javlja član $(\Delta x)^m$, onda je pogreška približno manja od broja (vidi [Matematika 2, poglavlje 2.7.4](http://lavica.fesb.hr/mat2/predavanja/node44.html))

$$
E=\frac{\big(\frac{n}{2}\big)^m}{n^m-\big(\frac{n}{2}\big)^m}(I_n-I_{n/2}).$$

Predznak broja $E$ daje i predznak pogreške, odnosno, ako je $E>0$, tada je približno

$$
\int_a^b f(x)\, dx\in[I_n,I_n+E],$$

a ako je $E\leq 0$, tada je približno

$$
\int_a^b f(x)\, dx\in[I_n+E,I_n].$$

Kod dokazivanja ovih formula pretpostavljamo da postoji broj $\omega$ takav da je 

$$
\displaystyle I=I_n+E_n, \qquad E_n=\omega\,  (\Delta x)^m \tag{5}$$

za svaki $\Delta x$. Ova pretpostavka nije uvijek točno ispunjena, no u velikom broju slučajeva možemo smatrati da ona vrijedi. Na primjer, kod trapezne formule možemo uzeti 

$$\omega=-\frac{b-a}{12}\max_{x\in(a,b)}f''(x)$$

pa pretpostavka (5) znači da možemo uzeti (približno) isti $\omega$ za različite vrijednosti od $\Delta x$. Koristeći pretpostavku imamo

$$
\begin{aligned}
E_{n}&=I-I_n =\omega (\Delta x)^m =\omega \, \left(\frac{b-a}{n}\right)^m,\cr
E_{n/2}&=I-I_{n/2} =\omega (2\Delta x)^m =\omega\,  \left(\frac{b-a}{\frac{n}{2}}\right)^m.
\end{aligned}$$

Dakle,

$$
\displaystyle I_{n}-I_{n/2}=\omega \, (b-a)^m \left(\frac{1}{(n/2)^m}-\frac{1}{n^m}\right)$$

pa je

$$
\displaystyle \omega = \frac{n^m (n/2)^m}{(b-a)^m} \frac{I_{n}-I_{n/2}}{n^m-(n/2)^m},$$

što konačno daje

$$
\displaystyle E_{n}\approx \frac{(n/2)^m} {n^m-(n/2)^m} (I_{n}-I_{n/2})=E.$$

"""

# ╔═╡ e9054e4d-0e5a-4b26-a323-895ec36ddf73
function Trapez(f::Function,a::Number,b::Number,n::Int64=16)
    # n je broj intervala
    m=2
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx*(Y[1]/2+sum(Y[2:end-1])+Y[end]/2)
    # Richardsonova ekstrapolacija
    I₂=2*Δx*(Y[1]/2+sum(Y[3:2:end-2])+Y[end]/2)
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 0283137b-539e-43c2-b2e4-34f1890d78b9
function Simpson(f::Function,a::Number,b::Number,n::Int64)
    # n je broj intervala, djeljiv s 4
    m=4
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx/3*(Y[1]+4*sum(Y[2:2:end-1])+2*sum(Y[3:2:end-2])+Y[end])
    # Richardsonova ekstrapolacija
    I₂=2*Δx/3*(Y[1]+4*sum(Y[3:4:end-2])+2*sum(Y[5:4:end-4])+Y[end])
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 4baefc88-ba8b-43f7-8cec-48bf4d80c314
md"""
## Eliptički integral

Izračunajmo opseg elipse s polu-osima $2$ i $1$ (vidi  [Matematika 2, poglavlje 2.7.1](http://lavica.fesb.hr/mat2/predavanja/node41.html)). Parametarske jednadžbe elipse glase

$$
x=2\cos t,\quad y=\sin t,\quad t\in[0,\pi/2],$$

pa je četvrtina opsega jednaka

$$
\frac{O}{4}\int\limits_0^{\pi/2} \sqrt{(-2\sin t)^2+(\cos t)^2}\, dt
=2\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt.$$

Radi se o eliptičkom integralu druge vrste koji nije elementarno rješiv, ali se može naći u [tablicama](http://nvlpubs.nist.gov/nistpubs/jres/50/jresv50n1p43_A1b.pdf).

Vidimo da je 

$$O=8\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt \approx 8\cdot 1.21125.$$

"""

# ╔═╡ 612d35a0-558b-11eb-1609-551b746ec5b4
# 1/k iz tablice
1/(sqrt(3/4))

# ╔═╡ d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
begin
	f₁(x)=sqrt(1-(3.0)/4*cos(x)^2)
	Trapez(f₁,0,π/2,4)
end

# ╔═╡ c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
Trapez(f₁,0,pi/2,10)

# ╔═╡ e370021d-2f28-4ec2-85c6-6bb094d242d1
Trapez(f₁,0,pi/2,24)

# ╔═╡ faab641c-5a27-4e36-a4dd-c8e82272eae9
Simpson(f₁,0,π/2,4)

# ╔═╡ c58f0660-6dd6-4045-bb33-a74351c48f73
Simpson(f₁,0,π/2,16)

# ╔═╡ a7cb847f-bc7d-4299-9b13-37a698a1f28b
Simpson(f₁,0,π/2,24)

# ╔═╡ e5f42142-f569-4413-bacf-06ad8df859a7


# ╔═╡ aa07f01e-3bb9-4b47-b197-aba3c1310a2c
md"""
## Broj $\pi$

Vrijedi 

$$
\int_0^1 \frac{4}{1+x^2}\, dx=\pi.$$

Aproksimirajmo $\pi$ numeričkom integracijom i provjerimo pogrešku (vidi [Numerička matematika, poglavlje 7.3](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).

Pomoću trapezne formula možemo dobiti najviše pet točnih decimala. 
Simpsonova formula je točnija, ali je konvergencija spora.
"""

# ╔═╡ 44425308-ab3f-4793-b89e-32a18665bb80
plot(x->4/(1+x^2),0,1)

# ╔═╡ 311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# Za provjeru
BigFloat(π)

# ╔═╡ 1c13d010-afbf-45b4-8210-a7cfea8641a4
begin
	f₂(x)=4/(1+x^2)
	myπ=Trapez(f₂,0,1,10)
end

# ╔═╡ 94267bb3-624e-49b4-bbc8-2c24aacbd095
myπ[1]-BigFloat(π)

# ╔═╡ cb4df66a-1b00-403f-b49b-7a275dc2929f
myπ₁=Trapez(f₂,0,1,100)

# ╔═╡ 5cb0017b-1722-4f38-9f53-76ff788fce9e
myπ₁[1]-BigFloat(π)

# ╔═╡ d1b0fca4-b727-4670-ad25-3545dfb62629
myπ₂=Simpson(f₂,0,1,16)

# ╔═╡ 29e800e1-b38e-4dbf-a1f1-6084b35ef762
myπ₂[1]-BigFloat(π)

# ╔═╡ ad0d05f0-8e66-4410-956b-bf6802ada45a
myπ₃=Simpson(f₂,0,1,64)

# ╔═╡ 11c3d407-0923-4d65-85a2-b30911fd75d2
myπ₃[1]-BigFloat(π)

# ╔═╡ e7b2a060-7ae9-4418-a434-1693f73d04ac
md"""
# Gaussova kvadratura

Slično kao u formuli (1), integral aproksimiramo sumom umnožaka vrijednosti funkcije i odgovarajućih težina:

$$
\int_{a}^b \omega(x) f(x)\, dx=\sum_{k=1}^n \omega_k f(x_k),$$

gdje je $\omega(x)$ __težinska funkcija__.

Točke $x_k$ su nul-točke odgovarajućeg ortogonalnog polinoma $P_{n}(x)$ reda $n$, na primjer, __Legendreovog polinoma__ za $[a,b]=[-1,1]$ i $\omega(x)=1$ ili __Čebiševljevog polinoma__ za 
$[a,b]=[-1,1]$ i $\omega(x)=\displaystyle\frac{1}{\sqrt{1-x^2}}$ (vidi bilježnicu [NA12 Ortogonalni polinomi.ipynb](https://nbviewer.jupyter.org/github/ivanslapnicar/Numericka_analiza/blob/master/src/Jupyter/NA12%20Ortogonalni%20polinomi.ipynb)).

Težine su jednake

$$
\omega_k=\int_a^b \omega(x) \prod_{{i=1}\atop {i\neq k}}^n\frac{x-x_i}{x_k-x_i} \, dx.$$

__Pogreška__ je dana s

$$
E=\frac{f^{(2n)}(\xi)}{(2n)!}\int_a^b \omega(x) P_n^2(x)\, dx.$$

Za detalje vidi [Numerical Analysis, poglavlje 7.3](https://books.google.hr/books?id=kPDtAp3UZtIC&hl=hr&source=gbs_book_other_versions).

__Napomena__. Legendreovi i Čebiševljevi polinomi su definirani na intervalu $[-1,1]$ pa koristimo supstituciju

$$
\int_{a}^b \omega(x) f(x)\, dx = \frac{b-a}{2} \int_{-1}^1 \omega\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) f\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) dt.$$
"""

# ╔═╡ 656b359d-9504-4fa3-a327-d0977e10d642
mapnodes(x,a,b)=(b-a)*x/2 .+(a+b)/2

# ╔═╡ 427bfd71-5f3d-4b0b-b6f8-3753b5844a87
md"""
# Paketi

Profesionalne rutine za numeričku integraciju su složene, a većina programa ima ugrađene odgovarajuće naredbe ili funkcije. Tako, na primjer, 

* Ṁatlab ima naredbu `quad` koja koristi adaptivnu Simpsonovu formulu, a 
* Julia u paketu [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) ima funkciju `quadgk()` i računa integral u $O(n^2)$ operacija.
* Julia također ima i paket  [`FastGaussQuadrature.jl`](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) koji brzo računa točke i težine za zadani $n$ i razne težinske funkcije pa se pomoću točaka i težina lako izračuna integral u $O(n)$ operacija.
"""

# ╔═╡ cb4a8760-9ba7-4bb2-beea-9876a680829f
# ?quadgk

# ╔═╡ 4feab23e-b51b-4818-8343-a73deb67c942
# 1/8 opsega elipse, default je order=7
quadgk(f₁,0,π/2,order=10)

# ╔═╡ fa5a8850-c845-4328-94ca-cef66b4b52dd
# Broj π
quadgk(f₂,0,1,order=7)

# ╔═╡ ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# Granice mogu biti i beskonačne
quadgk(x->exp(-x),0,Inf)

# ╔═╡ 4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
varinfo(FastGaussQuadrature)

# ╔═╡ 5fbe0581-65b8-44b4-94f5-35950d579c92
# Na primjer
methods(gausschebyshev)

# ╔═╡ b9a58890-f3a2-49bb-836c-95a2971c9da7
gausschebyshev(16)

# ╔═╡ 4ae529de-6df8-4d38-917f-7b972904204d
# Sada računajmo integrale. U našem slučaju je ω(x)=1 
# pa nam treba Legendreov polinom.
ξ,ω=gausslegendre(64)

# ╔═╡ 988b6a6c-7732-47cd-915d-a26f44e1b3f0
# 1/8 opsega elipse, a=0, b=π/2
(π/2-0)/2*dot(ω,map(f₁,mapnodes(ξ,0,π/2)))

# ╔═╡ 5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# Broj π, a=0, b=1
(1-0)/2*dot(ω,map(f₂,mapnodes(ξ,0,1)))

# ╔═╡ 75fe9d8a-5e70-4d31-98f2-0ea82c511698
md"""
# Clenshaw-Curtisova kvadratura

Uz supstituciju $x=\cos\theta$, vrijedi

$$
I\equiv \int\limits_{-1}^1f(x)\, dx =\int\limits_0^\pi f(\cos\theta)\sin\theta \, d\theta.$$

Integral na desnoj strani se računa integriranjem Fourierovog reda parnog proširenja podintegralne funkcije:

$$
I\approx a_0+\sum_{k=1}^n \frac{2a_{2k}}{1-(2k)^2},$$

pri čemu se koeficijenti $a_k$ računaju formulom

$$
a_k=\frac{2}{\pi}\int\limits_0^\pi f(\cos\theta)\cos(k\theta)\,d\theta.$$

Koeficijenti $a_k$ se mogu računati numeričkom integracijom ili korištenjem brze Fourierove transformacije (FFT), što je puno brže.
Za detalje vidi [Clenshaw-Curtis Quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature).

Ukoliko se integrira na intervalu $[a,b]$, prebacivanje u interval $[-1,1]$ vrši se kao u formuli (3).
"""

# ╔═╡ 3d03b422-76c1-4f9b-ad69-12c07f2d4403
function ClenshawCurtis(f::Function,a::Number,b::Number,n::Int64)
    # Implementacija koristeći numeričku integraciju
    z=Vector{Float64}(undef,n)
    g(x)=f(mapnodes(x,a,b))
    for i=1:n
        h(x)=g(cos(x))*cos(2*(i-1)*x)
        z[i]=2*quadgk(h,0,pi)[1]/pi
    end
    return (z[1]+2*sum([z[i]/(1-4*(i-1)^2) for i=2:n]))*(b-a)/2
end

# ╔═╡ cc0b62ab-821c-4b16-9c73-ec131176e296
# 1/8 opsega elipse
ClenshawCurtis(f₁,0,pi/2,8)

# ╔═╡ 6ab103b9-e138-4c57-8730-83fe6313bbb0
# Broj π
ClenshawCurtis(f₂,0,1,8)

# ╔═╡ 35dbcd72-9f1d-43a2-9134-4116eefacc97
 π

# ╔═╡ 307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# "Nepravi" integral
ClenshawCurtis(x->exp(-x),0,1000,70)

# ╔═╡ 748ace02-d80f-4952-883b-de220083d58a
function ClenshawCurtisFFT(f::Function,a::Number,b::Number,n::Int64)
    # Brza implementacija pomoću fft(), 2^n je broj točaka
    g(x)=f(mapnodes(x,a,b))
    w=map(x->g(cos(x)),range(0,stop=2*pi,length=2^n))
    w[1]=(w[1]+w[end])/2
    z=real(fft(w))
    z/=2.0^(n-1)
    return (z[1]+2*sum([z[i]/(1-(i-1)^2) for i=3:2:2^(n-1)]))*(b-a)/2
end

# ╔═╡ 44406f97-1a7e-490b-9e75-3ab42e7f353f
ClenshawCurtisFFT(f₁,0,pi/2,4)

# ╔═╡ 41654b34-29c9-44d9-a0a2-d30a77d1b29c
ClenshawCurtisFFT(f₂,0,1,16),pi

# ╔═╡ 15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
ClenshawCurtisFFT(x->exp(-x),0,1000,18)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"

[compat]
FFTW = "~1.4.5"
FastGaussQuadrature = "~0.4.9"
Plots = "~1.25.4"
PlutoUI = "~0.7.27"
QuadGK = "~2.4.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9faf218ea18c51fcccaf956c8d39614c9d30fe8b"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

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
git-tree-sha1 = "d711603452231bad418bd5e0c91f1abd650cba71"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.3"

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

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

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

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "58d83dd5a78a36205bdfddb82b1bb67682e64487"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.9"

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
git-tree-sha1 = "b9a93bcdf34618031891ee56aad94cfff0843753"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f97acd98255568c3c9b416c5a3cf246c1315771b"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.0+0"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

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
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

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

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

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

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "68604313ed59f0408313228ba09e79252e4b2da8"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "71d65e9242935132e71c4fbf084451579491166a"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "fed057115644d04fba7f4d768faeeeff6ad11a60"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.27"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

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
git-tree-sha1 = "de9e88179b584ba9cf3cc5edbb7a41f26ce42cda"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

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
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

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
# ╠═7d2eceff-5dc8-4088-8504-8fa469948983
# ╠═3c1c274f-a299-4f97-a94f-5ac190e00f76
# ╟─c42b0321-1003-4935-857e-a8422fabd485
# ╟─5071e7f3-d353-4429-8885-ca9431edf5a9
# ╟─04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
# ╟─7a11ccd0-68c1-46d3-9249-3cc9f867ef65
# ╠═e9054e4d-0e5a-4b26-a323-895ec36ddf73
# ╠═0283137b-539e-43c2-b2e4-34f1890d78b9
# ╟─4baefc88-ba8b-43f7-8cec-48bf4d80c314
# ╠═612d35a0-558b-11eb-1609-551b746ec5b4
# ╠═d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
# ╠═c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
# ╠═e370021d-2f28-4ec2-85c6-6bb094d242d1
# ╠═faab641c-5a27-4e36-a4dd-c8e82272eae9
# ╠═c58f0660-6dd6-4045-bb33-a74351c48f73
# ╠═a7cb847f-bc7d-4299-9b13-37a698a1f28b
# ╠═e5f42142-f569-4413-bacf-06ad8df859a7
# ╟─aa07f01e-3bb9-4b47-b197-aba3c1310a2c
# ╠═44425308-ab3f-4793-b89e-32a18665bb80
# ╠═311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# ╠═1c13d010-afbf-45b4-8210-a7cfea8641a4
# ╠═94267bb3-624e-49b4-bbc8-2c24aacbd095
# ╠═cb4df66a-1b00-403f-b49b-7a275dc2929f
# ╠═5cb0017b-1722-4f38-9f53-76ff788fce9e
# ╠═d1b0fca4-b727-4670-ad25-3545dfb62629
# ╠═29e800e1-b38e-4dbf-a1f1-6084b35ef762
# ╠═ad0d05f0-8e66-4410-956b-bf6802ada45a
# ╠═11c3d407-0923-4d65-85a2-b30911fd75d2
# ╟─e7b2a060-7ae9-4418-a434-1693f73d04ac
# ╠═656b359d-9504-4fa3-a327-d0977e10d642
# ╟─427bfd71-5f3d-4b0b-b6f8-3753b5844a87
# ╠═cb4a8760-9ba7-4bb2-beea-9876a680829f
# ╠═4feab23e-b51b-4818-8343-a73deb67c942
# ╠═fa5a8850-c845-4328-94ca-cef66b4b52dd
# ╠═ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# ╠═4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
# ╠═5fbe0581-65b8-44b4-94f5-35950d579c92
# ╠═b9a58890-f3a2-49bb-836c-95a2971c9da7
# ╠═4ae529de-6df8-4d38-917f-7b972904204d
# ╠═988b6a6c-7732-47cd-915d-a26f44e1b3f0
# ╠═5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# ╟─75fe9d8a-5e70-4d31-98f2-0ea82c511698
# ╠═3d03b422-76c1-4f9b-ad69-12c07f2d4403
# ╠═cc0b62ab-821c-4b16-9c73-ec131176e296
# ╠═6ab103b9-e138-4c57-8730-83fe6313bbb0
# ╠═35dbcd72-9f1d-43a2-9134-4116eefacc97
# ╠═307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# ╠═456ae14f-d92c-4561-bcfd-d1b9849b8244
# ╠═748ace02-d80f-4952-883b-de220083d58a
# ╠═44406f97-1a7e-490b-9e75-3ab42e7f353f
# ╠═41654b34-29c9-44d9-a0a2-d30a77d1b29c
# ╠═15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
